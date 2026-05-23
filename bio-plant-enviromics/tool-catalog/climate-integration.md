# Climate Data Integration (气候数据整合)

**Goal:** Acquire, process, and integrate gridded climate/weather data with
plant breeding trial metadata

## Prerequisites
- R 4.0+, packages: `nasapower`, `geodata`, `raster`, `sf`, `tidyverse`
- Command-line: `cdo` (Climate Data Operators), GDAL
- Trial metadata: location (lat/lon), planting date, harvest date

## Data Sources

| Source | Resolution | Variables | R Package |
|--------|-----------|-----------|-----------|
| NASA POWER | 0.5° (~55km) | Tmin, Tmax, Precip, Rad | `nasapower` |
| WorldClim | 1 km² (30yr avg) | Bioclimatic vars | `geodata` |
| CHIRPS | 0.05° (~5km) | Precipitation | `chirps` |
| SoilGrids | 250m | Soil texture, pH, C, N | `geodata` |
| CRU TS | 0.5° | Monthly climate | `utils::download.file` |
| AgERA5 | 0.1° (~11km) | Agro-climate vars | Copernicus CDS API |
| Daymet (N. America) | 1km | Daily Tmin, Tmax, Precip | `daymetr` |

## Basic Usage

### NASA POWER: Quick Daily Weather Retrieval

```r
library(nasapower)

# Define trial locations and dates
locations <- data.frame(
  site = c("BJ_2023", "HN_2023", "HB_2023"),
  lat = c(39.9, 18.3, 30.6),
  lon = c(116.4, 109.5, 114.3),
  planting = as.Date(c("2023-04-15", "2023-03-01", "2023-04-01")),
  harvest = as.Date(c("2023-06-15", "2023-05-20", "2023-06-10"))
)

# Retrieve daily weather for each location
fetch_weather <- function(lat, lon, start_date, end_date) {
  get_power(
    community = "ag",
    lonlat = c(lon, lat),
    pars = c("T2M_MAX", "T2M_MIN", "PRECTOTCORR", "ALLSKY_SFC_SW_DWN"),
    dates = c(start_date, end_date),
    temporal_api = "daily"
  )
}

# Fetch for all locations
weather_list <- locations %>%
  rowwise() %>%
  mutate(weather = list(fetch_weather(lat, lon, planting, harvest)))
```

### WorldClim: 30-year Climate Normals

```r
library(geodata)

# Download bioclimatic variables at 2.5 arc-minute resolution
bioclim <- worldclim_global(
  var = "bio",
  res = 2.5,
  path = tempdir()
)

# Extract at trial locations
locations_sf <- st_as_sf(locations, coords = c("lon", "lat"), crs = 4326)
bioclim_values <- extract(bioclim, locations_sf)
# BIO1 = Annual Mean Temperature
# BIO12 = Annual Precipitation
# BIO5 = Max Temperature of Warmest Month
# BIO15 = Precipitation Seasonality
```

### Computing Environmental Covariates (W-series)

```r
# After fetching daily weather, compute environmental covariates
# Based on the concept from EnvRtype / W-series

compute_env_covariates <- function(daily_weather, planting_date, harvest_date) {
  daily_weather %>%
    mutate(
      date = as.Date(YYYYMMDD, format = "%Y%m%d"),
      GDD = pmax((T2M_MAX + T2M_MIN) / 2 - 10, 0),  # base temp = 10°C
      # Accumulate over growing season
      GDD_cum = cumsum(GDD),
      # Window-specific variables
      days_from_planting = as.numeric(date - planting_date),
      phase = case_when(
        days_from_planting <= 30 ~ "vegetative",
        days_from_planting <= 60 ~ "reproductive",
        TRUE ~ "grain_filling"
      )
    ) %>%
    group_by(phase) %>%
    summarise(
      Tmean = mean((T2M_MAX + T2M_MIN) / 2, na.rm = TRUE),
      Tmax_mean = mean(T2M_MAX, na.rm = TRUE),
      Tmin_mean = mean(T2M_MIN, na.rm = TRUE),
      Precip_sum = sum(PRECTOTCORR, na.rm = TRUE),
      Rad_mean = mean(ALLSKY_SFC_SW_DWN, na.rm = TRUE),
      GDD_sum = sum(GDD, na.rm = TRUE)
    )
}
```

### SoilGrids: Soil Properties

```r
library(geodata)

# Download soil data
soil_ph <- soil_world(
  var = "phh2o",    # pH in water
  depth = 5,        # 0-5 cm
  path = tempdir()
)

# Extract at locations
soil_values <- extract(soil_ph, locations_sf)
```

## Using CDO for NetCDF/Gridded Climate Data

```bash
# CDO (Climate Data Operators) for large gridded datasets

# Extract grid points for trial locations from NetCDF
cdo -remapnn,lon=116.4_lat=39.9 input.nc output_BJ.nc

# Compute growing season aggregates
cdo -yearsmon -select,startdate=2023-04-15,enddate=2023-06-15 \
    input.nc season_subset.nc

# Calculate growing degree days
cdo -yearsum -expr,'GDD=max((Tmax+Tmin)/2-10,0)' input.nc gdd.nc
```

## Plant Relevance

- **Crop-specific base temperatures for GDD**:
  - Wheat, barley: 0°C or 5°C
  - Maize, sorghum: 10°C
  - Rice: 10°C
  - Soybean: 10°C
  - Cotton: 15.6°C
- **Critical windows vary by species**: Always consult crop-specific
  phenology before computing environmental indices. The grain-filling
  period may be 20 days (wheat) or 45 days (maize).
- **Elevation matters**: Mountainous regions have large climate variation
  over short distances. 1km gridded data may be insufficient. Use
  CHELSA (30 arc-second) or downscaled products.
- **Data quality**: NASA POWER is satellite-derived and may have biases
  in cloudy regions or complex terrain. If ground station data is
  available, prefer it or at least cross-validate.

## Common Errors

| Error | Cause | Fix |
|-------|-------|-----|
| Wrong time period extracted | Misalignment of planting date format | Verify date parsing; use consistent format |
| Missing values at coastal locations | Grid cell falls in ocean | Buffer inland or use nearest non-NA cell |
| Very low precipitation values | Units mismatch (mm vs kg/m²) | Check dataset documentation for units |
| High GDD but low yield | Heat stress not captured by mean GDD | Add Tmax > threshold (e.g., Tmax > 35°C) as stress indicator |
