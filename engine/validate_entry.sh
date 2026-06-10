#!/bin/bash
# validate_entry.sh — verify a knowledge entry meets B+C standards
# Enhanced version with params coverage, plant-specific params, and tool_id checks
# Usage:
#   validate_entry.sh <entry_dir>                  # full 4-module validation (general layer)
#   validate_entry.sh --layer user <entry_dir>     # minimal validation (user layer)
set -euo pipefail

LAYER="general"
ENTRY_DIR=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --layer)
            if [ -z "${2:-}" ] || [[ "$2" =~ ^-- ]]; then
                echo "ERROR: --layer requires a value (general or user)"
                exit 2
            fi
            LAYER="$2"
            shift 2
            ;;
        *)
            ENTRY_DIR="$1"
            shift
            ;;
    esac
done

ENTRY_DIR="${ENTRY_DIR:-.}"

# Validate LAYER value
if [ "$LAYER" != "general" ] && [ "$LAYER" != "user" ]; then
    echo "ERROR: --layer must be \"general\" or \"user\", got \"$LAYER\""
    exit 2
fi
ERRORS=0
WARNINGS=0

# Configurable thresholds
PARAMS_COVERAGE_THRESHOLD=${PARAMS_COVERAGE_THRESHOLD:-60}

echo "=== Validating: $ENTRY_DIR ==="

# Check required files — different requirements per layer
if [ "$LAYER" = "user" ]; then
    # User layer: meta.yaml + notebook.md required; rules.yaml/consult-guide.md/analysis-primer.md optional
    for f in meta.yaml notebook.md; do
        if [ -f "$ENTRY_DIR/$f" ]; then
            echo "  [OK] $f"
        else
            echo "  [MISSING] $f (required for user layer)"
            ERRORS=$((ERRORS + 1))
        fi
    done
    for f in rules.yaml consult-guide.md analysis-primer.md; do
        if [ -f "$ENTRY_DIR/$f" ]; then
            echo "  [OK] $f (optional)"
        else
            echo "  [INFO] $f not present (optional for user layer)"
        fi
    done

    # Validate meta.yaml structure (replaces rules.yaml validation for user layer)
    if [ -f "$ENTRY_DIR/meta.yaml" ]; then
        python3 -c "
import yaml, sys

with open('$ENTRY_DIR/meta.yaml') as f:
    meta = yaml.safe_load(f)

required = ['name', 'triggers', 'inputs', 'outputs']
for key in required:
    if key not in meta or not meta[key]:
        print(f'  [MISSING] meta.yaml: {key} is empty or missing')
        sys.exit(1)
    else:
        val = meta[key]
        if isinstance(val, list):
            print(f'  [OK] meta.yaml: {key} ({len(val)} items)')
        else:
            print(f'  [OK] meta.yaml: {key} = {val}')

# Check extends field
if 'extends' in meta:
    print(f'  [OK] meta.yaml: extends = {meta[\"extends\"]}')
else:
    print('  [INFO] meta.yaml: no extends field (standalone entry)')

# Check status
if 'status' in meta:
    status = meta['status']
    if status in ('draft', 'confirmed', 'matured'):
        print(f'  [OK] meta.yaml: status = {status}')
    else:
        print(f'  [WARN] meta.yaml: unknown status \"{status}\", expected draft/confirmed/matured')
else:
    print('  [INFO] meta.yaml: no status field')

print(f'VALIDATION_PASSED=true')
"
        [ $? -ne 0 ] && ERRORS=$((ERRORS + 1))
    fi
else
    # General layer: all 4 modules required (original behavior)
    for f in rules.yaml notebook.md consult-guide.md analysis-primer.md; do
        if [ -f "$ENTRY_DIR/$f" ]; then
            echo "  [OK] $f"
        else
            echo "  [MISSING] $f"
            ERRORS=$((ERRORS + 1))
        fi
    done
fi

# Check minimum line counts for new module files
if [ -f "$ENTRY_DIR/consult-guide.md" ]; then
    LINES=$(wc -l < "$ENTRY_DIR/consult-guide.md" | tr -d ' ')
    if [ "$LAYER" = "user" ]; then
        echo "  [OK] consult-guide.md: $LINES lines (optional for user layer)"
    elif [ "$LINES" -ge 30 ]; then
        echo "  [OK] consult-guide.md: $LINES lines (>= 30)"
    else
        echo "  [FAIL] consult-guide.md: $LINES lines (< 30 required)"
        ERRORS=$((ERRORS + 1))
    fi
fi

if [ -f "$ENTRY_DIR/analysis-primer.md" ]; then
    LINES=$(wc -l < "$ENTRY_DIR/analysis-primer.md" | tr -d ' ')
    if [ "$LAYER" = "user" ]; then
        echo "  [OK] analysis-primer.md: $LINES lines (optional for user layer)"
    elif [ "$LINES" -ge 20 ]; then
        echo "  [OK] analysis-primer.md: $LINES lines (>= 20)"
    else
        echo "  [FAIL] analysis-primer.md: $LINES lines (< 20 required)"
        ERRORS=$((ERRORS + 1))
    fi
fi

# Validate rules.yaml structure (required for general layer, optional for user layer)
if [ -f "$ENTRY_DIR/rules.yaml" ]; then
    if [ "$LAYER" = "user" ]; then
        echo "  [INFO] rules.yaml present (optional for user layer — performing deep check)"
    fi
    python3 -c "
import yaml, sys

with open('$ENTRY_DIR/rules.yaml') as f:
    rules = yaml.safe_load(f)

required = ['name', 'description', 'triggers', 'inputs', 'outputs']
for key in required:
    if key not in rules:
        print(f'  [MISSING] rules.yaml: {key}')
        sys.exit(1)
    else:
        print(f'  [OK] rules.yaml: {key}')

# Check at least one decision section exists
decision_sections = [k for k in rules if k not in required + ['version', 'analysis', 'inputs_required', 'design_gates']]
if not decision_sections:
    print('  [MISSING] rules.yaml: no decision sections found')
    sys.exit(1)
print(f'  [OK] rules.yaml: {len(decision_sections)} decision sections')

# ============================================================
# NEW: Check rule_id and params coverage
# ============================================================
rule_count = 0
rules_with_params = 0
tool_related_rules = 0
tool_related_with_tool_id = 0
plant_specific_params = 0
missing_rule_ids = []

# Sections that typically involve tool calls (need tool_id)
tool_sections = ['method', 'tool', 'step', 'pipeline', 'enrichment', 'align',
                 'quant', 'assembly', 'call', 'detect', 'annotation', 'annot']
# Sections that are logic/control (don't need tool_id)
logic_sections = ['scenario', 'gate', 'decision', 'quality', 'format',
                  'wgcna', 'multiple', 'filter', 'threshold']

# Plant-specific parameter keywords
plant_keywords = [
    'plant', 'arabidopsis', 'rice', 'maize', 'wheat', 'soybean', 'tomato',
    'polyploid', 'subgenome', 'homeolog', 'tss_region', 'promoter_length',
    'gtr+g+asc', 'asc_corr', 'lfc_threshold', 'gap_penalty'
]

for section in decision_sections:
    if isinstance(rules[section], list):
        for rule in rules[section]:
            # Check rule_id
            if 'rule_id' not in rule and 'id' not in rule:
                missing_rule_ids.append(section)
                continue
            rule_count += 1

            # Check params
            if 'params' in rule and rule['params']:
                rules_with_params += 1
                # Check for plant-specific params
                params_str = str(rule['params']).lower()
                for kw in plant_keywords:
                    if kw in params_str:
                        plant_specific_params += 1
                        break

            # Check tool_id - distinguish tool-related vs logic rules
            section_lower = section.lower()
            is_tool_section = any(t in section_lower for t in tool_sections)
            is_logic_section = any(t in section_lower for t in logic_sections)

            # Also check if this is a fallback/delegate rule (doesn't need tool_id)
            rule_id_val = rule.get('rule_id') or rule.get('id', '')
            is_fallback = 'fallback' in str(rule_id_val).lower() or 'delegate' in str(rule.get('action', '')).lower()

            if is_tool_section and not is_logic_section and not is_fallback:
                tool_related_rules += 1
                if 'tool_id' in rule and rule['tool_id']:
                    tool_related_with_tool_id += 1

# Report rule_id check
if missing_rule_ids:
    print(f'  [FAIL] rules.yaml: missing rule_id in sections: {missing_rule_ids}')
    sys.exit(1)
print(f'  [OK] rules.yaml: {rule_count} rules with rule_id or id')

# ============================================================
# NEW: Report params coverage
# ============================================================
if rule_count > 0:
    coverage = (rules_with_params / rule_count) * 100
    if coverage >= $PARAMS_COVERAGE_THRESHOLD:
        print(f'  [OK] params coverage: {rules_with_params}/{rule_count} ({coverage:.1f}%)')
    else:
        print(f'  [WARN] params coverage: {rules_with_params}/{rule_count} ({coverage:.1f}%) < $PARAMS_COVERAGE_THRESHOLD% threshold')

    # Report plant-specific params
    if plant_specific_params > 0:
        print(f'  [OK] plant-specific params: {plant_specific_params} rules')
    else:
        print(f'  [INFO] no plant-specific params detected (consider adding plant defaults)')

    # Report tool_id coverage (only for tool-related rules)
    if tool_related_rules > 0:
        tool_coverage = (tool_related_with_tool_id / tool_related_rules) * 100
        if tool_coverage >= 80:
            print(f'  [OK] tool_id coverage: {tool_related_with_tool_id}/{tool_related_rules} ({tool_coverage:.1f}%) [tool-related rules]')
        elif tool_coverage >= 50:
            print(f'  [WARN] tool_id coverage: {tool_related_with_tool_id}/{tool_related_rules} ({tool_coverage:.1f}%) [tool-related rules, target: 80%]')
        else:
            print(f'  [FAIL] tool_id coverage: {tool_related_with_tool_id}/{tool_related_rules} ({tool_coverage:.1f}%) [tool-related rules, target: 80%]')
        print(f'  [INFO] logic/control rules: {rule_count - tool_related_rules} (tool_id not required)')
    else:
        print(f'  [INFO] no tool-related rules found')

# ============================================================
# NEW: Check design_gates
# ============================================================
if 'design_gates' in rules:
    gates = rules['design_gates']
    if isinstance(gates, list) and len(gates) > 0:
        print(f'  [OK] design_gates: {len(gates)} gates defined')
    else:
        print('  [WARN] design_gates: empty or invalid')
else:
    print('  [INFO] no design_gates section (recommended for quality control)')

print(f'VALIDATION_PASSED=true')
"
    [ $? -ne 0 ] && ERRORS=$((ERRORS + 1))
fi

# ============================================================
# NEW: Check tool-registry references (if applicable)
# ============================================================
TOOL_REGISTRY="tool-registry"
if [ -d "$TOOL_REGISTRY" ]; then
    # Extract tool_ids from rules.yaml and check if corresponding docs exist
    if [ -f "$ENTRY_DIR/rules.yaml" ]; then
        python3 -c "
import yaml, os

with open('$ENTRY_DIR/rules.yaml') as f:
    rules = yaml.safe_load(f)

tool_ids = set()
for key, value in rules.items():
    if isinstance(value, list):
        for item in value:
            if isinstance(item, dict) and 'tool_id' in item:
                tool_ids.add(item['tool_id'])

if tool_ids:
    missing_tools = []
    for tid in tool_ids:
        doc_file = f'tool-registry/{tid}.md'
        if not os.path.exists(doc_file):
            missing_tools.append(tid)

    if missing_tools:
        print(f'  [INFO] tool_id without registry doc: {missing_tools[:5]}')  # Show first 5
    else:
        print(f'  [OK] all tool_ids have registry docs ({len(tool_ids)} checked)')
"
    fi
fi

echo ""
echo "=== Summary ==="
if [ $ERRORS -gt 0 ]; then
    echo "  ERRORS: $ERRORS"
    echo "  FAILED"
    exit 1
else
    echo "  PASSED"
    exit 0
fi
