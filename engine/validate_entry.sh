#!/bin/bash
# validate_entry.sh — verify a knowledge entry meets B+C standards
set -euo pipefail

ENTRY_DIR="${1:-.}"
ERRORS=0

echo "=== Validating: $ENTRY_DIR ==="

# Check required files
for f in rules.yaml notebook.md; do
    if [ -f "$ENTRY_DIR/$f" ]; then
        echo "  [OK] $f"
    else
        echo "  [MISSING] $f"
        ERRORS=$((ERRORS + 1))
    fi
done

# Validate rules.yaml structure
if [ -f "$ENTRY_DIR/rules.yaml" ]; then
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
decision_sections = [k for k in rules if k not in required + ['version', 'analysis', 'inputs_required']]
if not decision_sections:
    print('  [MISSING] rules.yaml: no decision sections found')
    sys.exit(1)
print(f'  [OK] rules.yaml: {len(decision_sections)} decision sections')

# Check every rule has rule_id or id
rule_count = 0
for section in decision_sections:
    if isinstance(rules[section], list):
        for rule in rules[section]:
            if 'rule_id' not in rule and 'id' not in rule:
                print(f'  [MISSING] rule_id or id in section {section}')
                sys.exit(1)
            rule_count += 1
print(f'  [OK] rules.yaml: {rule_count} rules with rule_id or id')
"
    [ $? -ne 0 ] && ERRORS=$((ERRORS + 1))
fi

echo "=== $ERRORS errors ==="
exit $ERRORS
