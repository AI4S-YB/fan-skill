#!/bin/bash
# test_decision_rules.sh — L2: 验证每条规则可被触发
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SKILL_DIR="$(dirname "$SCRIPT_DIR")"

echo "=== L2: Decision Rule Reachability Test ==="

python3 << 'PYEOF'
import yaml, sys, os

skill_dir = os.environ.get('SKILL_DIR', '.')
with open(f'{skill_dir}/test/test_cases.yaml') as f:
    cases = yaml.safe_load(f)['test_cases']

with open(f'{skill_dir}/decision-matrix.yaml') as f:
    matrix = yaml.safe_load(f)

passed = 0
failed = 0

def match_condition(condition, profile):
    """检查 profile 是否满足 condition"""
    if condition == "none_matched":
        return False
    for key, val in condition.items():
        if key not in profile:
            return False
        actual = str(profile[key])
        if isinstance(val, list):
            # list value means "actual must be one of these"
            if actual not in [str(v) for v in val]:
                return False
        elif isinstance(val, str) and (val.startswith('>=') or val.startswith('<=') or val.startswith('>') or val.startswith('<')):
            # numeric comparison
            try:
                if val.startswith('>='):
                    if not (float(actual) >= float(val[2:])):
                        return False
                elif val.startswith('<='):
                    if not (float(actual) <= float(val[2:])):
                        return False
                elif val.startswith('>'):
                    if not (float(actual) > float(val[1:])):
                        return False
                elif val.startswith('<'):
                    if not (float(actual) < float(val[1:])):
                        return False
            except ValueError:
                return False
        else:
            if str(val) != actual:
                return False
    return True

for case in cases:
    section = case['rule_section']
    rules = matrix.get(section, [])
    matched = None
    best_priority = -1

    for rule in rules:
        condition = rule.get('condition', {})
        priority = rule.get('priority', 0)
        if match_condition(condition, case['input_profile']) and priority > best_priority:
            matched = rule
            best_priority = priority

    # Fallback rule lookup: check the fallback section when no rule matched
    if matched is None:
        fallback_rules = matrix.get('fallback', [])
        for rule in fallback_rules:
            if rule.get('action') == 'delegate_to_expert':
                matched = rule
                break

    expected = case['expected']
    if matched is None and expected.get('action') == 'delegate_to_expert':
        print(f"  [PASS] {case['name']} → fallback (expected)")
        passed += 1
    elif matched is None:
        print(f"  [FAIL] {case['name']} → no rule matched (expected {expected.get('rule_id', '?')})")
        failed += 1
    else:
        match_ok = True
        for key, val in expected.items():
            if str(matched.get(key, '')) != str(val):
                print(f"  [FAIL] {case['name']} → {key}: got '{matched.get(key)}', expected '{val}'")
                match_ok = False
                failed += 1
                break
        if match_ok:
            print(f"  [PASS] {case['name']} → {matched.get('rule_id', '?')}")
            passed += 1

print(f"\n=== Results: {passed} passed, {failed} failed ===")
sys.exit(0 if failed == 0 else 1)
PYEOF
