#!/usr/bin/env python3
"""Fan-Skill Rule Engine — C-layer condition matching."""
import yaml, sys

def match_condition(condition, profile):
    if condition == "none_matched":
        return False
    for key, val in condition.items():
        if key not in profile:
            return False
        actual = str(profile[key])
        if isinstance(val, list):
            if actual not in [str(v) for v in val]:
                return False
        elif isinstance(val, str) and (val.startswith('>=') or val.startswith('<=') or val.startswith('>') or val.startswith('<')):
            try:
                op, threshold = val[0:2] if val[1] == '=' else (val[0], val[1:])
                if op in ('>=', '=>'):
                    threshold = val[2:] if val.startswith('>=') else val[1:]
                actual_f = float(actual)
                threshold_f = float(threshold)
                if val.startswith('>='):
                    if not (actual_f >= threshold_f): return False
                elif val.startswith('<='):
                    if not (actual_f <= threshold_f): return False
                elif val.startswith('>'):
                    if not (actual_f > threshold_f): return False
                elif val.startswith('<'):
                    if not (actual_f < threshold_f): return False
            except ValueError:
                return False
        else:
            if str(val) != actual:
                return False
    return True

def evaluate_rules(rules_section, profile):
    best_match = None
    best_priority = -1
    for rule in rules_section:
        condition = rule.get('condition', {})
        priority = rule.get('priority', 0)
        if match_condition(condition, profile) and priority > best_priority:
            best_match = rule
            best_priority = priority
    return best_match

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: rule_engine.py <rules.yaml> <profile.yaml>")
        sys.exit(1)
    with open(sys.argv[1]) as f:
        rules = yaml.safe_load(f)
    with open(sys.argv[2]) as f:
        profile = yaml.safe_load(f)

    for section_name, section_rules in rules.items():
        if section_name in ('name', 'description', 'triggers', 'inputs', 'outputs', 'version', 'analysis'):
            continue
        if isinstance(section_rules, list):
            result = evaluate_rules(section_rules, profile)
            if result:
                action = result.get('recommend', '') or result.get('action', '') or str(result.get('interpretation', ''))
                print(f"[{section_name}] rule_id={result.get('rule_id', '?')} → {action}")
