#!/usr/bin/env python3
"""
[DEPRECATED] Chain discovery via Python keyword matching.
Superseded by LLM-driven semantic matching in SKILL.md Phase 2 (fix #1, 2026-05-23).

The Agent now reads knowledge-base entries directly and reasons about
chain feasibility based on inputs/outputs declarations. This script
is kept for reference and may be revived if hard constraint checking
is needed in the future.
"""
"""
Fan-Skill Chain Discovery Engine
Given user goal + available data, find all feasible analysis paths.
"""
import yaml, os, sys, glob

def load_knowledge_base(kb_dir="knowledge-base"):
    entries = {}
    for rules_file in sorted(glob.glob(f"{kb_dir}/*/rules.yaml")):
        name = os.path.basename(os.path.dirname(rules_file))
        with open(rules_file) as f:
            entries[name] = yaml.safe_load(f)
    return entries

def match_entries(user_goal, entries):
    """Match user intent to knowledge entries via triggers and description."""
    relevant = {}
    goal_lower = user_goal.lower()
    for name, entry in entries.items():
        triggers = entry.get("triggers", [])
        desc = entry.get("description", "").lower()
        score = 0
        for t in triggers:
            if t.lower() in goal_lower:
                score += 10
        if any(word in goal_lower for word in desc.split()[:10]):
            score += 3
        if score > 0:
            relevant[name] = {"entry": entry, "score": score}
    return dict(sorted(relevant.items(), key=lambda x: x[1]["score"], reverse=True))

def check_inputs(required_inputs, user_data):
    """Check which required inputs are available in user data."""
    missing = []
    for inp in required_inputs:
        found = False
        for key in user_data:
            if inp.lower() in key.lower() or key.lower() in inp.lower():
                found = True
                break
        if not found:
            missing.append(inp)
    return missing

def discover_chains(user_goal, user_data, entries):
    """Discover feasible analysis chains."""
    relevant = match_entries(user_goal, entries)

    ready = []
    data_needed = []

    for name, info in relevant.items():
        entry = info["entry"]
        required = entry.get("inputs", [])
        missing = check_inputs(required, user_data)
        outputs = entry.get("outputs", [])

        if not missing:
            ready.append({
                "name": name,
                "description": entry.get("description", ""),
                "score": info["score"],
                "outputs": outputs,
            })
        else:
            data_needed.append({
                "name": name,
                "description": entry.get("description", ""),
                "score": info["score"],
                "missing_data": missing,
                "outputs": outputs,
            })

    # Discover chains: if one entry's outputs satisfy another's inputs
    chains = []
    for r in ready:
        chain = [r["name"]]
        extended = user_data.copy()
        extended.update({o: True for o in r["outputs"]})
        for n, info in relevant.items():
            if n == r["name"]:
                continue
            if n not in [c for c in chain]:
                if not check_inputs(info["entry"].get("inputs", []), extended):
                    chains.append(chain + [n])

    return {
        "ready": ready,
        "data_needed": data_needed,
        "chains": chains,
    }

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: discover_chains.py '<user_goal>' [user_data_json]")
        sys.exit(1)

    goal = sys.argv[1]
    user_data = {}
    if len(sys.argv) > 2:
        user_data = yaml.safe_load(sys.argv[2]) if isinstance(sys.argv[2], str) else {}

    entries = load_knowledge_base()
    result = discover_chains(goal, user_data, entries)
    print(yaml.dump(result, allow_unicode=True, default_flow_style=False))
