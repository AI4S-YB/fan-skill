# Contributing to Fan-Skill

## Adding a New Analysis

1. Copy `templates/rules-template.yaml` → `knowledge-base/<your-analysis>/rules.yaml`
2. Copy `templates/notebook-template.md` → `knowledge-base/<your-analysis>/notebook.md`
3. Fill in the templates with your domain expertise
4. Add tool docs to `tool-registry/` if using new tools
5. Run `engine/validate_entry.sh knowledge-base/<your-analysis>/` — must pass
6. Submit a PR

## Quality Standards

- [ ] rules.yaml has name, description, triggers, inputs, outputs
- [ ] Every decision rule has rule_id, priority, condition, recommend/action, reason
- [ ] At least one fallback rule (priority 0, action: delegate_to_expert)
- [ ] notebook.md covers method selection, common pitfalls, plant-specific notes
- [ ] All tool_id values reference existing tool-registry entries or new ones you've added
- [ ] Validation passes: `engine/validate_entry.sh` returns 0

## B+C Architecture

Fan-Skill uses a dual-mode architecture:
- **C-layer (rules.yaml)**: Structured, deterministic, testable decision rules
- **B-layer (notebook.md)**: Narrative expert guidance for flexible reasoning

New entries should follow this pattern. See `knowledge-base/gwas/` for a reference implementation.
