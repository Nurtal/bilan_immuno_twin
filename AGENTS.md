# Agent skills

### Issue tracker

Issues and specs live as GitHub issues in `Nurtal/bilan_immuno_twin`, via the `gh` CLI. See `docs/agents/issue-tracker.md`.

### Triage labels

The five canonical triage labels (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`) as default names. See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: root `CONTEXT.md` + `docs/adr/`. See `docs/agents/domain.md`.

### CLI output

The `bilan` CLI (ADR-0004) defaults to JSON. `--format csv` emits one or more
flat tables, each headed, separated by a **blank line**, so downstream tooling
must split on blank lines before parsing. A subcommand with several outputs
emits them as separate tables (e.g. `simulate --perturb` emits the baseline and
perturbed trajectory tables, then comparison and score tables). When a command
parses a multi-table stream it should treat a blank line as the table
separator. See `docs/spec.md` (US #21/#22).
