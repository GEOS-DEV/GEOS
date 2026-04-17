This directory contains runtime CI configuration consumed by `.github/workflows/ci_tests.yml`.

Organization-specific config files live in `.github/ci/orgs/<org>.json`, where `<org>` matches
`github.repository_owner`. The resolver also checks the lowercase form of the owner name.
Each file contains the config entry for that org directly; it does not need a top-level wrapper keyed
by the org name. These files are the source of truth for org-specific CI behavior.

Each org JSON file must define:

- `storage_provider`
- `sccache_profile`
- `integrated_tests_artifact_bucket_path`
- `runners`

Current runner keys used by the workflow:

- `default`
- `cpu_heavy`
- `integrated_tests`
- `code_coverage`
- `cuda`

Selection order:

1. Checked-in org config file in `.github/ci/orgs/`
2. Optional per-field overrides from GitHub Actions vars such as `CI_STORAGE_PROVIDER`

Provider behavior still comes from built-in defaults plus `vars.CI_PROVIDER_CONFIG_JSON`.
Secrets still come from GitHub org/repo secrets.
