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

Optional org JSON fields:

- `runner_ca_bundle_host_paths`
  A map keyed by runner label, not by workflow role. Each value is a host-side certificate bundle path
  to bind into the container as `/certs/ca-bundle.crt`.
  Example:
  `{"streak2": "/etc/pki/tls/certs/ca-bundle.crt"}`
  The resolver first looks for an exact runner-label match. If none is found, it falls back to the
  prefix before the first `-`, so a runner label like `streak2-32core` will reuse the `streak2` entry.
  Use this only for runners that need an LLNL-style outbound cert bundle; leave it empty or omit it for
  machines that do not.

Current runner keys used by the workflow:

- `default`
- `cpu_heavy`
- `integrated_tests`
- `code_coverage`
- `cuda`

Each runner value must be a single runner label string, for example `"ubuntu-22.04"` or `"streak2"`.

Selection order:

1. Checked-in org config file in `.github/ci/orgs/`
2. Optional per-field overrides from GitHub Actions vars such as `CI_STORAGE_PROVIDER`

Provider behavior still comes from built-in defaults plus `vars.CI_PROVIDER_CONFIG_JSON`.
Secrets still come from GitHub org/repo secrets.
