This directory contains runtime CI configuration consumed by `.github/workflows/ci_tests.yml`.

Two checked-in files drive behavior for each CI run:

1. `.github/ci/orgs/<org>.json` — per-organization choices (which provider, which runners, etc.).
2. `.github/ci/providers/<storage_provider>.json` — the shell commands and URI scheme for that
   provider. Reused by any org that selects this provider.

Every field in both files is required. There are no defaults: the full behavior of a CI run is
always visible by reading these two JSON files and nothing else.

## Org files (`.github/ci/orgs/<org>.json`)

`<org>` matches `github.repository_owner`. The resolver also checks the lowercase form of the
owner name. Each file contains the config entry for that org directly; there is no top-level
wrapper keyed by the org name.

Required fields (all of them, always):

- `storage_provider`
  Name of a file under `.github/ci/providers/`. Must match `^[A-Za-z0-9_-]+$`. Selects the upload
  CLI and URI scheme.

- `sccache_profile`
  Name of the sccache profile to pull from the inherited `SCCACHE_PROFILES_JSON` secret. Does not
  need to match `storage_provider`, though it typically does.

- `integrated_tests_artifact_bucket_path`
  Path under the bucket where integrated-test artifacts are uploaded. For providers whose URI
  scheme is `<scheme>://<bucket>/<path>`, the first path segment is the bucket name.

- `integrated_tests_baseline_fallback_public_url_prefix`
  Public URL prefix used by integrated tests when the requested baseline archive is missing from
  both the runner-local baseline cache and this org's public artifact URL prefix. Use an empty
  string to disable fallback. The workflow appends the baseline archive filename from
  `.integrated_tests.yaml`, e.g. `baseline_integratedTests-pr3959-16478-faf1698.tar.gz`.

- `artifact_public_url_base`
  Public URL root for uploaded artifacts. The consumer builds each final URL as
  `${artifact_public_url_base}/${integrated_tests_artifact_bucket_path}/${filename}`, so the base
  must already include any bucket/host prefix required for a valid public URL. For GCS this is
  typically `https://storage.googleapis.com`. For Cloudflare R2 public dev URLs it is the
  `https://pub-<hash>.r2.dev` domain for the bucket.

- `artifact_public_url_bucket_scoped`
  Boolean. Set `true` when `artifact_public_url_base` is rooted *inside* a specific bucket and
  cannot traverse up to the account root (e.g. Cloudflare R2 public dev URLs of the form
  `https://pub-<hash>.r2.dev`, which always resolve to one bucket). With the flag set, the
  consumer strips the first path segment (the bucket name) from
  `integrated_tests_artifact_bucket_path` before appending it to the base, so a base of
  `https://pub-<hash>.r2.dev` and a path of `geosx/integratedTests` yields
  `https://pub-<hash>.r2.dev/integratedTests/...` rather than doubling the bucket. Set `false`
  for providers whose public root sits above the bucket (GCS, custom domains mapped to the
  account, etc.).

- `runner_ca_bundle_host_paths`
  Map keyed by runner label, not by workflow role. Each value is a host-side certificate bundle
  path to bind into the container as `/certs/ca-bundle.crt`.
  Example: `{"streak2": "/etc/pki/tls/certs/ca-bundle.crt"}`.
  The resolver first looks for an exact runner-label match. If none is found, it falls back to
  the prefix before the first `-`, so a runner label like `streak2-32core` will reuse the
  `streak2` entry. Use `{}` for runners that do not need an LLNL-style outbound cert bundle.

- `runner_cuda_architectures`
  Map keyed by runner label, not by workflow role. Each value is forwarded to CMake as
  `CMAKE_CUDA_ARCHITECTURES` for CUDA jobs on that runner. The same exact-match then prefix
  fallback used by `runner_ca_bundle_host_paths` applies. Example: `{"streak2": "86"}`.

- `runners`
  Map from runner role to runner label. Required roles:
  - `default`
  - `cpu_heavy`
  - `integrated_tests`
  - `code_coverage`
  - `cuda`

  Each value is a single runner label string, e.g. `"ubuntu-22.04"` or `"streak2"`.

## Provider files (`.github/ci/providers/<storage_provider>.json`)

Required fields (all of them, always):

- `artifact_upload_command`
  Shell command that uploads a single file. `$UPLOAD_SRC` and `$UPLOAD_DST` are exported by the
  consumer before execution. The string is passed to `eval` — see the Security section below.

- `artifact_upload_pre_command`
  Shell command that runs once per job before the first upload, typically to authenticate the
  CLI and validate credentials. `$ARTIFACT_UPLOAD_CREDENTIALS_FILE` is exported by the consumer.
  The string is passed to `eval`.

- `artifact_upload_uri_root`
  URI scheme prefix the upload CLI expects, e.g. `"gs://"` or `"s3://"`. The consumer composes
  `${artifact_upload_uri_root}${integrated_tests_artifact_bucket_path}/${filename}` as
  `$UPLOAD_DST`.

## Security

The `artifact_upload_command` and `artifact_upload_pre_command` strings in every provider file
are passed to `eval` on the runner at upload time. They are executable shell, not inert
configuration. A change to those strings is a change to code that runs with the job's
environment — including credentials files, sccache tokens, and on self-hosted runners, the
runner host itself.

Review edits to `.github/ci/providers/*.json` with the same scrutiny you would apply to a
workflow YAML change. Never populate these strings from PR- or user-controllable sources.

## Secrets

Secret payloads (artifact upload credentials, sccache config) are supplied through GitHub
repo/org secrets and consumed in the reusable workflow. They are not part of the JSON files
described here.
