# test-project
Tests build options in linux and macos. Each build option (uberenv-pkg, install, dev-build) is tested in linux, and currently only uberenv-pkg is tested in macos (see `.github/workflows`). The `--package-final-phase` option is also tested for dev-build and install, since it is a commonly used option in projects that utilize Uberenv.

# Testing locally
If you wish to test locally, look at `.github/workflows`, and follow a job's "Run Uberenv" commands (e.g. "build_uberenv_mode").
