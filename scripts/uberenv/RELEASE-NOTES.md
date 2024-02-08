# Uberenv Software Release Notes

Notes describing significant changes in each Uberenv release are documented
in this file.

The format of this file is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

The Uberenv project release numbers follow [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## Unreleased

### Added
- Adds `--skip-setup` option, which skips installing and setting up Spack and installs only
  (using pre-setup Spack). Useful for air-gapped builds.
- Allow `--package_final_phase` command to be used for install build option, which controls
  after which phase Spack should stop.
- Allow projects to force specifying `--prefix` on command line via new project option:
  `force_commandline_prefix`.  This is useful when the default `uberenv_libs` can not exist
  inside of your project's source repository.
- Adds support for Windows builds using [Vcpkg].
- Adds the `--triplet` command line argument for setting the Vcpkg build configuration.
- Adds the `--vcpkg-ports-path` command line argument for setting the path to the vcpkg ports directory.
- Adds support for the spack `--reuse` option for more relaxed tpl/dependency matching.
- Adds support for `--upstream` with `--setup-only`.
- Allow to prevent Uberenv from installing clingo setting `spack_setup_clingo`
  to `no` in `.uberenv_config.json`.
- Adds the `--spack-debug` option to run spack spec/install commands in debug mode.
- Adds the `--spack-allow-deprecated` option, to allow spack to build packages marked deprecated.

### Changed
- All spack specs are now expressed inside single quotes to protect the parsing of complex flags.
- Added ability to have multiple packages directories that will get copied into spack on top of
  each other via project configuration option: `spack_packages_path`
- Pretty print various options to screen for readability
- Allow `.uberenv_config.json` to live at the same level as `uberenv.py`
- No longer removes symlinks when using the directory of `uberenv.py`
- Reduce Spack's git history to a bare minimum
- Better error message for out-of-date `pip`, better documentation for `spack_concretizer` setting
- Uberenv now requires python version 3.3 or above.

### Fixed


[Vcpkg]: https://github.com/microsoft/vcpkg
