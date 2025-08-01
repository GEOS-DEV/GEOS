# GEOS and geosPythonPackages Integration Testing

This document describes the CI integration testing setup between GEOS and geosPythonPackages repositories to ensure compatibility before merging changes.

## Overview

The GEOS project uses tools from the geosPythonPackages repository. To prevent breaking changes in geosPythonPackages from affecting GEOS, we have implemented a cross-repository CI testing system.

## How It Works

### 1. GEOS CI Integration

GEOS CI now includes a Python integration test that:
- Runs alongside other GEOS CI jobs
- Tests that geosPythonPackages can be properly installed and used
- Validates key Python tools like `preprocess_xml`, `format_xml`, and `mesh-doctor`
- Ensures Python packages can be imported correctly

### 2. geosPythonPackages CI Requirements

When creating a PR in geosPythonPackages, a strict sequential process is enforced:

#### Sequential Workflow (Enforced Order):

**Step 1: Standard CI Tests Run First**
- All standard geosPythonPackages tests must pass first
- Includes: semantic check, build, lint, pytest
- GEOS integration will NOT start until these pass

**Step 2: Label Check (Only After Standard Tests Pass)**
- System checks for `test-geos-integration` label
- ❌ **FAILS** if label is missing
- ✅ **PROCEEDS** if label is present

**Step 3: GEOS Integration Testing (Automatic Trigger)**
- geosPythonPackages CI triggers GEOS workflows
- GEOS tests run using the PR branch
- Results are monitored and reported back

**Step 4: Final Validation**
- ✅ **VALID**: All steps pass (Standard CI + Label + GEOS tests)
- ❌ **INVALID**: Any step fails

#### Benefits of Sequential Approach:
- Developers fix basic issues before expensive GEOS testing
- Clear feedback at each step
- Efficient resource usage
- Protected GEOS repository

## For geosPythonPackages Contributors

### Required Steps for Every PR (Sequential Process):

1. **Create your PR** and let standard CI tests run first
2. **Fix any standard CI issues** (semantic check, build, lint, test) - GEOS integration will NOT start until these pass
3. **Add the `test-geos-integration` label** only AFTER standard tests pass
4. **Wait for GEOS integration tests** to complete automatically (15-30 minutes)
5. **Address any GEOS failures** if the integration tests fail
6. **Proceed with normal review** once all tests pass

### Sequential Workflow Enforcement:

**Step 1: Standard CI Tests (MUST pass first)**
- Semantic PR title check
- Build and install all packages
- Lint with yapf
- Test with pytest
- ⚠️ **GEOS integration will NOT start until these pass**

**Step 2: Label Check (after standard tests pass)**
- Check for `test-geos-integration` label
- ❌ **FAIL** if label is missing
- ✅ **PROCEED** if label is present

**Step 3: GEOS Integration Testing (triggered automatically)**
- Trigger GEOS CI with your branch
- Wait for GEOS tests to complete
- Report results back to PR

**Step 4: Final Validation**
- ✅ **VALID** = Standard tests ✅ + Label ✅ + GEOS tests ✅
- ❌ **INVALID** = Any step fails

### Understanding Test Results:

- **Standard Tests Failing**: Fix your code first, don't add the label yet
- **Label Missing**: After standard tests pass, add the `test-geos-integration` label
- **GEOS Tests Running**: You'll see a "Test GEOS Integration" workflow running
- **GEOS Tests Failed**: Check the workflow logs to understand what broke in GEOS
- **GEOS Tests Passed**: Your changes are compatible with GEOS

## For GEOS Contributors

### New CI Jobs:

1. **python_integration_test**: Runs as part of standard GEOS CI
   - Tests current geosPythonPackages main branch
   - Ensures basic Python package functionality

2. **test_geospythonpackages_integration**: Can be triggered manually or from geosPythonPackages
   - Tests specific geosPythonPackages branches/PRs
   - Used for cross-repository validation

## Troubleshooting

### geosPythonPackages PR Failing

If your geosPythonPackages PR fails GEOS integration:

1. **Check the GEOS workflow logs** (link provided in the failure message)
2. **Common issues**:
   - API changes that break GEOS usage
   - Missing dependencies
   - Import errors
   - Tool execution failures
3. **Fix the issues** in your geosPythonPackages branch
4. **Push updates** - this will re-trigger the GEOS integration test

### Missing Required Secrets

The geosPythonPackages repository needs a `GEOS_INTEGRATION_TOKEN` secret with permissions to trigger workflows in the GEOS repository.

### Manual Testing

You can manually trigger GEOS integration testing by:
1. Going to the GEOS repository's Actions tab
2. Finding the "Test GEOS Integration from geosPythonPackages" workflow
3. Clicking "Run workflow" and specifying your branch

## Files Modified/Created

### GEOS Repository:
- `.github/workflows/python_integration_test.yml` - New reusable workflow for Python integration testing
- `.github/workflows/test_geospythonpackages_integration.yml` - New workflow that can be triggered from geosPythonPackages
- `.github/workflows/ci_tests.yml` - Updated to include Python integration test
- `scripts/test_python_integration.sh` - New test script for Python package validation
- `src/docs/GEOS_PYTHON_INTEGRATION_TESTING.md` - Documentation for GEOS contributors

### geosPythonPackages Repository:
- `.github/workflows/python-package.yml` - **Updated with complete sequential workflow logic**
  - Standard CI validation
  - Label checking (only after standard tests pass)
  - GEOS integration triggering
  - Final validation with clear VALID/INVALID results
- `GEOS_INTEGRATION_TESTING.md` - User guide for geosPythonPackages contributors

### Key Implementation Details:
- **Sequential enforcement**: Standard tests → Label check → GEOS tests → Final validation
- **Cross-repository triggering**: Uses GitHub API to trigger GEOS workflows
- **Clear status reporting**: Explicit VALID/INVALID results with detailed explanations
- **Resource efficiency**: GEOS CI only runs when standard tests pass

## Benefits

1. **Prevents Breaking Changes**: geosPythonPackages changes are tested against GEOS before merging
2. **Automated Process**: No manual coordination needed between repositories
3. **Clear Feedback**: Contributors get immediate feedback if their changes break GEOS
4. **Modular Design**: GEOS Python testing is separate from other CI jobs
5. **Flexible**: Can test any geosPythonPackages branch against GEOS

## Future Enhancements

1. **Automatic label addition** based on which files are modified
2. **Integration with PR status checks** for clearer feedback
3. **Testing against multiple GEOS branches** if needed
4. **Performance optimization** to reduce test time
