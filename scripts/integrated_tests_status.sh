#!/bin/bash

if [[ $TRAVIS_PULL_REQUEST == "" || $TRAVIS_PULL_REQUEST == false ]]; then
  echo
  echo "SUCCESS: Not a pull request."
  echo 
  echo "Check passed!"

  exit 0
fi

echo
echo "Checking if Pull Request has the \"requires rebaseline\" label"
echo "or the integratedTests submodule has been updated..."
echo

needs_rebase=$(curl -H "Accept: application/vnd.github+json" \
  https://api.github.com/repos/GEOS-DEV/GEOS/pulls/$TRAVIS_PULL_REQUEST | \
  jq '[.labels[] | .name | contains("rebaseline")] | any')

echo
echo "\"Requires rebaseline\" label found in Pull Request? $needs_rebase"
echo

# Checks the first 100 modified files (Github max limit)
tests_modified=$(curl -H "Accept: application/vnd.github+json" \
  https://api.github.com/repos/GEOS-DEV/GEOS/pulls/$TRAVIS_PULL_REQUEST/files?per_page=100 | \
  jq '[.[] | .filename | contains("integratedTests")] | any')

echo
echo "integratedTests submodule modified by Pull Request? $tests_modified"
echo

if [[ $needs_rebase == false && $tests_updated == false ]]; then
  echo
  echo "SUCCESS: IntegratedTests were not modified."
  echo 
  echo "Check passed!"

  exit 0
else
  echo
  echo "IntegratedTests were changed and/or require rebaseline, checking if"
  echo "Pull Request is being tested by the head of \"next\" branch..."
  echo "(Checking \"next\" head commit message for PR's branch name"
  echo "or PR number)"
  echo
fi

commit_msg=$(curl -H "Accept: application/vnd.github+json" \
  https://api.github.com/repos/GEOS-DEV/GEOS/branches/next | \
  jq '.commit.commit.message')

echo
echo "Head commit message of \"next\" branch is: $commit_msg"
echo "    Pull Request branch is: $TRAVIS_BRANCH"
echo "    Pull Request number is: $TRAVIS_PULL_REQUEST"

if [[ "$commit_msg" == *"$TRAVIS_BRANCH"* || \
	  "$commit_msg" == *"$TRAVIS_PULL_REQUEST"* ]]; then
  echo
  echo "Head of \"next\" branch is testing this Pull Request."
  echo "Checking status of commit..."
else
  echo
  echo "FAILURE: Head of \"next\" branch has not tested this Pull Request."
  echo
  echo "Please notify an approved developer to merge your Pull Request into"
  echo "\"next\" branch for integration testing."

  exit 1
fi

commit_sha=$(curl -H "Accept: application/vnd.github+json" \
  https://api.github.com/repos/GEOS-DEV/GEOS/branches/next | \
  jq ".commit.sha")

echo
echo  "Head commit SHA of \"next\" branch is: $commit_sha"
echo

commit_state=$(curl -H "Accept: application/vnd.github+json" \
  https://api.github.com/repos/GEOS-DEV/GEOS/commits/$commit_sha/status | \
  jq ".state")


echo
echo  "Head commit state of \"next\" branch is: $commit_state"
echo

if [[ "$commit_state" == *"success"* ]]; then
  echo
  echo  "SUCCESS: Integrated tests succeeded on \"next\" branch."
  echo
  echo  "Check passed!"
  echo

  exit 0
else
  echo
  echo  "FAILURE: Integrated tests failed to pass on \"next\" branch."
  echo
  echo  "Please modify, rebaseline the integratedTests."
  echo  "Notify an approved developer to re-merge your updated Pull Request into"
  echo  "\"next\" branch for integration testing."
  echo

  exit 1
fi
