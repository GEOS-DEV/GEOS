#!/bin/bash

uncrustifyCommand=$1
changedFiles="$(git status -s --untracked-files=no | grep -e '.hpp' -e '.cpp' | cut -c 4- | xargs)"

echo $uncrustifyCommand
echo $changedFiles

