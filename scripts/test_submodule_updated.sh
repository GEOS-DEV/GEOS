#!/bin/bash

# Submodules not checking for
declare -ar exclusion_list=( "blt" "PVTPackage")
echo "Submodules that are excluded from sync test : ${exclusion_list[@]}"

# Initialize PR submodule hashes
declare -ar pr_hashes_array=( $(git submodule status | awk '{print $1}') )

# Initialize submodule paths
declare -ar paths_array=( $(git submodule status | awk '{print $2}') )

length=${#paths_array[@]}

# Returns exit code 0 if the hash for every submodule in the PR is equal
# to the hash of each submodule's main branch (develop or master).
# Note: See "exclusion_list" for submodules that are exempted.
exit_code=0
unsync_submodules=()

for (( i=0; i<$length; i++))
do
  # Just the submodule name
  module_name="$(basename ${paths_array[$i]})"
  
  # Check if submodule is excluded from check
  excluded=0
  for ex in "${exclusion_list[@]}"
  do
    if [ "$module_name" = "$ex" ]
    then
      excluded=1
    fi
  done

  # Check hashes
  if [ $excluded -eq 0 ]
  then
    # Do quick 1 second pull of submodule to get .git files.
    timeout 1s git submodule update --quiet --init ${paths_array[$i]}

    # Determine if develop or master is main branch of submodule.
    if $( cd ${paths_array[$i]} && \
      git rev-parse --quiet --verify origin/develop > /dev/null)
    then 
      main_branch="origin/develop"
    else
      main_branch="origin/master"
    fi

    # Submodule main hash
    main_hash="$( cd ${paths_array[$i]} && git rev-parse $main_branch )"

    # PR hash with prefixed character removed
    pr_hash="$( echo ${pr_hashes_array[$i]} | tr -cd [:alnum:] )"

    if [ $pr_hash == $main_hash ]
    then
      echo "PR branch and $main_branch have the same hashes for submodule"\
           "$module_name : $pr_hash"
    else
      echo "PR branch and $main_branch have different hashes for submodule"\
           "$module_name:"
      echo "---- PR branch has hash $pr_hash"
      echo "---- $main_branch branch has hash $main_hash"
      unsync_submodules+=( "$module_name" )
      exit_code=1
    fi
  fi
done

echo $'\n'
if [ $exit_code -eq 1 ]
then
  echo "##vso[task.logissue type=error]This PR has the following submodules"\
       "that are out of sync with master or develop : ${unsync_submodules[@]}"
  echo "##vso[task.logissue type=error]FAILURE : Please make sure your branch"\
       "is up to date with develop."\
       "Merge any submodule changes into the submodule's develop or"\
       "master branch before merging this PR with the main GEOSX repository."
else
  echo "SUCCESS : PR submodules are up to date!"
fi

exit $exit_code
