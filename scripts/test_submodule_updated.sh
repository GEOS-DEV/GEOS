#!/bin/bash

# Initialize associative array for PR branch
declare -a pr_hashes_array=( $(git submodule status | awk '{print $1}') )
declare -a pr_names_array=( $(git submodule status | awk '{print $2}') )
declare -A pr_name_hash_dict

pr_length=${#pr_names_array[@]}

for (( i=0; i<$pr_length; i++))
do
  key="$(basename ${pr_names_array[$i]})"
  value="$( echo ${pr_hashes_array[$i]} | tr -cd [:alnum:] )"
  pr_name_hash_dict[$key]=$value
done

# Initialize associative array for develop branch
git checkout --quiet develop
git pull --quiet
declare -a dev_hashes_array=( $(git submodule status | awk '{print $1}') )
declare -a dev_names_array=( $(git submodule status | awk '{print $2}') )
declare -A dev_name_hash_dict

dev_length=${#dev_names_array[@]}

for (( i=0; i<$dev_length; i++))
do
  key="$(basename ${dev_names_array[$i]})"
  value="$( echo ${dev_hashes_array[$i]} | tr -cd [:alnum:] )"
  dev_name_hash_dict[$key]=$value
done

# Check that hashes are the same between submodules in PR and develop branches.
# Returns exit code 0 if the hashes for every submodule in PR that exists in
# the develop branch are the same
exit_code=0
unsync_submodules=()

for key in "${!pr_name_hash_dict[@]}"
do
  if [ ${dev_name_hash_dict[$key]+_} ]
  then
    if [ "${pr_name_hash_dict[$key]}" == "${dev_name_hash_dict[$key]}" ]
    then
      echo "PR and develop branch have the same hashes for submodule"\
           "$key : ${pr_name_hash_dict[$key]}"
    else
      echo "PR and develop branch have different hashes for submodule $key:"
      echo "---- PR branch has hash ${pr_name_hash_dict[$key]}"
      echo "---- develop branch has hash ${dev_name_hash_dict[$key]}"
      unsync_submodules+=( "$key" )
      exit_code=1
    fi
  else
    echo "Submodule $key not found in develop branch"
  fi

done

echo $'\n'
if [ $exit_code -eq 1 ]
then
  echo "##vso[task.logissue type=error]This PR has the following submodules"\
       "that are out of sync with develop : ${unsync_submodules[@]}"
  echo "##vso[task.logissue type=error]FAILURE : Please make sure your branch"\
       "is up to date with develop."\
       "Merge any submodule changes into the submodule's develop or"\
       "master branch before merging this PR with the main GEOSX repository."
else
  echo "SUCCESS : PR submodules are up to date!"
fi

exit $exit_code

