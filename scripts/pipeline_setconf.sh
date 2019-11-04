#!/usr/bin/env bash

source pipeline_common.sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
HEADER_FORMAT='%-12s%-60s%-20s\n'
ROW_FORMAT='\e[1m%-12s\e[0m%-60s%-20s\n'
declare -A CONFIG_DESCR
declare -A CONFIG_FILES

usage() {
  cat <<EOF
${USELINE}

This script is used to list all available configuration files and copies
selected ones.

Usage: $0 COMMAND [CONFIG_NAME]

COMMAND        Mandatory. Has to be one of:
               list|ls : Lists all available pipeline-config files with short
                         descriptions. 
               copy|cp : Copies the config file with CONFIG_NAME
               help    : Prints this info

CONFIG_NAME    Name of the config as printed in the list;
               spaces are NOT allowed.

EOF
}

get_single_arg() {
  local arg_as_array=( $1 )
  local array_len=${#arg_as_array[@]}
  [ "${array_len}" -gt "1" ] && exit
  echo "${arg_as_array[0]}"
}


get_attr() {
  search_term="(?<=${2}=\\\").*(?=\\\")"
  # echo "Looking for ${2} in ${1}"
  # echo "${search_term}"
  echo $(grep -oP ${search_term} ${1})
}

get_all_config() {
  config_paths=($(ls ${DIR}/pipeline_conf-*.sh))
  for ((i=0;i<${#config_paths[@]};i=i+1)); do
    cur_path="${config_paths[${i}]}"
    cur_descr=$(get_attr "${cur_path}" CONF_TYPE)
    cur_name=$(get_attr "${cur_path}" CONF_NAME)
    # cur_file=${cur_path##*/}
    if [ ! -z ${cur_name} ]; then
      CONFIG_DESCR[${cur_name}]="${cur_descr}"
      CONFIG_FILES[${cur_name}]="${cur_path}"
    fi
  done
}

list_all_config() {
  printf ${HEADER_FORMAT} "Name" "Description" "Filename"
  printf ${HEADER_FORMAT} "----" "-----------" "--------"
  for config_name in "${!CONFIG_FILES[@]}"; do
    config_path="${CONFIG_FILES[${config_name}]}"
    fname=${config_path##*/}
    printf ${ROW_FORMAT} "${config_name}" "${CONFIG_DESCR[${config_name}]}" "${fname}"
  done
}

# Set up command
cur_command=$(get_single_arg "${1}")
cur_config=$(get_single_arg "${2}")
get_all_config

[ -z "${cur_command}" ] && usage && exit 1
case "${cur_command}" in
  list|ls)
    list_all_config
    ;;
  help)
    usage
    ;;
  copy|cp)
    if [ -z "${CONFIG_FILES[${cur_config}]}" ]; then
      echo "Could not find a config file corresponding to '${cur_config}'"
      exit 1
    else
      config_path="${CONFIG_FILES[${cur_config}]}"
      fname=${config_path##*/}
      cp "${config_path}" . && echo "Created '${fname}' for '${cur_config}'"
    fi
    ;;
  *)
    echo "Command can be one of [list, copy, help]"
    usage
    ;;
esac

exit 0
