#!/bin/bash

# WHITELIST - creates a whitelist of barcodes using umi_tools
# -----------------------------------------------------------

:<<'LICENSE'
"pipeline" is a collection of shell scripts that together provide
a configurable and semi-automated pipeline to trim, filter and map
large sequence files produced by Next Generation Sequencing platforms.

Copyright (C) 2015-2019  A. Bulak Arpat

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
LICENSE

source pipeline_common.sh

config_file="pipeline_conf.sh"
log_file="pipeline.WHITELIST.log"
whitelist_file="whitelist.pipeline"

while [[ $# > 1 ]]; do
  key="$1"

  case $key in
    -c|--config-file)
      config_file="$2"
      shift
      ;;
    -l|--log-file)
      log_file="$2"
      shift
      ;;
    -f|--file-name)
      whitelist_file_cli="$2"
      shift
      ;;
    -o|--out-dir)
      out_dir_cli="$2"
      shift
      ;;
    *)
    # unknown option
    logs ${META} ${WARNING}"Skipping unknown parameter $key"
    ;;
  esac
  shift   # past argument or value
done

source ${config_file}

[ ! -z ${UMI_WHITELIST} ] && whitelist_file="${UMI_WHITELIST}"
[ ! -z ${whitelist_file_cli} ] && whitelist_file="${whitelist_file_cli}"
[ ! -z ${DIR_NAME["umi"]} ] && out_dir="${DIR_NAME["umi"]}"
[ ! -z ${out_dir_cli} ] && out_dir="${out_dir_cli}"
log_file=${out_dir}/${log_file}

if command -v cutadapt >/dev/null 2>&1; then
  cut_ver=$(cutadapt --version 2>&1)
else
  logs ${WHITELIST} ${ERROR}"Could not find 'cutadapt'. WHITELIST terminates."
  exit 1
fi

if command -v umi_tools >/dev/null 2>&1; then
  umi_ver=$(umi_tools --version 2>&1)
else
  logs ${WHITELIST} ${ERROR}"Could not find 'umi_tools'. WHITELIST terminates."
  exit 1
fi

pipe="cutadapt -a ${ADAPTER} ${CUTADAPT_OPTS} --quiet - |\
 umi_tools whitelist ${UMI_EXTRACT_OPTS} --log2stderr\
 2>>${log_file} 1>${out_dir}/${whitelist_file}"

cat >${log_file} <<EOF
==== Logfile for WHITELIST ====

== cutadapt ==
${cut_ver}

== umi_tools ==
${umi_ver}

== whitelisting command line ==
${pipe}

EOF

eval ${pipe}
cur_err=$?
if (( ${cur_err} )); then
  logs ${WHITELIST} ${ERROR}"Whitelisting failed."
  echo "[ $(timestamp) ] ERROR - whitelisting failed." >> ${log_file}
else
  logs ${WHITELIST} "Whitelisting successfully finished."
  echo "[ $(timestamp) ] OK - finished with success" >> ${log_file}
fi

exit ${cur_err}
