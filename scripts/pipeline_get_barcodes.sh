#!/bin/bash

# WHITELIST - creates a whitelist of barcodes using umi_tools
# -----------------------------------------------------------

:<<'LICENSE'
"pipeline" is a collection of shell scripts that together provide
a configurable and semi-automated pipeline to trim, filter and map
large sequence files produced by Next Generation Sequencing platforms.

Copyright (C) 2015-2020  A. Bulak Arpat

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

usage(){
  cat <<EOF

${USELINE}

Usage: $0 [options] SAMPLE_FILE

-s|--seq	     number of max sequences that will be used
		     for whitelisting from each sample. If
		     there are multiple files for a sample
		     they will be merged read-by-read upto
		     this number. 0 implies no limit.
                     Default is 10 million.
-n|--concurrent      number of samples that will be
                     processed concurrently. default 4
-c|--config-file     configuration file to use. default
                     is 'pipeline_conf.sh'
-f|--file-name       base filename for whitelists. Final
                     names will be SAMPLEBASE_basefilename.
		     If not given, it will be read from
		     config-file (preferred way).
-l|--log-file        logfile that will be used for whitelising
		     logs. Final names will be
		     SAMPLEBASE_logfile.
		     Default is pipeline.WHITELIST.log
-r|--raw-dir	     path to directory where raw sequence
		     files are located. If not given, it will
		     be read from main config (preferred).
		     Currently, it is '${RAW_DIR}'
-o|--out-dir	     output dir. If not given read from
		     main config (preferred).

${LICENSE}

EOF
	exit 1
}

num_seq=10000000
config_file="pipeline_conf.sh"
log_file="pipeline.WHITELIST.log"
whitelist_file="whitelist.pipeline"
raw_dir="${RAW_DIR}"
limit=4

while [[ $# > 1 ]]; do
  key="$1"

  case $key in
    -o|--out-dir)
      out_dir_cli="$2"
      shift
      ;;
    -s|--seq)
      num_seq="$2"
      shift
      ;;
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
    -r|--raw-dir)
      raw_dir="$2"
      shift
      ;;
    -n|--concurrent)
      limit="$2"
      shift
      ;;
    *)
      # unknown option
      logs ${META} ${WARNING}"Skipping unknown parameter $key"
      ;;
  esac
  shift   # past argument or value
done

[ -z $1 ] && usage

samples=$1
      
source ${config_file}

[ ! -z ${UMI_WHITELIST} ] && whitelist_file="${UMI_WHITELIST}"
[ ! -z ${whitelist_file_cli} ] && whitelist_file="${whitelist_file_cli}" 
[ ! -z ${DIR_NAME["umi"]} ] && out_dir="${DIR_NAME["umi"]}"
[ ! -z ${out_dir_cli} ] && out_dir="${out_dir_cli}"
[ ! -z ${out_dir} ] && outdir="--out-dir ${out_dir}"

export raw_dir
export config_file
export whitelist_file
export log_file
export num_seq
export out_dir

cat ${samples} | xargs -n 1 -P ${limit} -I SAMPLE\
 bash -c "echo 'Processing SAMPLE ...' && "\
"merge_fastq --max-seq ${num_seq} ${raw_dir}/SAMPLE*.fastq.gz | pipeline_whitelist.sh "\
"--config-file ${config_file} --file-name SAMPLE_${whitelist_file} "\
"--log-file SAMPLE_${log_file} ${out_dir} && echo 'Finished processing SAMPLE'"
echo "barcode estimation is finished."
exit 0
