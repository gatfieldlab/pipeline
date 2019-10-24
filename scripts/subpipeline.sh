#!/usr/bin/env bash

# SUB - splits the jobs for each sample and calls MAPPER
# and performs STAR mapping and SAM2BAM for each sample

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

BASE="$1"
SAMPLEDB="$3"
FILEARRAY=""
starproc=6

source $CONFIG_FILE

logfile="${__dir}/${2}/logs${SUFFIX}/${BASE}.pipeline.SUB.log"
starlogfile="${__dir}/${2}/logs${SUFFIX}/${BASE}.pipeline.STAR.log"

declare -A errorlog=( ["Have files"]="not processed"
                      ["MAPPER status"]="not processed"
                      ["TOPHAT2 mapping"]="obsolete"
                      ["STAR mapping"]="not processed"
		      ["UMI deduplication"]="not processed")

log_and_exit() {
  if [ -f ${logfile} ]; then

    logs ${SUB} $WARNING"Found a log file, most likely from previous runs.\
 Overriding it!"

    rm -f ${logfile}
  fi
  for logtype in "${!errorlog[@]}"; do
    printf "%-16s : %s\n" "$logtype" "${errorlog[$logtype]}" >> ${logfile}
  done
  exit $1
}

# Type is set by BASE
: ${DEFAULT_TYPE:=$(GET_TYPE ${BASE})}
TY=$DEFAULT_TYPE

if [ "${SAMPLEDB}" == "skip" ]; then
  logs ${SUB} "<${BASE}> Skipping DB search by settings, calling MAPPER."
  errorlog["Have files"]="[ $(timestamp) ] OK - Skipped by --database argument"
else
  i=0
  raw_ext=${FILE_EXT["raw"]}
  logs ${SUB} "<${BASE}> Searching for files, using '${raw_ext}' in '${SAMPLEDB}'"
  while IFS= read -r file; do
      i=$(( $i+1 ))
      FILEARRAY=$(printf "$FILEARRAY\n${file:0:${#file}-${#raw_ext}}")
  done < <(grep -oE "${BASE}\S*${raw_ext}" $SAMPLEDB)
  if (( $i < 1 )); then
    logs ${SUB} $ERROR"<${BASE}> Could not find any files; terminating"
    errorlog["Have files"]="[ $(timestamp) ] ERROR - None found"
    log_and_exit 1
  fi
  logs ${SUB} "<${BASE}> Found $i files; calling MAPPER."
  errorlog["Have files"]="[ $(timestamp) ] OK - Found $i files in DB"
fi

# Spawn MAPPER processes
#-----------------------
cd $2
pipeline_bwt2_single.sh "${BASE}" "${FILEARRAY}"
if (( $? )); then

  logs ${SUB} $ERROR"<${BASE}> One or more of the MAPPER processes exited\
 with a non-zero error code; terminating"

  errorlog["MAPPER status"]="[ $(timestamp) ] ERROR -\
 One or more MAPPER processes failed"

  log_and_exit 1
else
  errorlog["MAPPER status"]="[ $(timestamp) ] OK - All finished"
fi

exit_code=0

  
# STAR genomic mapping
# --------------------
if (( $MAP_WITH_STAR )); then
  # Check if we have STAR and get its version
  star_ver=$( get_version "STAR" "--version" "${SUB}" )
  samtools_ver=$( get_version "samtools" "--version" "${SUB}" )  
  if [ -z "${star_ver}" ]; then
    exit_code=1
    errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - STAR program not found"
    logs ${SUB} ${ERROR}"STAR program not found [${BASE}]. STAR mapping will be\
 skipped!"
  fi
  if [ -z "${samtools_ver}" ]; then
    exit_code=1
    errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - samtools program not found"
    logs ${SUB} ${ERROR}"samtools program not found [${BASE}]. STAR mapping will be\
 skipped!"
  fi
  # Figure out what to use as input
  if (( ${#MAPTYPE[@]} == 0 )); then
    if (( ${#PRE_PROC_STEPS[@]} == 0 )); then
      FINALMAP=${INIT_PROC}
    else
      FINALMAP=${PRE_PROC_STEPS[${#PRE_PROC_STEPS[@]}-1]}
    fi
  else
    FINALMAP=${MAPTYPE[${#MAPTYPE[@]}-1]}
  fi
  star_input="fq/${BASE}${FILE_EXT[$FINALMAP]}"
  if [ ! -f ${star_input} ]; then
    exit_code=1
    errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - Input file not found"
    logs ${SUB} ${ERROR}"Input file '${star_input}' could not be found for\
 [${BASE}]. STAR mapping will be skipped!"
  fi
  if (( ! exit_code )); then
    # CPU reservation
    reserve_cpu ${starproc} "${SUB} Waiting for free CPU to map [${BASE}]\
 against genomic database using STAR"

    logs ${SUB} "<${BASE}> Entering now STAR genomic mapping"
    star_out=STAR"${SUFFIX}"/"${BASE}"/
    mkdir -p "${star_out}"
    readFileCmd=-
    use_gzip "${star_input}" && readFileCmd=zcat
    seedSearchStartLmax=50
    gtf_mode=""
    memory_mode="LoadAndKeep"
    [ ${filter_low[${TY}]+is_set} ] && seedSearchStartLmax="${filter_low[${TY}]}"
    [ ${STAR_GTF} ] && gtf_mode="--sjdbGTFfile ${STAR_GTF} --twopassMode Basic" && memory_mode="NoSharedMemory"
    logs ${SUB} "<${BASE}> Mapping against ${BOWTIE2X['genome']}\
 genomic database using STAR..."
    logs ${SUB} "<${BASE}> Using seedSearchStartLmax value = ${seedSearchStartLmax}"
    logs ${SUB} "<${BASE}> GTF support params: '${gtf_mode}'"
    star_cmd="STAR --runThreadN ${starproc} --genomeDir=${STAR_INDEX['genome']} \
    --readFilesIn ${star_input} --outFileNamePrefix ${star_out} \
    --readFilesCommand ${readFileCmd} --genomeLoad ${memory_mode} \
    --seedSearchStartLmax ${seedSearchStartLmax} ${gtf_mode} \
    --outSAMtype BAM SortedByCoordinate Unsorted --alignSJDBoverhangMin 1 \
    --alignIntronMax 1000000 --outFilterType BySJout --alignSJoverhangMin 8 \
    --limitBAMsortRAM 15000000000 --outReadsUnmapped Fastx"
    
    if (( LOGGING )); then
	cat >$starlogfile <<EOF
==== '${CONF_TYPE}' pipeline STAR log file ====

== STAR ==
${star_ver}

== samtools ==
${samtools_ver}

== Input file ==
${star_input}

== STAR variable params ==
runThreadN          ${starproc}
genomeDir           ${STAR_INDEX['genome']}
outFileNamePrefix   ${star_out}
readFilesCommand    ${readFileCmd}
seedSearchStartLmax ${seedSearchStartLmax}

== STAR mapper command line ==
${star_cmd}

EOF
    fi
    eval ${star_cmd}
    cur_err=$?
    if (( $cur_err )); then
      logs ${SUB} ${ERROR}"<${BASE}> Mapping against ${BOWTIE2X['genome']}\
 genomic database with STAR failed."
      errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - Mapping failed"
      if (( LOGGING )); then
        echo "
STAR mapping failed!" >>$starlogfile
      fi	  
    else
      logs ${SUB} "<${BASE}> Mapping against ${BOWTIE2X['genome']} genomic\
 database with STAR finished OK"
      echo STAR --genomeDir="${STAR_INDEX['genome']}" \
           --genomeLoad Remove > "${CLEANUP_HOOKS}"
      for logName in Log.out Log.final.out Log.progress.out; do
        cp ${star_out}${logName} ${star_out}Mapper.${logName}
      done
      logs ${SUB} "Indexing [${BASE}] genomic map"
      samtools index "${star_out}"Aligned.sortedByCoord.out.bam

      logs ${SUB} "Creating bedGraph files for [${BASE}]"

      star_cmd="STAR --runMode inputAlignmentsFromBAM --outWigType bedGraph \
      --inputBAMfile ${star_out}Aligned.sortedByCoord.out.bam \
      --outWigStrand Unstranded --runThreadN ${starproc} \
      --outFileNamePrefix ${star_out}"
      if (( LOGGING )); then
	cat >>$starlogfile <<EOF
== STAR bedGraph command line ==
${star_cmd}
 
EOF
      fi
      eval ${star_cmd}
      cur_err=$?
      if (( $cur_err )); then
        logs ${SUB} ${ERROR}"<${BASE}> bedGraph files could not be created"
        errorlog["STAR mapping"]="[ $(timestamp) ] ERROR - bedGraph failed"
      fi
    fi
    release_cpu ${starproc}
    if (( $cur_err )); then
      exit_code=1
    else
      errorlog["STAR mapping"]="[ $(timestamp) ] OK - All finished"
    fi
  fi 
else
  errorlog["STAR mapping"]="[ $(timestamp) ] OK - Skipped by configuration"
fi


# UMI deduplication
# -----------------
if (( ${#DEDUP[@]} )); then
  # Check if we have umi_tools and get its version
  umi_ver=$( get_version "umi_tools" "--version" "${SUB}" )
  samtools_ver=$( get_version "samtools" "--version" "${SUB}" )
  if [ -z "${umi_ver}" ]; then
    exit_code=1
    errorlog["UMI deduplication"]="[ $(timestamp) ] ERROR - umi_tools program not found"
    logs ${SUB} ${ERROR}"UMI-tools (umi_tools) not found [${BASE}]. UMI\
 deduplication will be skipped!"
  fi
  if [ -z "${samtools_ver}" ]; then
    exit_code=1
    errorlog["UMI deduplication"]="[ $(timestamp) ] ERROR - samtools program not found"
    logs ${SUB} ${ERROR}"samtools program not found [${BASE}]. UMI deduplication\
 will be skipped!"
  fi
  # loop through possible input types and associated dedup methods
  for maptype in "${!DEDUP[@]}"; do
    # Figure out what to use as input
    if [ "${maptype}" = "genome" ]; then
      umi_base="STAR${SUFFIX}"
      umi_input="${umi_base}/${BASE}/Aligned.sortedByCoord.out.bam"
      umi_log="${umi_base}/${BASE}/umi_tools.log"
      umi_out="${umi_base}/${BASE}/Aligned"
      deduplogfile="${umi_base}/${BASE}/dedup.log"
    else
      umi_base="${LIBEXTS[${maptype}]}"
      umi_input="sam/${BASE}.${umi_base}.sorted.bam"
      umi_log="logs${SUFFIX}/${BASE}.${umi_base}.umi_tools.log"
      umi_out="sam/${BASE}.${umi_base}"
      deduplogfile="logs${SUFFIX}/${BASE}.${umi_base}.dedup.log"
    fi
    if [ ! -f "${umi_input}" ]; then
      exit_code=1
      errorlog["UMI deduplication"]="[ $(timestamp) ] ERROR - Input file not found"
      logs ${SUB} ${ERROR}"Input file '${umi_input}.bam' could not be found for\
 [${BASE}]. Deduplication will be skipped!"
    fi
    # Setup DEDUP type + cell-barcode splitting option
    cur_umi_opts=($(echo ${DEDUP[${maptype}]} | tr "|;\-," " "))
    umi_type="${cur_umi_opts[0]}"
    cell_split=0
    if [ "${#cur_umi_opts[@]}" -eq 2 ]; then
      if [[ " s split " =~ " ${cur_umi_opts[1]} " ]]; then
	cell_split=1
      fi
    fi
    case "${umi_type}" in
      skip)
	umi_pipe=""
	umi_bam=""
	split_pipe="samtools view ${umi_input} | awk -v filebase=${umi_out} \
 '{split($1,barcode,\"_\"); print \$0 > filebase \"_split_\" barcode[2] \".sam\"}'"
     	;;
      dedup)
	umi_out="${umi_out}.dedup"
	umi_pipe="umi_tools dedup -I ${umi_input} --method=directional --per-cell \
 --read-length --log=${umi_log} --output-bam "
	umi_bam="--stdout ${umi_out}.bam"
	split_pipe="| samtools view -h - | awk -v filebase=${umi_out} \
 '{split($1,barcode,\"_\"); print \$0 > filebase \"_split_\" barcode[2] \".sam\"}'"
	;;
      group)
	umi_out="${umi_out}.group"
	umi_pipe="umi_tools group -I ${umi_input} --method=directional --per-cell \
 --read-length --log=${umi_log} --output-bam "
	umi_bam="--stdout ${umi_out}.bam"
	split_pipe="| samtools view -h - | awk -v filebase=${umi_out} \
 '{split($1,barcode,\"_\"); print \$0 > filebase \"_split_\" barcode[2] \".sam\"}'"
	;;
      filter)
	umi_out="${umi_out}.filter"
	umi_pipe="umi_tools group -I ${umi_input} --method=directional --per-cell \
 --read-length --log=${umi_log$ --output-bam | samtools view -h | filterUmiFromSam \
 2>UMILOG"
	umi_bam="| samtools view -bS - > ${umi_out}.bam"
	split_pipe="| awk -v filebase=${umi_out} \
 '{split($1,barcode,\"_\"); print \$0 > filebase \"_split_\" barcode[2] \".sam\"}'"
	;;
      *)
	# anything else, including empty options raise error
	exit_code=1
	errorlog["UMI deduplication"]="[ $(timestamp) ] ERROR - improper DEDUP options"
	logs ${SUB} ${ERROR}"DEDUP option was not set properly for\
 [${BASE}]. Deduplication will be skipped!"
	;;
    esac
    logs ${SUB} "<${BASE}> Deduplicating ${umi_base}"
    if (( $cell_split )); then
      umi_command="${umi_pipe}${split_pipe}"
    else
      umi_command="${umi_pipe}${umi_bam}"
    fi
    if (( LOGGING )); then
	cat >$deduplogfile <<EOF
==== '${CONF_TYPE}' pipeline DEDUPLICATION log file ====

== umi_tools ==
${umi_ver}

== samtools ==
${samtools_ver}

== Input file ==
${umi_input}

== umi_tools params ==
function           ${umi_type}
barcode splitting  ${cell_split}

== deduplication command line ==
${umi_command}

EOF
    fi
    eval ${umi_command}
    cur_err=$?
    if (( $cur_err )); then
      logs ${SUB} ${ERROR}"<${BASE}> Deduplication failed."
      errorlog["UMI deduplication"]="[ $(timestamp) ] ERROR - Deduplication failed"
      if (( LOGGING )); then
        echo "
UMI deduplication failed!" >>$deduplogfile
      fi	  
    else
      logs ${SUB} "<${BASE}> UMI deduplication finished successfully"
      if (( $cell_split )); then
	logs ${SUB} "<${BASE}> Splitting BAM by adaptor barcodes"
	find . -type f -name "${umi_out}_?*.sam" | \
	  xargs -i bash -c \
		'cat "${1%_*}_.sam" {} | samtools view -bS - > "${1%.sam}.bam"' - '{}'
	cur_err=$?
	if (( $cur_err )); then
	  logs ${SUB} ${ERROR}"<$BASE> Splitting by barcodes failed!"
	  errorlog["UMI deduplication"]="[ $(timestamp) ] ERROR - Barcode splitting failed"
	  if (( LOGGING )); then
	    echo "
Splitting by adaptor barcodes failed!" >>$deduplogfile
	  fi
	else
	  logs ${SUB} "<${BASE}> Barcode splitting finished OK"
	fi
      fi
    fi
    if (( $cur_err )); then
      exit_code=1
    else
      errorlog["UMI deduplication"]="[ $(timestamp) ] OK - All finished"
    fi
  done
fi

cd ..

log_and_exit ${exit_code}
