#!/usr/bin/env bash

:<<'LICENSE'
"pipeline" is a collection of shell scripts that together provide
a configurable and semi-automated pipeline to trim, filter and map
large sequence files produced by Next Generation Sequencing platforms.

Copyright (C) 2015  A. Bulak Arpat

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
# Type of experiment (not too long, human readable, only for logs etc)
CONF_TYPE="Regular Ribosome Footprint Profiling"

# Sample configuration
# ====================

DEFAULT_TYPE=""
#function GET_TYPE() { get_read_type_from_first_split "0" "$@"; }
function GET_TYPE() { get_read_type_from_first_chrs 2 "$@"; }

# Should we be logging? 1 is yes, 0 is no. No other possibilities. 0 not tested
LOGGING=1

# Pre-processing configuration
# ============================

# Where to start, normally 'raw'
INIT_PROC="raw"

# Pre-processing steps to do in order
# trim     : Adaper removal
# s-filter : Filter by read size
# q-filter : Filter by read quality
declare -a PRE_PROC_STEPS=("trim" "s-filter" "q-filter")

# Trimming configuration - if 'trim' is enabled
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
CUTADAPT_OPTS="--match-read-wildcards"

# Size filtering configuration - if 's-filter' is enabled
declare -A filter_low=( ["RP"]="26" ["TR"]="21" )
declare -A filter_high=( ["RP"]="35" ["TR"]="60" )

# Quality filtering configuration - if 'q-filter' is enabled
FASTQ_QUALITY_FILTER_OPT="-Q33 -q 30 -p 90"

# File extensions if output is enabled - beware of gzip or not !!
declare -A FILE_EXT=( ["raw"]=".fastq.gz"
                      ["trim"]="_trimmed.fastq.gz"
                      ["s-filter"]="_s_filtered.fastq.gz"
                      ["q-filter"]="_q_filtered.fastq.gz" )
# Log file extensions - will always log if logging is enabled
declare -A LOG_EXT=(  ["raw"]=""
                      ["trim"]="_trimming.log"
                      ["s-filter"]="_size_filter.log"
                      ["q-filter"]="_qual_filter.log" )
# Directory structure - ugly as we can't export hash arrays in bash
declare -A DIR_NAME=( ["raw"]=${RAW_DIR}
                      ["trim"]=${TRIMMED_DIR}
                      ["s-filter"]=${FILTERED_DIR}
                      ["q-filter"]=${FILTERED_DIR} )

# Mapping configuration
# =====================

# map types AND order !
declare -a MAPTYPE=("mouse-rrna" "human-rrna" "mouse-trna"
                    "drosophila-cdna" "mouse-cdna")

# follow the sequential mapping with tophat2 against 'genome'?
# 0 is no, 1 is yes. No other possibilities.
MAP_WITH_TOPHAT2=1

# library extensions - will be used in output files, has to be unique
# attention, uniqueness will NOT be checked!
declare -A LIBEXTS=(  ["mouse-rrna"]="mouse_rRNA"
                      ["human-rrna"]="human_rRNA"
                      ["mouse-trna"]="mouse_tRNA"
                      ["mouse-cdna"]="mouse_cDNA" )

# bowtie2 indexes which should be in "bowtie2 index directory" kept in
# $BOWTIE2_INDEXES environmental variable
declare -A BOWTIE2X=( ["mouse-rrna"]="Mmusculus_v38_rRNA.nr"
                      ["human-rrna"]="human_rRNA-ba"
                      ["mouse-trna"]="Mmusculus_v38_all-tRNA"
                      ["mouse-cdna"]="Mmusculus.GRCm38.84.cdna.ensembl"
                      ["genome"]="Mmusculus.GRCm38.84.dna.ensembl" )

# keep or discard the SAMs
declare -A KEEPSAM=(  ["mouse-rrna"]="no"
                      ["human-rrna"]="no"
                      ["mouse-trna"]="no"
                      ["mouse-cdna"]="yes" )

# use filter_sam to filter SAMs and if so with which options
# possible values are
#   no : do not use filter_sam,
#   passtru : SAM will pass through filter_sam but not filtered at all,
#   basic : unmapped (bit flag 0x4) reads will be filtered out,
#   rrna : special filtering for rRNA mapping,
#   filter : will keep mapped, sense reads with best alignment scores
# If after these commands, a 'k' or 'keep' is added after a seperator (|;-, )
# then filtered out alignments will be kept in a separate BAM.
# examples: "rrna,keep" "filter" "filter-k" "filter;keep" etc
# In all cases, SAM will be converted to BAM with samtools
declare -A TOMERGE=(  ["mouse-cdna"]="filter" )

# sort the BAMs using samtools?
declare -A SORTSAM=(  ["mouse-cdna"]="yes" )

# if BAMs are sorted, keep the unsorted version?
# if not set, it will be interpreted as 'no'!
declare -A KEEPUNSORTED=(  ["mouse-cdna"]="yes" )

# keep the output unmmaped fastq files
# beware this also includes Pre-processing fastq files
declare -A KEEPFASTQ=(  ["trim"]="no"
                        ["s-filter"]="no"
                        ["q-filter"]="yes"
                        ["mouse-rrna"]="no"
                        ["human-rrna"]="no"
                        ["mouse-trna"]="no"
                        ["mouse-cdna"]="yes" )

# output unmapped fastq suffixes
          FILE_EXT+=( ["mouse-rrna"]="non_mouse-rRNA.fastq.gz"
                      ["human-rrna"]="non_rRNA.fastq.gz"
                      ["mouse-trna"]="non_mouse-tRNA.fastq.gz"
                      ["mouse-cdna"]="non_cdna.fastq.gz" )

# which GTF DB to use for -GTF option of tophat2
GTF_DB="/local/databases/mouse/gtf/Mus_musculus.GRCm38.84.sorted.gtf"

# default bowtie2 params
declare -A BOWTIE_PARAMS=( ["default"]="-p 2 -L 15 -k 20 --no-unal -q" )
#                           ["mouse-rrna"]="-p 2 -L 20 -k 5 --no-unal"
#                           ["a map step"]="..." )
#                           -f: fasta, -q: fastq
