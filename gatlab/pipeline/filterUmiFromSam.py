#!/usr/bin/env python3

"""
filterUmiFromSam

Dedups SAM output of UMI-tools group/count tool
"""

######################################################################
### This script filere mapped reads that were grouped using the
### UMI-tools count function. It deals with the fact that reads
### normally map to multiple locations (eg. transcripts) and it
### retains always the same read per group. This maintains the
### multiple mapping information that might be lost if a read is
### randombly choosen for each group.

from __future__ import print_function
import sys
import argparse
import re

__author__ = "René Dreos"
__copyright__ = "Copyright 2019, René Dreos"
__license__ = "GPLv3"
__version__ = "1.0.0"
__maintainer__ = "René Dreos"
__email__ = ""
__status__ = "Development"


def parserFunction():
    # parse the command line arguments:
    parser = argparse.ArgumentParser(description='This script filters reads '
                                     'in a SAM file (from STDIN) that belong '
                                     'to the same group as defined by '
                                     'UMI-tools making sure to keep the same '
                                     'read representative for each group. Note'
                                     'that SAM file must be sorted by read ID.')
    parser.add_argument('-d','--debug',
                        help='Produce a verbose output, useful for debugging '
                        'purposes.',
                        action="store_true")
    return parser.parse_args()

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def compress(string):
    ids = string.split( "_" )
    name = ids[0].split( ':' )
    result = name[4] + ':' +  name[5] + ':' + name[6]
    return result

def main():
    args = parserFunction()
    dbg = args.debug
    
    gReps = dict()      # Link group ID with read representative
    readHash = dict()   # Hash that store read IDs with new ID
    #gIDs = dict()       # Might be deleted?
    readStatus = dict() # 1: read is a representative; 0: read is PCR
                        # artefact
    readGroups = dict() # Is the group of the actual line
    count = dict()      # Counter container
    allGroups = dict()  # All groups that a read belong to
    allReads = dict()   # All reads that belong to a group 
    alignment = ''      # Is the SAM line
    oldReadId = ''      # Id of previous line
    oldGroupId = ''     # Group ID of previous line
    lines = 1           # Number of line read
 
    count[ "_total_aligments" ] = 0
    count[ "_total_reads" ] = 0
    count[ "_total_groups" ] = 0
    count[ "_printed_alignments" ] = 0
    count[ "_singletons_multi" ] = 0
    count[ "_branched_groups" ] = 0
    
    for inline in sys.stdin:
        inline = inline.rstrip()
        # check if read is aligned
        if inline.startswith( '@' ):
            print(inline)
        else:
            count[ "_total_aligments" ] += 1
            lineFields = inline.split( "\t" )
            #newReadId = compress(lineFields[0])
            #readHash[newReadId] = lineFields[0]
            readId = lineFields[0]
            mGroup = re.search('UG:i:([0-9]+)', inline)
            if mGroup:
                readGroup = int(mGroup.group(1))
            else:
                eprint('Read', readId,
                       'does not belong to any group. Quitting',
                       sep=' ')
                break
            ##if dbg: print(readId, ':', readGroup, sep=' ')

            if readId not in allGroups:
                allGroups[readId] = [readGroup]
            else:
                allGroups[readId].append(readGroup)
                
            if readGroup not in allReads:
                allReads[readGroup] = [readId]
            else:
                allReads[readGroup].append(readId)

            #if readGroup not in gIDs: gIDs[readGroup] = True

            ## When group ID change (first line of a new group):
            if readGroup != oldGroupId:
                if oldGroupId != '':
                    count[ "_total_groups" ] += 1
                    count[ "_printed_alignments" ] += 1
                    # there are cases that a group of reads that were
                    # not representatives in any group are now
                    # togheter. Give them a new representative
                    if oldGroupId not in gReps:
                        gReps[oldGroupId] = oldReadId
                        globalGroups = set()
                        for read in allReads[oldGroupId]:
                            globalGroups = globalGroups.union(allGroups[read])
                            if dbg: print(read, 'global groups', globalGroups)
                        firstGroup = sorted(globalGroups)[0]
                        firstRep = gReps[firstGroup]
                        # add this group to the representative
                        allGroups[firstRep].append(oldGroupId)
                        ## alignment = oldLine.replace(oldReadId, firstRep)
                        alignment = oldLine.replace(oldReadId,
                                                    firstRep)
                        # check if all reads in the group belong to
                        # the same groups as the representative
                        checkReads = True
                        counter = 0
                        for read in allReads[oldGroupId]:
                            groupInt = set(allGroups[read]).intersection(set(allGroups[firstRep]))
                            if groupInt != set(allGroups[read]):
                                counter += 1

                        if counter > 1:
                            checkReads = False
                                
                        if dbg:
                            print(oldGroupId, 'did not have a representative',
                                  'taking the last read id:', oldReadId,
                                  sep=' ')
                            print(oldReadId, '->', firstGroup, '->', firstRep)
                            if not checkReads:
                                print('WARNING: this is a branching group')

                        if checkReads:
                            count[ "_singletons_multi" ] += 1
                        else:
                            count[ "_branched_groups" ] += 1
                            
                    if not dbg:
                        print(alignment)

                # if this read was already in a group:
                if readId in readGroups:
                    # check if this read was a representative of that
                    # group
                    if gReps[readGroups[readId]] == readId:
                    #if readStatus[readId]:
                        if dbg:
                            print(readId, 'REP1',
                                  allGroups[readId], sep='\t')
                        gReps[readGroup] = readId
                        readStatus[readId] = True
                        alignment = inline
                    else:
                        readStatus[readId] = False
                        if dbg:
                            print(readId, 'PCR',
                                  allGroups[readId], sep='\t')
                else: # (if this read was not in a group make the REP
                      # of this group)
                    if dbg:
                        print(readId, 'REP2',
                              allGroups[readId], sep='\t')
                    gReps[readGroup] = readId
                    readStatus[readId] = True
                    alignment = inline
                    
                # readGroups[readId] = readGroup
                # oldGroupId = readGroup
                
            else: # if I am inside the group members:
                # if the read was a REP in another group:
                if readId in readStatus and readStatus[readId]:
                    if dbg: print(readId, 'REP3',
                                  allGroups[readId], sep='\t')
                    alignment = inline
                    ## impose this read as representative (clean the
                    ## readStatus dictionary)
                    if readGroup in gReps:
                        readStatus[gReps[readGroup]] = False
                    readStatus[readId] = True
                    gReps[readGroup] = readId
                    readGroups[readId] = readGroup
                else: # (the read was not a representative, then it is
                      # a PCR)
                    readStatus[readId] = False
                    readGroups[readId] = readGroup
                    # check if read group does not have a rep, if the
                    # read was present in another group and if this
                    # has a representative
                    if readGroup not in gReps:
                        for group in allGroups[readId]:
                            if group in gReps:
                                firstRep = gReps[group]
                                alignment = inline.replace(readId,
                                                           firstRep)
                                count[ "_singletons_multi" ] += 1
                                gReps[readGroup] = firstRep
                                break
                    
                    if dbg: print(readId, 'PCR',
                                  allGroups[readId], sep='\t')
                    
            oldLine = inline
            oldReadId = readId
            readGroups[readId] = readGroup
            oldGroupId = readGroup
        
    # Print summary:
    print("Total alignments:              ",
          count[ "_total_aligments" ], file=sys.stderr)
    print("Total reads:                   ",
          len(readStatus), file=sys.stderr)
    print("Total unique reads:            ",
          sum(readStatus.values()), file=sys.stderr)
    print("Average duplication level:     ",
          round(len(readStatus)/sum(readStatus.values()), 3),
          file=sys.stderr)    
    print("Merged reads groups:           ",
          count[ "_singletons_multi" ], file=sys.stderr)
    print("Branched reads groups:         ",
          count[ "_branched_groups" ], file=sys.stderr)
    print("Total deduplicated alignments: ",
          count[ "_total_groups" ], file=sys.stderr)
    print("Printed alignements:           ",
          count[ "_printed_alignments" ], file=sys.stderr)
