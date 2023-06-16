#!/bin/bash
# Takes the output from cluster identifier (from Scramble) and creates a fasta file with the read information and consensus sequence of clipped reads
# Example usage:
#	ConvertClusterIdentifierOutputToFasta.sh cluster_identifier_output.txt > cluster_identifier_output.fasta

input_cluster_data=${1}
output_clip_consensus=$(echo $input_cluster_data | rev | cut -d '.' -f 2- | rev)".clipConsensus.fasta"
output_anchor_consensus=$(echo $input_cluster_data | rev | cut -d '.' -f 2- | rev)".anchoredConsensus.fasta"

echo ${output_clip_consensus}
echo ${output_anchor_consensus}
echo ${input_cluster_data}

awk -F'\t' '{print ">"$1"_"$2"_"$3"_ClipConsensus"; print $4}' $1 > ${output_clip_consensus}
awk -F'\t' '{print ">"$1"_"$2"_"$3"_AnchoredReadSeq"; print $5}' $1 > ${output_anchor_consensus}
