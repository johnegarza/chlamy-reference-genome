#!/usr/bin/env bash

#SBATCH -n 1
#SBATCH --mem=16000

python end_indexer.py fosmid_pairs.txt ../full_node_list.txt > indexed_ends.txt
