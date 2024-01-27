#!/usr/bin/bash
filename="$1"
while read -r line
do
    name="$line"
    echo "Name read from file - $name"
    cat "../../1KP/gene_trees/filtered/FAA/raxmlboot.$name.filterlen33/RAxML_bipartitions.final" >> tc.tre
done < "$filename"
