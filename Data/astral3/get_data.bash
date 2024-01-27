#!/usr/bin/bash
filename="$1"
while read -r line
do
    name="$line"
    echo "Name read from file - $name"
    cat "../../1KP/gene_trees/unfiltered/FNA2AA/raxmlboot.$name/RAxML_bipartitions.final" >> tc.tre
done < "$filename"
