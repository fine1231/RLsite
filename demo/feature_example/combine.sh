#!/bin/bash
split -l 11 -d -a 3 ../SpaceEnergy+LN/1uud_001.ism ./temp_
for f in temp_*; do
    cat ../SASA/1uud_B/features/ASA/seq_asa.pre | awk 'NR==1 {print $1}' >> "$f"
    cat ../closeness/Outputs/1uud_closeness.txt | awk 'NR==1 {print $1}' >> "$f"
    cat ../degree/Outputs/1uud_degree.txt | awk -F',' '{print $1}' >> "$f"
    cat ../SpaceEnergy+LN/LN.dat | awk 'NR==1 {for (i=1; i<=NF; i++) print $i}' | sed '1d' >> "$f"
    cat ../plmdca/1uud_B/plmDCAout/1uud_B_msa_evo.txt | awk 'NR==1 {for (i=1; i<=NF; i++) print $i}' >> "$f"
    sed -i 's/\r//g' $f
done
paste ./temp_* > 1uud_001.ism
rm ./temp_*
