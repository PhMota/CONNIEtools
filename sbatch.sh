#!/bin/bash

export CONNIE_INPUT
files=$(ls -tr /share/storage2/connie/nu_processing/scripts/ProcCat/cut_scn_osi_raw_gain_catalog_data_*[^skim1].root)
for i in $files
do
    CONNIE_INPUT=$i
    echo $CONNIE_INPUT
    sbatch script.slurm
done
exit