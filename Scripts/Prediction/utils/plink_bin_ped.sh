#!/bin/bash

ml PLINK/1.07-x86_64

# Args
while getopts ":d:n:" opt; do
    case $opt in
    d)
        ped_map_dir=$OPTARG;;
    n)
        ped_map_name=$OPTARG;;
    esac
done

cmd="cd ${ped_map_dir}; plink --file ${ped_map_name} --out ${ped_map_name} --make-bed --all --noweb"
# echo "${cmd}"
eval "${cmd}"
