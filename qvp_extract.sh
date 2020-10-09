#!/bin/bash
dir=''
#'/home/users/earmlu/ncas_radar_gws_v1/joshua/data/xband/chilbolton/cfradial/calib_v1/sur/'
directorylist=$(find '/home/users/earmlu/ncas_radar_gws_v1/joshua/data/xband/chilbolton/cfradial/calib_v1/sur/20*' -type d)
echo directorylist

for directory in $directorylist
do
  echo ${directory}
  echo ${directory#${dir}}
  #echo python qvp_extraction.py /home/users/earmlu/ncas_radar_gws_v1/joshua/data/xband/chilbolton/cfradial/calib_v1/sur/ 20 ../data/ -s ${directory#${dir}}T000000 -e ${directory#${dir}}T235959
  echo python qvp_extraction.py ../../Insects/Biodar/ 2 ../../Insects/Biodar/ -s ${directory#${dir}}T000000 -e ${directory#${dir}}T235959 -c 134  -a 45,185
  
  #$(basename ${filename})
done
