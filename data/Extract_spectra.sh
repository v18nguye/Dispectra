#!/bin/bash

#
# sshsf acces
#

#mkdir -p REJEU_OLERON
#sshfs datarmor:/home/shom_simurep/private/fleckler/OLERON REJEU_OLERON

#
# get point names
#

tmp1=$(ls REJEU_OLERON/WW3_POINTS/ww3.NODE*_200001_spec.nc)
tmp2=$(basename -a $tmp1)
tmp3=${tmp2//_200001_spec.nc/}
pnt_names=${tmp3//ww3./}


#
# get spectra and extract corresponding timeseries from field file
#

STRDATE='201001'

for pnt_name in $pnt_names; do

  spc_file1="REJEU_OLERON/WW3_POINTS/ww3.${pnt_name}_${STRDATE}_spec.nc"
  fld_file1="REJEU_OLERON/WW3_FIELDS//ww3.${STRDATE}.nc"

  spc_file2="WW3/ww3.${pnt_name}_${STRDATE}_spec.nc"
  fld_file2="WW3/ww3.${pnt_name}_${STRDATE}_para.nc"

  tmp1=${pnt_name//NODE/}
  pnt_num=$(echo $tmp1 | sed 's/^0*//')

  echo tmp1=$tmp1
  echo pnt_num=$pnt_num

  # copy spectra AS IS
  cp -v $spc_file1 $spc_file2 
  
  # extract Timeseries from WW3 field file at corresponding node for corresponding date
  ncks --alphabetize -d node,$(($pnt_num-1)),$(($pnt_num-1)),1 -O -o $fld_file2 -v longitude,latitude,pnr,phs0,phs1,phs2,phs3,ptp0,ptp1,ptp2,ptp3,pdir0,pdir1,pdir2,pdir3 $fld_file1

done
