#!/bin/bash
###############################################################################
#Program: y.zero.sh
#Programmer: Steve Saleeby
#Description: For use with testing the SCM initial conditions. The SCM.IN data
#  from RAMS for the fields below are already initialized and contain non-zero
#  values. This script replaces those with zero fields so we can test the SCM
#  initialization and see that it matches RAMS.
###############################################################################

#Text file with initial zero data to use
zerotxt="y.zero.txt"

#Number of variables below
numvars="16"

#List of variables to replaces with zeros
var[1]=cn1np.txt
var[2]=cn1mp.txt
var[3]=cn2np.txt
var[4]=cn2mp.txt
var[5]=md1np.txt
var[6]=md1mp.txt
var[7]=md2np.txt
var[8]=md2mp.txt
var[9]=salt_film_np.txt
var[10]=salt_film_mp.txt
var[11]=salt_jet_np.txt
var[12]=salt_jet_mp.txt
var[13]=salt_spum_np.txt
var[14]=salt_spum_mp.txt
var[15]=cifnp.txt
var[16]=dustfrac.txt

#Loop through all variable files in directory "scm.in"
#and replace the onces in the list above with zeroes.
i=1
while [ $i -le $numvars ];
do
 if [ -f scm.in/${var[$i]} ]; then
  cp $zerotxt scm.in/${var[$i]}
 fi
 i=`echo "$i + 1" | bc`
done
