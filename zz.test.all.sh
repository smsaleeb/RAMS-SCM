#!/bin/bash

version="scm"
indir="AA.scm.in.out.data"

#Make fresh I/O directories for SCM
if [ -d scm.in  ]; then rm -f scm.in/*; fi
if [ -d scm.out ]; then rm -f scm.out/*; fi
if [ ! -d scm.in   ]; then mkdir scm.in; fi
if [ ! -d scm.out  ]; then mkdir scm.out; fi
outf="z.out.diffs.txt"
echo "DOING DIFFS" > $outf

#******************************************************************
#Testing DUST 3 later timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.dust3 ]; then rm -fr scm.out.dust3; fi
cp $indir/RAMSIN.full.dust3 RAMSIN
cp $indir/scm.in.dust3/* scm.in
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.dust3" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.dust3 >> $outf
cp -r scm.out scm.out.dust3

#******************************************************************
#Testing DUST 3 initialization and 1st timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.dust3.0 ]; then rm -fr scm.out.dust3.0; fi
cp $indir/RAMSIN.full.dust3 RAMSIN
cp $indir/scm.in.dust3.0/* scm.in
y.zero.sh
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.dust3.0" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.dust3.0 >> $outf
cp -r scm.out scm.out.dust3.0

#******************************************************************
#Testing SPL CLN later timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.spl.cln ]; then rm -fr scm.out.spl.cln; fi
cp $indir/RAMSIN.full.spl.cln RAMSIN
cp $indir/scm.in.spl.cln/* scm.in
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.spl.cln" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.spl.cln >> $outf
cp -r scm.out scm.out.spl.cln

#******************************************************************
#Testing SPL CLN initialization and 1st timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.spl.cln.0 ]; then rm -fr scm.out.spl.cln.0; fi
cp $indir/RAMSIN.full.spl.cln RAMSIN
cp $indir/scm.in.spl.cln.0/* scm.in
y.zero.sh
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.spl.cln.0" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.spl.cln.0 >> $outf
cp -r scm.out scm.out.spl.cln.0

#******************************************************************
#Testing Supercell with lofted dust later timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.supercell.aeroloft ]; then rm -fr scm.out.supercell.aeroloft; fi
cp $indir/RAMSIN.supercell.aeroloft RAMSIN
cp $indir/scm.in.supercell.aeroloft/* scm.in
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.supercell.aeroloft" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.supercell.aeroloft >> $outf
cp -r scm.out scm.out.supercell.aeroloft

#******************************************************************
#Testing Supercell with lofted dust off timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.supercell.aeroloft.dtrad ]; then rm -fr scm.out.supercell.aeroloft.dtrad; fi
cp $indir/RAMSIN.supercell.aeroloft.dtrad RAMSIN
cp $indir/scm.in.supercell.aeroloft.dtrad/* scm.in
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.supercell.aeroloft.dtrad" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.supercell.aeroloft.dtrad >> $outf
cp -r scm.out scm.out.supercell.aeroloft.dtrad

#******************************************************************
#Testing Supercell with dust profile later timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.supercell.aeroprof ]; then rm -fr scm.out.supercell.aeroprof; fi
cp $indir/RAMSIN.supercell.aeroprof RAMSIN
cp $indir/scm.in.supercell.aeroprof/* scm.in
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.supercell.aeroprof" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.supercell.aeroprof >> $outf
cp -r scm.out scm.out.supercell.aeroprof

#******************************************************************
#Testing Supercell with dust profile initialization and 1st timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.supercell.aeroprof.0 ]; then rm -fr scm.out.supercell.aeroprof.0; fi
cp $indir/RAMSIN.supercell.aeroprof RAMSIN
cp $indir/scm.in.supercell.aeroprof.0/* scm.in
y.zero.sh
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.supercell.aeroprof.0" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.supercell.aeroprof.0 >> $outf
cp -r scm.out scm.out.supercell.aeroprof.0

#******************************************************************
#Testing Supercell with absorbing carbon profile later timestep
if [ -f scm.in -o -f scm.out ]; then rm -f scm.in/* scm.out/*; fi
if [ -d scm.out.supercell.aeroabsc ]; then rm -fr scm.out.supercell.aeroabsc; fi
cp $indir/RAMSIN.supercell.aeroabsc RAMSIN
cp $indir/scm.in.supercell.aeroabsc/* scm.in
rams-$version
echo "------------------------------------------------------" >> $outf
echo "Diffs scm.out.supercell.aeroabsc" >> $outf
y.diff.sh 2 scm.out $indir/scm.out.supercell.aeroabsc >> $outf
cp -r scm.out scm.out.supercell.aeroabsc

#******************************************************************
#Final removal of temporary scm.in scm.out directories
if [ -d scm.in -o -d scm.out ]; then rm -fr scm.in scm.out; fi
if [ -f RAMSIN ]; then rm -f RAMSIN; fi

###############################################################################
###############################################################################

#Delete all the output directories. Could comment this out if you need to
#examine the output vertical profiles for troubleshooting.

rm -fr scm.out.dust* scm.out.spl* scm.out.supercell*
