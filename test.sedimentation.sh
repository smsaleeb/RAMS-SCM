#!/bin/bash
RAMSNAMELIST='test.sedimentation.RAMSIN'
scmin='scm.sedimentation-test.in'
scmout='scm.sedimentation-test.out'

#Make fresh I/O directories for SCM
if [ -d scm.in  ]; then rm -f scm.in/*; fi
if [ -d scm.out ]; then rm -f scm.out/*; fi
if [ ! -d scm.in  ]; then mkdir scm.in; fi
if [ ! -d scm.out ]; then mkdir scm.out; fi

# Initial state set in file "mic_init_scm.f90" subroutine "init_custom"

# Copy specific RAMSIN and input files into default "RAMSIN" and "scm.in"
cp $RAMSNAMELIST RAMSIN
cp $scmin/* scm.in
# Run the SCM
rams-scm
# Make the specific output directory
if [ ! -d $scmout ]; then
 mkdir $scmout
fi
# Move files from default directories to test specific ones
mv scm.out/* $scmout
rm -f scm.in/* scm.out/*
rmdir scm.in scm.out
rm -f RAMSIN

#Make some plots with python if system is set up to do so. Change as needed here.
conda activate General
python3 test.PlotVerticalProfiles.py
test.mp4-from-images.sh 7 scm.sedimentation-test.out
