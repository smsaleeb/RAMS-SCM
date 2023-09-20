#!/bin/bash
###################################################################################################
#Using ffmpeg installed into python environment.
#Harded coded to look for first file like "000.png". The "000" is the
#first file indicator. The type of image is set below in "FILE_SUF".
###################################################################################################

#Path to ffmpeg
ffmpeg="/home/smsaleeb/miniconda3/envs/General/bin/ffmpeg"
#Image file type suffic
FILE_SUF="png"

#Input directory path for analysis directories and prefix
if [ -n "$1" -a -n "$2" ]; then
 frate=$1 # Input framerate
 fdir=$2  # Input directory path to list of files to animate
 echo ""
 echo "Making Animations Now..."
else
 echo "Input arguments with frame rate and directory path to images to create mp4 video"
 echo "EX. Pybash.mp4-from-images.sh 15 ./PNGS/z.XY.single_filled_contours"
 exit
fi

#LOOPING THRU VARIOUS PLOT DIRECTORIES
if [ -d $fdir ]; then
  cd $fdir
  #LOOPING THRU LIST OF PLOTS TYPES
  for f in *.t000.$FILE_SUF
  do
    #Do this if the initial 000 file exists
    if [ -f $f ]; then
     #Get base file names
     file1=$(basename $f)
     imagename=${file1: : -9}

     #DETERMINE IF ANIMATED IMAGES ALREADY MADE
     animname="anim."$imagename".mp4"
     animgif="anim."$imagename".gif"
     if [ -f $animname ]; then rm -f $animname; fi
     if [ -f $animgif ];  then rm -f $animgif;  fi
     echo "Making new animation: $fdir --- $animname"
     imagenamestart=$imagename".t000."$FILE_SUF
     #echo "  $imagename  :  $imagenamestart"

     #Find last file in series
     lastim=-1
     for fl in $imagename*
     do
       lastim=`echo "$lastim + 1" | bc`
     done
     if [ $lastim -lt 100 ]; then lastim="0"$lastim; fi
     if [ $lastim -lt 10  ]; then lastim="0"$lastim; fi
     #echo "Last time: $lastim"

     #START THE ANIMATION MAKING IF THE FIRST FILE EXISTS
     if [ -f $imagenamestart ]; then
       echo "File: $imagename.*.$FILE_SUF"

       #Make MP4 video
       ############################################################################################
       # 1. "-pix_fmt yuv420p" will ensure a chroma subsampling that is compatible for all players
       # 2. "-crf 0" will give best quality and can be from 0 to 51
       # 3. "-b 2500k" is the bitrate of the movie

       $ffmpeg -framerate $frate -pattern_type glob -i "$imagename.*.$FILE_SUF" \
        -pix_fmt yuv420p -b 2500k -vcodec libx264 \
        -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" $animname 1>/dev/null 2>&1

     fi #endif check for file
    fi #endif check for initial file
  done
fi #endif check for directory
echo " "
echo "FINISHED"
