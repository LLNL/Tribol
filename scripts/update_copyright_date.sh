#!/usr/bin/env bash

# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC.

#=============================================================================
# Attempt to update  the copyright date in all source-controlled files 
# that contain the text: "Lawrence Livermore National Security, LLC.".
#
# IMPORTANT: Since this file is not modified (it is running the shell
# script commands), you must EDIT THE COPYRIGHT DATES ABOVE MANUALLY.
#
# Edit the 'find' command below to change the set of files that will be
# modified.
#
# Change the 'sed' command below to change the content that is changed
# in each file and what it is changed to.
#
#=============================================================================
#
# If you need to modify this script, you may want to run each of these
# commands individual from the command line to make sure things are doing
# what you think they should be doing. This is why they are separated into
# steps here.
#
#=============================================================================
#
# Note: The following files do not fit this pattern and require hand editing
# (or integration into this script):
#
#   < add any such files here >
#=============================================================================

#=============================================================================
# First find all the files we want to modify
#=============================================================================
git grep -l "Lawrence Livermore National Security, LLC." | grep -v update_copyright > files2change

#=============================================================================
# Replace the old copyright dates with new dates
#=============================================================================
for i in `cat files2change`
do
    echo $i
    cp $i $i.sed.bak
    sed "s/Copyright (c) 2017-2022/Copyright (c) 2017-2023/" $i.sed.bak > $i
done

#=============================================================================
# Remove temporary files created in the process
#=============================================================================
find . -name \*.sed.bak -exec rm {} \;
rm files2change
