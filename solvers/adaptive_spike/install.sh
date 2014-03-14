#!/bin/sh

#*******************************************************************************
#                              INTEL CONFIDENTIAL
#   Copyright(C) 2008 Intel Corporation. All Rights Reserved.
#   The source code contained  or  described herein and all documents related to
#   the source code ("Material") are owned by Intel Corporation or its suppliers
#   or licensors.  Title to the  Material remains with  Intel Corporation or its
#   suppliers and licensors. The Material contains trade secrets and proprietary
#   and  confidential  information of  Intel or its suppliers and licensors. The
#   Material  is  protected  by  worldwide  copyright  and trade secret laws and
#   treaty  provisions. No part of the Material may be used, copied, reproduced,
#   modified, published, uploaded, posted, transmitted, distributed or disclosed
#   in any way without Intel's prior express written permission.
#   No license  under any  patent, copyright, trade secret or other intellectual
#   property right is granted to or conferred upon you by disclosure or delivery
#   of the Materials,  either expressly, by implication, inducement, estoppel or
#   otherwise.  Any  license  under  such  intellectual property  rights must be
#   express and approved by Intel in writing.
#*******************************************************************************

# Installation script of Intel(R) Adaptive Spike-Based Solver

while true; do
  clear
  echo "Welcome to the Intel(R) Adaptive Spike-Based Solver Installation"
  echo
  echo "Please make your selection by entering an option from the choices below:"
  echo
  echo "	1. Install"
  echo "	2. Installation Guide"
  echo "	h. Help"
  echo "	x. Exit"
  echo
  read -p "Please type a selection:  " choice
  case $choice in
    1)
      # Ask the user to agree to EULA before proceeding installation
      
      echo "--------------------------------------------------------------------------------"
      echo "Please read the following license agreement carefully.  Prior to installing the"
      echo "software, you will be asked to agree to the terms and conditions of the"
      echo "following license agreement."
      echo "--------------------------------------------------------------------------------"
      echo "Please press Enter to continue."
      read  dummy &> /dev/null

      more spikeEULA.txt
      while true; do
        echo
        echo "Enter 'accept' to continue, 'reject' to return to the main menu"
        read ans &> /dev/null
        case $ans in
          [Aa][Cc][Cc][Ee][Pp][Tt])
            break 2;;
          [Rr][Ee][Jj][Ee][Cc][Tt])
            continue 2;;
        esac 
      done 
      break
      ;;
    2)
      more Install.txt
      echo -n "Press <Enter> to continue..."
      read  dummy &> /dev/null
      ;;
    h)
      echo "---------------------------------------------------------------------------"
        echo "Installation Options"
        echo " "
        echo "-- Option 1 --"
        echo "Select this option to proceed with the software installation. The installation will prompt you for an installation directory."
        echo " "
        echo "-- Option 2 --"
        echo "Select this option to read the Installation Guide."
        echo " "
        echo "-- Option h --"
        echo "Selecting option h displays this help message."
        echo " "
        echo "-- Option x --"
        echo "Selecting option x exits the software installation."
        echo "---------------------------------------------------------------------------"
        echo -n "Press <Enter> to continue..."
        read  dummy &> /dev/null
        ;;
    x)
      exit 0
      ;;
  esac
done


# Get the installation directory.  Default is $HOME/SPIKE/1.0

while true; do
  install_dir=$HOME/SPIKE/1.0
  read -p "Where do you want to install to?  Specify directory starting with '/'. [$install_dir]: " input_dir
  
  if [ ! -z $input_dir ]; then
    install_dir=$input_dir
  fi

  echo "Chosen: $install_dir"
  if echo $install_dir | grep -q -s ^/ ; then
    if [ -d $install_dir ]; then
      echo "WARNING: Destination directory is already exist."
      while true; do
        ans="No"
        read -p "Continuation can lead to loss of stored data. Would you like to overwrite this directory? ( Yes/No ) [ $ans ]: " ans
        case $ans in
          [Yy][Ee][Ss])
            break 2;;
          [Nn][Oo])
            continue 2;;
        esac
      done
    else
      mkdir -p $install_dir
      if [ $? == 0 ]; then
        break
      else
        echo "WARNING: Cannot create installation directory."
      fi
    fi
  fi
done  
# The installation directory is ready to use.
echo "Installing package..."

tar xvfz data.tar.gz -C $install_dir

# Use the initialization script template to generate the 
# initialization scripts with SPIKEROOT set to the installation directory
for scfile in $install_dir/tools/environment/*; do
  outfile=$install_dir/tools/environment/`basename $scfile .template`  
  sed 's:SPIKEROOT=.*:SPIKEROOT='"$install_dir:" < $scfile > $outfile
  rm $scfile
done

if [ $? = "0" ]; then
  echo "Intel(R) Adaptive Spike-Based Solver was successfully installed."
fi
