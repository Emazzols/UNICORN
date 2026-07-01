#!/bin/bash

# Marco Pizzichemi 1.10.2020 marco.pizzichemi@cern.ch

# This script compiles the simulation toolkit and all the support programs/dictionaries
# Run it from main git directory, with
#
# ./deploy.sh [OPTION]
#
# Argument OPTION is optional. If not given, machine compiler is used.
# The only other possible option is lxplus
#
# ./deploy.sh lxplus
#
# which will source LGC97python3 environment.
# If you want to add other options, you need to properly modify the SET ENV VARIABLES section


# exit when any command fails
set -e
# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "Last command \"${last_command}\" returned with exit code $?."' EXIT


### MAKE BUILD FOLDER (IF NOT THERE ALREADY)
mkdir -p build

CMAKE_ARGS=""

### SET ENV VARIABLES
if [ -z "$1" ]
  then
    echo "No argument supplied, so compiling with machine compiler..."
else
  if [ $1 = "lxplus" ]; then
    set --
    echo "Sourcing compiler for lxplus..."
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
    CMAKE_ARGS="--verbose -DCMAKE_CXX_COMPILER=`which g++` -DCMAKE_C_COMPILER=`which gcc`"
    ### old sources, left here for reference
    #source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.06.08/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
    #source /cvmfs/sft.cern.ch/lcg/external/gcc/4.9.1/x86_64-slc6/setup.sh
    #source /cvmfs/geant4.cern.ch/geant4/10.3/x86_64-slc6-gcc49-opt/bin/geant4.sh
  else
    echo "Invalid argument $1 - You can either provide no arg (and the script will use machine compiler), or lxplus, and the script will source the lxplus variables"
    exit 1
  fi
fi

### COMPILE ROOT DICTIONARIES
echo "Compiling ROOT dictionaries..."
rootcling -v4 -f ./include/dictionaries.cxx  -rmf ./build/dictionaries.rootmap -rml ./build/dictionaries.so  ./include/LinkDef.h
cp ./include/dictionaries_rdict.pcm ./build/
cp ./include/struct.hh ./build/
# cp ./include/struct.hh ./

### COMPILE PROGRAMS
#g++ -o ./build/simpleLoop simpleLoop.cpp ./include/dictionaries.cxx `root-config --cflags --glibs`
#g++ -o ./build/simToTree simToTree.cpp ./include/dictionaries.cxx `root-config --cflags --glibs`
# g++ -o ./build/simToTreeComplex simToTreeComplex.cpp ./include/dictionaries.cxx `root-config --cflags --glibs`
g++ -o ./build/extractConfiguration extractConfiguration.cpp `root-config --cflags --glibs` -lHist -lCore -lMathCore
g++ -o ./build/extractSimToTree extractSimToTree.cpp `root-config --cflags --glibs` -lHist -lCore -lMathCore

### COMPILE SIMULATION TOOLKIT
echo "Compiling g4matrix..."
cd build
echo "Running cmake command as: cmake $CMAKE_ARGS ../"
cmake $CMAKE_ARGS ../
make -j8
cd ..
echo "Done."
