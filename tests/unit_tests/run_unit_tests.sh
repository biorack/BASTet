#!/bin/bash
#Script to run the static code analysis for the OpenMSI Toolkit

# 1) Set the PYTHONPATH so that we can find the openmsi-tk installation
export PYTHONPATH=$PWD/../../
echo "Set PYTHONPATH to: " $PYTHONPATH

# 2) Execute all unit tests
python -m unittest discover -v omsi_tests/
