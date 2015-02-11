#!/bin/bash
#Script to run the static code analysis for the OpenMSI Toolkit

OUTPUT_DIR=$PWD"/_results"
STATISTICS_DIR=$PWD"/_statistics"
if [ ! -d "$OUTPUT_DIR" ]; then
  echo "Creating output directory: " $OUTPUT_DIR
  mkdir $OUTPUT_DIR
fi
if [ ! -d "$STATISTICS_DIR" ]; then
  echo "Creating statistics directory: " $STATISTICS_DIR
  mkdir $STATISTICS_DIR
fi


echo "Removing previous results"
rm $OUTPUT_DIR/pylint*.txt

echo "Executing pylint:"
pylint --rcfile=openmsi_tk_pylint_settings.rc ../../omsi/

echo "Moving pylint output to results directory"
mv pylint*.txt $OUTPUT_DIR

echo "Recording statistics"
cp $OUTPUT_DIR/pylint_global.txt $STATISTICS_DIR/$(date +'%m-%d-%Y__%T')_openmsi_tk_code_stats.txt
