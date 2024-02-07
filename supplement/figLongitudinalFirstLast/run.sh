#!/bin/bash

echo ""
echo "----------"
echo "Alternative longitudinal transition analysis"
echo "----------"

echo ""
echo "Determining longitudinal FTP staging for ADNI..."
python longitudinal_staging_adni.py
echo "Done."

echo ""
echo "Creating longitudinal staging heatmap..."
python longitudinal_staging_heatmap.py
echo "Done."

echo ""
echo "Running permutation test for longitudinal staging..."
python longitudinal_permutation_test.py
echo "Done."