#!/bin/bash

echo ""
echo "----------"
echo "LONGITUDINAL STAGING"
echo "----------"

echo ""
echo "Running survival plots for ADNI..."
SCRIPT="adni_cdr_survival.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running survival plots for OASIS..."
SCRIPT="oasis_cdr_survival.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

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