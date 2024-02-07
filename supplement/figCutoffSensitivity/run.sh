#!/bin/bash

echo ""
echo "----------"
echo "Sensitivity analysis for W-scores"
echo "----------"

echo ""
echo "Comparing ordering across W-score thresholds..."
Rscript adni_ordering_sensitivity.R
echo "Done."

echo ""
echo "Comparing staging across W-score thresholds..."
Rscript adni_staging_sensitivity.R
echo "Done."