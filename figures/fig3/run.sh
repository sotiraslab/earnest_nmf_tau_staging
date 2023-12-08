#!/bin/bash

echo ""
echo "----------"
echo "CROSS-SECTIONAL STAGING"
echo "----------"

echo ""
echo "Running bootstrapping of W-scores..."
SCRIPT="wscore_bootstrap.R"
Rscript "${SCRIPT}"
echo "Done."

echo ""
echo "Running barplots for ADNI..."
SCRIPT="adni_stage_associations.R"
Rscript "${SCRIPT}"
echo "Done."

echo ""
echo "Running barplots for OASIS..."
SCRIPT="oasis_stage_associations.R"
Rscript "${SCRIPT}"
echo "Done."

echo ""
echo "Running stage brain map..."
python stage_map.py
echo ''