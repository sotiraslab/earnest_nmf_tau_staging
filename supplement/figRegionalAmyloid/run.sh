#!/bin/bash

echo ""
echo "----------"
echo "Interregional tau & amyloid"
echo "----------"

echo ""
echo "For ADNI..."
Rscript adni_amyloid_regression.R
echo "Done."

echo ""
echo "For OASIS..."
Rscript oasis_amyloid_regression.R
echo "Done."
