#!/bin/bash

echo ""
echo "----------"
echo "Analysis of neuropsychological domains"
echo "----------"

echo ""
echo "For ADNI..."
Rscript adni_neuropsych_associations.R
echo "Done."

echo ""
echo "For OASIS..."
Rscript oasis_neuropsych_associations.R
echo "Done."
