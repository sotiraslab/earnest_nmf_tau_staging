#!/bin/bash

echo ""
echo "----------"
echo "Analysis of non-stageable individuals"
echo "----------"

echo ""
echo "For ADNI..."
Rscript adni_explore_ns.R
echo "Done."

echo ""
echo "For OASIS..."
Rscript oasis_explore_ns.R
echo "Done."