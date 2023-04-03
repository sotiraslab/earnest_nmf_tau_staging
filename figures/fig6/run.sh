#!/bin/bash

echo "------------------"
echo "TABLE 2"
echo "------------------"

echo ""
echo "WARNING: This takes a relatively long time to run!"
echo "         Especially on first time when expression data is download"
echo ""

echo ""
echo "Getting abagen gene expression data..."
python get_abagen_expression_dkt.py
echo "Done."

echo ""
echo "Getting PTCs ready for spatial correlation comparison..."
python prep_ptcs.py
echo "Done."

echo ""
echo "Running spatial correlation with neuromaps..."
python run_wightman_correlations.py
echo "Done."

echo ""
echo "Collecting significant results..."
python collect_results.py
echo "Done."

echo ""
echo "Creating gene expression brain figures..."
python gene_expression_figures.py
echo "Done."

echo "Create association table image..."
python gene_table.py
echo "Done."