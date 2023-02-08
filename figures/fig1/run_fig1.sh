#!/bin/bash

echo ""
echo "Running reconstruction error plots..."
python compare_recon_error.py
echo "Done."

echo ""
echo "Running intrasample reproducibility plots..."
python intrasample_reproducibility.py
echo "Done."

echo ""
echo "Running intrasample reproducibility plots..."
python intersample_reproducibility.py
echo "Done."
