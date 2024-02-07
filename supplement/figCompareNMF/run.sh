#!/bin/bash

echo ""
echo "----------"
echo "8 PTCS"
echo "----------"

echo ""
echo "Creating ADNI/OASIS figure comparison..."
python create_8ptc_figures.py
echo "Done."

echo ""
echo "Creating winner take all assignment figures..."
python create_winner_take_all_brains.py
echo "Done."

echo ""
echo "Running spin test comparison..."
python spintest.py
echo "Done."
