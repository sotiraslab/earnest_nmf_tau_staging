#!/bin/bash

echo ""
echo "----------"
echo "FIGURE 4"
echo "----------"

echo ""
echo "Running bootstrapping of W-scores..."
SCRIPT="wscore_bootstrap.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running barplots for ADNI..."
SCRIPT="adni_stage_barplots.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."

echo ""
echo "Running barplots for OASIS..."
SCRIPT="oasis_stage_barplots.R"
Rscript -e "source('${SCRIPT}', echo=T)"
echo "Done."