#!/bin/bash

echo ""
echo "----------"
echo "Boostrap staging for OASIS"
echo "----------"

SCRIPT="wscore_bootstrap_OASIS.R"
Rscript -e "source('${SCRIPT}', echo=T)"