#!/bin/bash

echo ""
echo "----------"
echo "CDR Survival with NS"
echo "----------"

SCRIPT="adni_cdr_survival_includeNS.R"
Rscript "${SCRIPT}"

SCRIPT="oasis_cdr_survival_includeNS.R"
Rscript "${SCRIPT}"