##########################################################
# PRINT
##########################################################

echo ""
echo "--------------------"
echo "NMF-tau"
echo "--------------------"
echo ""
echo "Author: Tom Earnest"
echo "Date: February 2023"

##########################################################
# PARSE ARGUMENTS
##########################################################

usage() {
	echo "Usage:	run.sh [options]"
	echo ""
	echo "Description:"
	echo "    Run the dataset creation for the NMF-tau project."
	echo ""
	}

CREATE_FIGURES="1"

while getopts hn arg
do
	case $arg in
	h)	usage
		exit 0;;
	?)	echo ""
		echo "Unknown arguments passed; exiting."
		echo ""
		usage;
		exit 1;;
	esac
done

##########################################################
# SET VARIABLES
##########################################################

MAIN_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

LOOKSCRIPT='run_build_datasets.sh'
if [[ ! -f ${MAIN_DIR}/${LOOKSCRIPT} ]]
then
	# default to trying current working directory
	echo ""
	echo "Warning: Cannot find main '${LOOKSCRIPT}' script at ${MAIN_DIR}."
	echo "Setting to current working directory (PWD): ${PWD}"
	MAIN_DIR=$PWD
fi

##########################################################
# CHECK NECESSARY FILES
##########################################################

RAWDATA=${MAIN_DIR}/'rawdata'

NEEDED='oasis_amyloid.csv oasis_flortaucipir.csv OASIS3_data_files PTDEMOG.csv'

echo ""
echo "Check for raw data files..."
echo "---------------"
echo ""

for item in $NEEDED
do
	path=${RAWDATA}/${item}
	if [[ ! -e $path ]]
	then
		echo ""
		echo "Error; cannot find at least one required raw data files: ${item}.  Check the readme."
		echo "Exiting."
		echo ""
		exit 1
    else
        echo "Found required raw data: ${item}"
	fi
done

echo ""
echo "Success!"

##########################################################
# CALL R SCRIPT HELPER
##########################################################

function callR {
	SCRIPT=$1
	COMMAND=(Rscript -e "source('${SCRIPT}', echo=T)")

	echo ""
	echo "Rscript command:"
	echo "    ${COMMAND[@]}"

	echo ""
	echo "********* COMMAND OUTPUT START *********"
	"${COMMAND[@]}"
	echo "********* COMMAND OUTPUT END   **********"	
}	

##########################################################
# ADNI
##########################################################

PREP=${MAIN_DIR}/'prep'
DERIVATIVES=${MAIN_DIR}/'derivatives'

PREP_ADNI=${PREP}/'adni'
DERIVATIVES_ADNI=${DERIVATIVES}/'adni'

PREP_OASIS=${PREP}/'oasis3'
DERIVATIVES_OASIS=${DERIVATIVES}/'oasis3'

echo
echo "ADNI"
echo "---------------"

mkdir -p ${DERIVATIVES_ADNI}

# 1. Create main dataframe
SCRIPT=${PREP_ADNI}/create_main_dataframe.R
echo ""
echo "1. Creating main ADNI dataset..."
callR $SCRIPT


# 2. Create NMF matrices
SCRIPT=${PREP_ADNI}/create_nmf_matrix.R
echo ""
echo "2. Creating NMF inputs..."
callR $SCRIPT

# 3. Project NMF
SCRIPT=${PREP_ADNI}/project_adni_8ptc.R
echo ""
echo "3. Getting projection onto NMF components..."
callR $SCRIPT

echo ""
echo "Preparing ADNI dataset, complete!"