# Learning flortaucipir patterns with NMF

This respository contains bublic code for Tom's project (applying non-negative matrix factroization to flortaucipir images).  It is currently a work in progress - I am moving some of the messier code and removing ADNI/OASIS3 data.

## Downloading raw data

This repository does not include data from ADNI/OASIS which requires approval to access.  This section will document the instructions for accessing the data necessary for this project.

**All files described here (for both ADNI & OASIS) should be placed in the `rawdata` folder.**   Note that some files will need to be renamed.

### ADNI

- **ADNIMERGE R package**
  - https://ida.loni.usc.edu/ > ADNI > Download > Study Data > Study Info > Data & Database > ADNIMERGE - Key ADNI tables merged into one table - Packages for R [ADNI1,GO,2]
  - Follow instructions for installation here: https://adni.bitbucket.io/
  - Data for this project were download on November 23, 2022
- **PTDEMOG.csv**
  - https://ida.loni.usc.edu/ > ADNI > Download > Study Data > Subject Characteristics > Subject Demographics > Subject Demographics [ADNI1,GO,2,3]
  - This file was specifically downloaded to find some demographic information missing from ADNIMERGE R.

### OASIS-3

- **Amyloid PUP results**
  - https://central.xnat.org/ > Browse > My Projects > OASIS3 > Add Tab > PUPs > Options > Edit Columns > Add column "PET_fSUVR_rsf_TOT_CORTMEAN" > Options > Spreadsheet
  - **RENAME THIS FILE:** `oasis_amyloid.csv`
- **Flortaucipir PUP results**
  - https://central.xnat.org/ > Browse > My Projects > OASIS3_AV1451 > Add Tab > PUPs > Options > Edit Columns > Add all columns containing "fSUVR" (the ones containing "rsf" are not needed, but can be included) > Options > Spreadsheet
  - **RENAME THIS FILE:** `oasis_flortaucipir.csv`
- **OASIS3 Data Files**
  - https://central.xnat.org/ > Browse > My Projects > OASIS3 > On the "Subjects" tab, click "0AS_data_files" > Click "MR Session" next to "OASIS3_data_files" > Select all entries ("demo through DUT") > Bulk Actions: Download.
  - Unzip the downloaded folder, and place it in the `rawdata` folder.  This folder should be called "OASIS3_data_files" (it should not have to be renamed).

