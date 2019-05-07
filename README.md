# Preprocess case and control cohorts before running CAF

Note this is a first commit of what I'm currently working on right now at CGR in Astrazeneca ...

## Objective
The idea of this pipeline is to preprocess data before running the collapsing analysis (CAF, https://confluence.astrazeneca.net/display/CFGR/Collapsing+Analysis). 

Given a list of cases and a list of controls, the **aim** here is to harmonise as much as possible and filter out samples with consent withdrawn, contamination or mismatch between self reported and CGR estimated gender, for instance. It also harmonise samples so as to match gender or age distributions in both lists. Basically, the list of controls is sampled depending on the distribution of what the user asks based on the distribution along the cases.

## Input

- Configuration file -- features/parameters to analyse together with their thresholds, paths to the files including the list of cases and controls, path to the QC metadata


## Output

- New file with the new list of cases for running CAF
- New file with the new list of controls for running CAF
- log output file including the features and thresholds criteria used in the preprocessing and the CGRSequenceIDs affected with different filters

## Filters

This is the full list of the filters we do (some of them depend on the thresholds defined in the configuration file).

- Contamination: we filter out samples having a VerifyBAMId-freemix value >= threshold (by default 0.04)
- Reported vs genetic gender: we filter out samples with gender mismatch comparing reported vs estimated


## Harmonisation

This is the criteria we can ask when harmonising data in both case and control cohorts (i.e. having similar distributions in both subsets):

- By age
- By gender
- By ancestry
- By smoking status
- By BMI

**NOTE** this implies including more columns (i.e. age, gender, ancestry, ...) in the input case/control cohort file.

## How to run
There is a config file where the user can specify the paths to the files containing the lists of cases and controls. There are 2 modes of running the preprocessing: genomic and clinic

### Genomic filtering
This is based on the QC metadata information. Together with the features and thresholds that the user defines in the configuration file some samples are filtered out.

### Clinical or phenotype filtering
This is based on the control and case cohort distributions. The user needs to provide extra columns to the list of cases and controls, depending on the extra filters she/he wants to do: _age_, _gender_, _smoking condition_


### 1. Work on the config file and use it to point to all the necessary files
Use the example [config](https://github.com/AstraZeneca-CGR/preprocess_CAF/blob/master/configuration_files/preprocessing_config.cfg) file as a template to define the locations of all the files and variables.

### 2. Load the virtual environment
```
source venv/bin/activate

```

### 3. Run preprocessing pipeline
```
python preprocess_CAF/running_preprocessing_CAF.py --cfg configuration_files/preprocessing_config.cfg

```

