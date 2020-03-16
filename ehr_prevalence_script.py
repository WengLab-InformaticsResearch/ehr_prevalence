from ehr_prevalence import *
from datetime import datetime
import gc


# Settings
# Directory where data are stored
data_dir = r'/home/user_xxxx/cohd/data'
# Directory where results will be stored
results_dir = r'/home/user_xxxx/cohd/results'
# File names of database dump files
file_person = 'person.txt'
file_concept_patient = 'unique_patient_concept_pairs_date.txt'
file_concept = 'concepts.txt'
file_descendants_dir = 'concept_descendants_direct.txt'
file_descendants_all = 'concept_descendants_all_observed.txt'
# Where were the database dump files generated from? 
# Affects how the dump files are read. Options:
#     ssms: Microsoft SQL Server Management Studio
#     mysql: MySQL
database = 'ssms'  
# Range of years to include to calculate the 5-year dataset (inclusive)
range_5year = (2013, 2017)
# Range of years to include to calculate the  lifetime dataset (inclusive)
range_lifetime = (1985, 2018)
# Randomize
randomize = True
# Minimum count for a concept to be included (inclusive)
min_count = 11

# Timestamp to label output files with
timestamp = '_' + datetime.now().strftime("%Y%m%d-%H%M")

# Set up logging
logging_setup(results_dir)

# Load the data
file_concept = os.path.join(data_dir, file_concept)
concepts = load_concepts(file_concept, database)
file_descendants_dir = os.path.join(data_dir, file_descendants_dir)
descendants_dir = load_descendants(file_descendants_dir, database)
file_descendants_all = os.path.join(data_dir, file_descendants_all)
descendants_all = load_descendants(file_descendants_all, database)
file_person = os.path.join(data_dir, file_person)
patient_info = load_patient_data(file_person, database)
file_concept_patient = os.path.join(data_dir, file_concept_patient)
cp_data = load_concept_patient_data(file_concept_patient, database, patient_info)

# Attempt to free up some memory
gc.collect()

# # Basic quality analysis
# quality_analysis(results_dir, cp_data, concepts, min_count=0)


# 5-year dataset
cp_data_5year = merge_concepts_years(cp_data, range_5year[0], range_5year[1])
cp_data_5year_hier = merge_ranged_concept_descendants(cp_data_5year, concepts, descendants_dir)
concepts_5year = single_concept_ranged_counts(results_dir, cp_data_5year_hier, randomize, min_count,
                                              additional_file_label='hierarchical' + timestamp)
single_concept_yearly_deviation(results_dir, cp_data, concepts_5year, range_5year, randomize=True,
                                file_label='hierarchical' + timestamp)
concept_pairs_5y = paired_concept_ranged_counts(results_dir, cp_data_5year_hier, randomize, min_count,
                                                additional_file_label='hierarchical' + timestamp)
paired_concept_yearly_deviation(results_dir, cp_data, concept_pairs_5y, range_5year, randomize=True,
                                file_label='hierarchical' + timestamp)

# Attempt to free up some large chunks of memory
del cp_data_5year
del cp_data_5year_hier
gc.collect()


# Lifetime dataset
cp_data_life = merge_concepts_years(cp_data, range_lifetime[0], range_lifetime[1])
cp_data_life_hier = merge_ranged_concept_descendants(cp_data_life, concepts, descendants_dir)
concepts_life = single_concept_ranged_counts(results_dir, cp_data_life_hier, randomize, min_count,
                                             additional_file_label='hierarchical' + timestamp)
range_lifetime_deviation = (1986, 2017)  # Only use data from full years for variance
single_concept_yearly_deviation(results_dir, cp_data, concepts_life, range_lifetime_deviation, randomize=True,
                                file_label='hierarchical' + timestamp)
concept_pairs_life = paired_concept_ranged_counts(results_dir, cp_data_life_hier, randomize, min_count,
                                                  additional_file_label='hierarchical' + timestamp)
paired_concept_yearly_deviation(results_dir, cp_data, concept_pairs_life, range_lifetime_deviation, randomize=True,
                                file_label='hierarchical' + timestamp)

# Attempt to free up some large chunks of memory
del cp_data_life
del cp_data_life_hier
gc.collect()


# Annual counts
# single_concept_yearly_counts(results_dir, cp_data, randomize, min_count)
# paired_concept_yearly_counts(results_dir, cp_data, randomize, min_count)

