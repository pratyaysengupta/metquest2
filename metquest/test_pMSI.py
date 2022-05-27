from MSI import calculate_pairwiseMSI

# path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/sbml-v2/test'
# path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/sbml/narrowset/HC'
# path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/guaymas_microbiome/Guaymas_metmodels_v2/allmedia_guaymas'
# path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/carved/ocular_narrow_set/BKCR'
path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/guaymas_microbiome/Guaymas_metmodels_v2/test'

# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/media/seed_metabolites_test.txt'
# sd_file = 'HC_minmedia_ce.txt'
sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/guaymas_microbiome/new_JW1_ce.txt'
# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/sbml-v2/narrowset/HC/Run2/HC_minmedia_ce.txt'

calculate_pairwiseMSI(path, sd_file)
