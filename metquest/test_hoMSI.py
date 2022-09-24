from metquest import calculate_higherorderMSI

path = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS'
# path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/carved/ocular_narrow_set/BKSW'
# path = r'C:\Users\dines\IIT-Madras(IC&SR)\IITM-LVPEI - IITM - IITM\occular_microbiome\met_models\sbml-v2\test'
# path = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/guaymas_microbiome/\
# Guaymas_metmodels_v2/allmedia_guaymas/'

# sd_file='/data/shreyansh/MSI/JW1_media.txt'
# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/media/\
# seed_metabolites_test.txt'
# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/media/seed_mets_test.txt'
# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/guaymas_microbiome/JW1_media.txt'
# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/sbml-v2/narrowset/combined/Run2/BK_minmedia_ce.txt'
# sd_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/occular_microbiome/met_models/carved/ocular_narrow_set/combined_BK/carvedCOMBK_minmambo.txt'
# sd_file = r'C:\Users\dines\IIT-Madras(IC&SR)\IITM-LVPEI - IITM - IITM\occular_microbiome\met_models\sbml-v2\narrowset\HC\Run2\HC_minmedia_ce.txt'
sd_file = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS\MOPS_Glu.txt'
cluster_file = 'individual_clusters'
# # cluster_file = 'C:/Users/DINESH KUMAR/Desktop/STUDY/MS/gut/guaymas_microbiome/guaymas__31clusters.txt'
# cluster_file = r'C:\Users\dines\IIT-Madras(IC&SR)\IITM-LVPEI - IITM - IITM\occular_microbiome\met_models\sbml-v2\test\clusters.csv'
# cluster_file = '15clusters.txt'
# # cluster_file='/data/shreyansh/MSI/guaymas__31clusters.txt'

calculate_higherorderMSI(path, sd_file, cluster_file)
