from metquest import calculate_pairwiseMSI


# path = r'C:\Users\dines\IIT-Madras(IC&SR)\IITM-LVPEI - IITM - IITM\occular_microbiome\met_models\sbml-v2\test'
# path = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS'

# sd_file = r'C:\Users\dines\IIT-Madras(IC&SR)\IITM-LVPEI - IITM - IITM\occular_microbiome\met_models\sbml-v2\narrowset\HC\Run2\HC_minmedia_ce.txt'
# sd_file = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS\MOPS_Glu.txt'


path = r'C:\Users\dines\PycharmProjects\share_code\msi\example\models_folder'
sd_file = r'C:\Users\dines\PycharmProjects\share_code\msi\example\seed_metabolites_file.txt'
calculate_pairwiseMSI(path, sd_file)
