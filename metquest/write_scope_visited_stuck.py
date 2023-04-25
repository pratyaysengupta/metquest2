# Writes all the scope, visited rxns and stuck rxns of all microbes in separate files
import metquest as mq
import cobra
import glob
import os
import pandas as pd

path = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS'
seedmet_file = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS\plant_minmedia.txt'
os.chdir(path)

file_names = glob.glob('*.xml')
if not file_names:
    file_names = glob.glob('*.sbml')

if not file_names:
    print("There are no sbml files. Please check the path")

print("Filenames", file_names)

model = []
for i in range(len(file_names)):
    model.append(cobra.io.read_sbml_model(file_names[i]))

org_info, scope, full_name_map, vis = mq.pairwiseMSI.find_stuck_rxns(model, file_names, seedmet_file, 1)

org_info_decrypt = mq.pairwiseMSI.decrypt_org_info(org_info, full_name_map)
vis_decrypt = mq.pairwiseMSI.decrypt_org_info(vis, full_name_map)

df1 = pd.DataFrame.from_dict(org_info_decrypt, orient="index")
df1.T.to_csv('./results/'+os.path.basename(seedmet_file).replace('.txt', '_stuckrxns.csv'))
df2 = pd.DataFrame.from_dict(scope, orient="index")
df2.T.to_csv('./results/'+os.path.basename(seedmet_file).replace('.txt', '_scope.csv'))
df3 = pd.DataFrame.from_dict(vis_decrypt, orient="index")
df3.T.to_csv('./results/'+os.path.basename(seedmet_file).replace('.txt', '_visitedrxns.csv'))
