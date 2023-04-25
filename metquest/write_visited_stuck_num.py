# Writes the number of visited rxns and number of stuck rxns of all microbes in a file
import metquest as mq
import cobra
import glob
import os
import pandas as pd

path = r'C:\Users\dines\Documents\plant-microbiome\models\test'
seedmet_file = r'C:\Users\dines\Documents\plant-microbiome\models\MOPS\MOPS_Glu.txt'
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

org_info, _, full_name_map, vis = mq.pairwiseMSI.find_stuck_rxns(model, file_names, seedmet_file, 1)

organism_stats = {}
for org in org_info:
    organism_stats[org] = []
    organism_stats[org].append(len(list(set(org_info[org]))))
    organism_stats[org].append(len(list(set(vis[org]))))

df1 = pd.DataFrame.from_dict(organism_stats, orient="index")
df1.columns = ['Stuck rxns', 'Visited rxns']
df1.to_csv('./results/organism_stats_'+os.path.basename(seedmet_file).replace('.txt', '.csv'))
