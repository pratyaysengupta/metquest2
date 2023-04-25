import os
import glob
import sys
import warnings
import itertools
import re
import pandas as pd
import cobra

import metquest as mq
from metquest import find_transport_rxns

warnings.filterwarnings("ignore")


def find_relievedrxns(model, org_info, org_info_pert):
    relieved = {}
    detailed_rel_rxns = {}
    rel_rxns_name = {}
    for i in org_info_pert:
        relieved[i] = list(set(org_info_pert[i]) - set(org_info[i]))

    for i in model:
        j = i.id
        detailed_rel_rxns[j] = []
        rel_rxns_name[j] = []
        if len(relieved[j]):
            rxn_ids = []
            for r in i.reactions:
                rxn_ids.append(r.id)
            for rel in relieved[j]:
                rel_rxn = i.reactions[rxn_ids.index(rel)].reaction
                detailed_rel_rxns[j].append(rel_rxn)
                rel_rxns_name[j].append(i.reactions[rxn_ids.index(rel)].name)

    return relieved, detailed_rel_rxns, rel_rxns_name


def preprocess_seedmets(seedmet_file, model):
    f = open(seedmet_file, 'r')
    seedmets, temp_seedmets = [], []
    while True:
        l = f.readline().strip()
        if l == '': break
        temp_seedmets.append(l)
    f.close()

    exc_rxns = model[0].exchanges
    met_id = list(exc_rxns[0].metabolites.keys())[0].id
    pattern = re.compile(r'[_[][a-z]\d*[]]*')
    exc_marker = pattern.findall(met_id)

    for m in model:
        for i in temp_seedmets:
            seedmets.append(m.id + ' ' + i + exc_marker[0])
            seedmets.append(m.id + ' ' + i + exc_marker[0].replace('e', 'c'))

    return set(seedmets)


def find_stuck_rxns(model, community, seedmet_file, no_of_orgs):
    # Constructing graphs

    warnings.filterwarnings("ignore")
    G, full_name_map = mq.construct_graph.create_graph(community, no_of_orgs)
    if not os.path.exists('results'):
        os.makedirs('results')
    f = open(seedmet_file, 'r')

    seedmets = preprocess_seedmets(seedmet_file, model)

    all_possible_combis = list(itertools.combinations(list(range(len(community))), int(no_of_orgs)))
    if no_of_orgs > 1 and sorted(community)[0][0] == '0':
        all_possible_combis = all_possible_combis[:len(community) - 1]
    org_info = {}
    scope = {}
    print('No. of graphs constructed: ', len(G))

    # This loop finds all the stuck reaction

    for i in range(len(all_possible_combis)):
        lbm, sd, s = mq.guided_bfs.forward_pass(G[i], seedmets)
        for j in range(len(all_possible_combis[i])):
            stuck = []
            rxnNode = []
            model1 = model[all_possible_combis[i][j]].id
            visited = list(sd.keys())
            for r in G[i].nodes:
                if r.find(model1) >= 0:
                    rxnNode.append(r)
            for rxn in rxnNode:
                if rxn in visited:
                    continue
                elif rxn.find('ERR') >= 0:
                    continue
                elif rxn.find('Org') >= 0:
                    if (rxn[len(model1) + 5] == 'I') or (rxn[len(model1) + 5] == 'R'):
                        stuck.append(rxn)
            org_info[model1] = stuck
            scope[model1] = s
    return org_info, scope, full_name_map


def decrypt_org_info(org_info, namemap):
    """
    This function decrypts the rxn ids using the data in corresponding namemaps
    :param org_info:
    :param namemap:
    :return:
        org_info: An dictionary of decrypted rxn ids for each community
    """
    for i in org_info:
        for j in range(len(org_info[i])):
            org_info[i][j] = namemap[org_info[i][j]]
    return org_info


def make_perturbed_community(rem_org, pert_models, pert_community):
    pert_model_ids = [i.id for i in pert_models]
    for i in rem_org:
        if i in pert_model_ids:
            pert_models.remove(pert_models[pert_model_ids.index(i)])
            pert_community.remove(pert_community[pert_model_ids.index(i)])
            pert_model_ids.remove(i)

    return pert_models, pert_community, pert_model_ids


def perform_task(cluster_file, sd_file, model, transport_rxns, pert_community,
                 org_info_wo_trans_rxn, rem_org_list, n):
    org_info_pert, scope_pert, namemap_pert = \
        find_stuck_rxns(model, pert_community, sd_file, len(pert_community))
    org_info_pert = decrypt_org_info(org_info_pert, namemap_pert)
    org_info_pert_wo_trans_rxn = {}
    for i in org_info_pert:
        org_info_pert_wo_trans_rxn[i] = list(set(org_info_pert[i]) - set(transport_rxns))

    g = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') + '_' +
             os.path.basename(sd_file).replace('.txt', '') + '/Community_without_clus' + str(n) + '.csv', 'w')
    for m in org_info_pert_wo_trans_rxn:
        g.write(m + ',' + str(len(org_info_pert_wo_trans_rxn[m])) + '\n')
    g.close()
    stuck_com = 0
    stuck_pert_com = 0
    for i in org_info_wo_trans_rxn:
        if i not in rem_org_list:
            stuck_com += len(org_info_wo_trans_rxn[i])
    for i in org_info_pert_wo_trans_rxn:
        stuck_pert_com += len(org_info_pert_wo_trans_rxn[i])
    msi = 1 - (stuck_com / stuck_pert_com)
    print(n, 'th cluster')
    return org_info_pert, org_info_pert_wo_trans_rxn, msi


def write_relieved_rxns(g, relieved, detailed_rel_rxns, rel_rxns_name):
    g.write('acceptor\trelieved reactions\n')

    for i in relieved:
        g.write(i + '\t')
        rel_rxns = list(set(relieved[i]))
        det_rel_rxns = list(set(detailed_rel_rxns[i]))
        rel_rxn_nam = list(set(rel_rxns_name[i]))
        for j in rel_rxns:
            g.write(j + '\t')
        g.write('\n')
        g.write('\t')
        for d in rel_rxn_nam:
            g.write(d + '\t')
        g.write('\n')
        g.write('\t')
        for k in det_rel_rxns:
            g.write(k + '\t')
        g.write('\n')


def write_relieved_rxn_metadata(h, org_info_wo_trans_rxn, org_info_pert_wo_trans_rxn):
    nrelieved = {}
    for i in org_info_pert_wo_trans_rxn:
        nrelieved[i] = len(org_info_pert_wo_trans_rxn[i]) - len(org_info_wo_trans_rxn[i])
        if nrelieved[i]:
            h.write(i + ',' + str(len(org_info_wo_trans_rxn[i])) + ',' + str(
                len(org_info_pert_wo_trans_rxn[i])) + ',' + str(nrelieved[i]) + '\n')


def calculate_higherorderMSI(path, sd_file, cluster_file):
    os.chdir(path)
    file_names = glob.glob('*.xml')
    if not file_names:
        file_names = glob.glob('*.sbml')

    if not file_names:
        print("There are no .xml files. Please check the path")
    print("Filenames", file_names)
    sys.path.append(path)
    community = file_names
    community.sort()
    transport_rxns, model = find_transport_rxns(community)
    org_info, scope, namemap = find_stuck_rxns(model, community, sd_file, len(community))
    org_info = decrypt_org_info(org_info, namemap)
    org_info_wo_trans_rxn = {}
    for i in org_info:
        org_info_wo_trans_rxn[i] = list(set(org_info[i]) - set(transport_rxns))

    if not os.path.exists('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') +
                          '_' + os.path.basename(sd_file).replace('.txt', '') + '/data_analysis'):
        os.makedirs('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') + '_' +
                    os.path.basename(sd_file).replace('.txt', '') + '/data_analysis')
    f = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') + '_' +
             os.path.basename(sd_file).replace('.txt', '') + '/community_unperturbed.csv', 'w')
    for i in org_info_wo_trans_rxn:
        f.write(i + ',' + str(len(org_info_wo_trans_rxn[i])) + '\n')
    f.close()

    if cluster_file == 'individual_clusters':
        # file_names = glob.glob('*.xml')
        # if not file_names:
        #     file_names = glob.glob('*.sbml')
        rem_org_list1 = {}
        rem_org_list2 = {}

        for i in range(len(model)):
            rem_org_list1[i] = [model[i].id]
            rem_org_list2[i] = [model[i].id]

    else:
        cluster_data = pd.read_csv(cluster_file, sep=',')
        rem_org_list1 = cluster_data.set_index('Cluster').T.to_dict('list')
        for n in rem_org_list1:
            rem_org_list1[n] = [j for j in rem_org_list1[n] if pd.isna(j) is False]
        # model_ids = [i.id for i in model]
        for n in rem_org_list1:
            rem_org_list1[n] = [cobra.io.read_sbml_model(i).id for i in rem_org_list1[n]]
            # rem_org_list1[n] = [model_ids[model_ids.index(i)] for i in rem_org_list1[n]]

        rem_org_list2 = rem_org_list1.copy()

    for nclus in rem_org_list2:
        rem_org_list2[nclus] = [x.replace('.xml', '') for x in rem_org_list2[nclus]]

    f = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') + '_' +
             os.path.basename(sd_file).replace('.txt', '') + '/higher_order_msi.csv', 'w')
    for n in rem_org_list1:
        os.chdir(path)
        new_models = model.copy()
        new_community = glob.glob('*.xml')
        if not new_community:
            new_community = glob.glob('*.sbml')
        new_community.sort()

        pert_models, pert_community, pert_model_ids = make_perturbed_community(rem_org_list1[n], new_models, new_community)

        org_info_pert, org_info_pert_wo_trans_rxn, msi = perform_task(cluster_file, sd_file, pert_models, transport_rxns,
                                                                      pert_community, org_info_wo_trans_rxn,
                                                                      rem_org_list2[n], n)
        for i in rem_org_list2[n]:
            f.write('Comm,clus_' + str(n) + '#' + i + ',' + str(msi) + '\n')

        if msi:
            relieved, detailed_rel_rxns, rel_rxns_name = find_relievedrxns(pert_models, org_info, org_info_pert)
            g = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '')
                     + '_' + os.path.basename(sd_file).replace('.txt', '') +
                     '/data_analysis/relieved_rxns_Comm--clus' + str(n) + '.tsv', 'w')
            write_relieved_rxns(g, relieved, detailed_rel_rxns, rel_rxns_name)
            g.close()

            h = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') + '_' +
                     os.path.basename(sd_file).replace('.txt', '') + '/data_analysis/Comm--clus' + str(n) + '.csv', 'w')
            h.write('Comm--clus' + str(n) + '\n')
            for i in rem_org_list2[n]:
                h.write(i + '\n')
            h.write('num of rxns relieved in the below orgs in the presence of clust' + str(n) + '\n')
            h.write('org,unpert,clust_' + str(n) + 'KO,rxns relieved\n')
            write_relieved_rxn_metadata(h, org_info_wo_trans_rxn, org_info_pert_wo_trans_rxn)
            h.close()
            print('Comm--clus' + str(n))

        os.chdir(path)
        new_models = model.copy()
        new_community = glob.glob('*.xml')
        if not new_community:
            new_community = glob.glob('*.sbml')
        new_community.sort()
        ko_models, ko_community, model_ids = make_perturbed_community(pert_model_ids, new_models, new_community)
        ko_org_list = [x for x in pert_model_ids]
        if len(ko_org_list) < len(model):
            org_info_pert, org_info_pert_wo_trans_rxn, msi = perform_task(cluster_file, sd_file, ko_models,
                                                                          transport_rxns, ko_community,
                                                                          org_info_wo_trans_rxn, ko_org_list, n)
            for i in ko_community:
                f.write('clus_' + str(n) + '#' + i + ',Comm,' + str(msi) + '\n')

            if msi:
                relieved, detailed_rel_rxns, rel_rxns_name = find_relievedrxns(ko_models, org_info, org_info_pert)
                g = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '')
                         + '_' + os.path.basename(sd_file).replace('.txt', '') +
                         '/data_analysis/relieved_rxns_clus--Comm' + str(n) + '.tsv', 'w')
                write_relieved_rxns(g, relieved, detailed_rel_rxns, rel_rxns_name)
                g.close()

                h = open('results/clusterKO_' + os.path.basename(cluster_file).replace('.csv', '') + '_' +
                         os.path.basename(sd_file).replace('.txt', '') + '/data_analysis/clus' + str(n) + '--Comm.csv',
                         'w')
                h.write('clus' + str(n) + '--Comm\n')
                for i in ko_org_list:
                    h.write(i + '\n')
                h.write('num of rxns relieved in the below orgs in the presence of Comm')
                h.write('org,unpert,commKO,rxns relieved\n')
                write_relieved_rxn_metadata(h, org_info_wo_trans_rxn, org_info_pert_wo_trans_rxn)
                h.close()
                print('clus' + str(n) + '--Comm')

    f.close()
