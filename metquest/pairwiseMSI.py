import os
import sys
import glob
import warnings
import itertools
import metquest as mq
from metquest import find_transport_rxns


def find_relieved_rxn(model, seed_name, org_info_single, org_info_pair):
    """
    This function extracts and writes the relieved rxns into a tsv file
    :param model:
    :param seed_name: name of the media used (identifer to know what media is used when analysis is done using multiple media)
    :param org_info_single: Dictionary containing stuck reactions of all microbes in the community
    :param org_info_pair: Dictionary containing stuck reactions of all microbes in the community
    :return: None
    """
    relieved = {}
    for org1 in model:
        for org2 in model:
            if org1.id + '_' + org2.id in org_info_pair.keys():
                relieved[org1.id + '_' + org2.id] = []
                temp = list(set(org_info_single[org1.id + '_' + org1.id]) - set(org_info_pair[org1.id + '_' + org2.id]))
                for j in temp:
                    relieved[org1.id + '_' + org2.id].append(j)
            else:
                continue

    rel_rxns_name = {}
    detailed_rel_rxns = {}
    for i in model:
        rxn_ids = []
        for r in i.reactions:
            rxn_ids.append(r.id)
        for j in model:
            org1 = i.id
            org2 = j.id
            if org1 + '_' + org2 in relieved.keys():
                detailed_rel_rxns[org1 + '_' + org2] = []
                rel_rxns_name[org1 + '_' + org2] = []
                for rel in relieved[org1 + '_' + org2]:
                    rel_rxn = i.reactions[rxn_ids.index(rel)].reaction
                    detailed_rel_rxns[org1 + '_' + org2].append(rel_rxn)
                    rel_rxns_name[org1 + '_' + org2].append(i.reactions[rxn_ids.index(rel)].name)

    relieved_rxn_output_file = 'results/relieved_rxns_' + seed_name + '_w_excrxns.tsv'
    with open(relieved_rxn_output_file, 'w') as g:
        header = 'acceptor\tdonor\trelieved reactions\n'
        g.write(header)
        for i in model:
            for j in model:
                org1 = i.id
                org2 = j.id
                if org1 + '_' + org2 in relieved.keys():
                    g.write(org1 + '\t' + org2 + '\t')
                    rel_rxns = list(set(relieved[org1 + '_' + org2]))
                    det_rel_rxns = list(set(detailed_rel_rxns[org1 + '_' + org2]))
                    rel_rxn_nam = list(set(rel_rxns_name[org1 + '_' + org2]))
                    for x in rel_rxns:
                        g.write(x + '\t')
                    g.write('\n')
                    g.write('\t' + '\t')
                    for d in rel_rxn_nam:
                        g.write(d + '\t')
                    g.write('\n')
                    g.write('\t' + '\t')
                    for k in det_rel_rxns:
                        g.write(k + '\t')
                    g.write('\n')
    print('relieved reactions are written at:\n', relieved_rxn_output_file)


def find_stuck_rxns(model, community, seedmet_file, no_of_orgs):
    """
    Constructs graphs using MetQuest and finds all stuck reactions in the cellular compartment
    :param model:
    :param community: list of GSMM files
    :param seedmet_file: path to txt file containing seed metabolites
    :param no_of_orgs: number of organisms in a community
    :return:
        org_info: Dictionary containing stuck reactions of all microbes in the community
        scope: Dictionary containing all the metabolites that can be produced by the microbes in the community
        namemap: Dictionaru containing all the decrypted rxn ids
    """

    warnings.filterwarnings("ignore")
    G, full_name_map = mq.construct_graph.create_graph(community, no_of_orgs)
    if not os.path.exists('results'):
        os.makedirs('results')
    f = open(seedmet_file, 'r')

    # Reading seed metabolites

    seedmets, temp_seedmets = [], []
    while True:
        l = f.readline().strip()
        if l == '': break
        temp_seedmets.append(l)
    f.close()
    for m in model:
        for i in temp_seedmets:
            seedmets.append(m.id + ' ' + i)
    seedmets = set(seedmets)

    all_possible_combis = list(itertools.combinations(list(range(len(community))), int(no_of_orgs)))
    if no_of_orgs > 1 and sorted(community)[0][0] == '0':
        all_possible_combis = all_possible_combis[:len(community) - 1]
    org_info = {}
    scope = {}
    vis = {}
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
            model2 = model[all_possible_combis[i][j - 1]].id
            org_info[model1 + '_' + model2] = stuck
            scope[model1 + '_' + model2] = s
            vis[model1 + '_' + model2] = visited
    return org_info, scope, full_name_map, vis


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


def pMSI(community, sd_file):
    """
    Calculates MSI for CarveMe models
    Extracts and writes relieved reactions in every pair
    :param community: list of GSMM files
    :param sd_file: path to txt file containing seed metabolites
    :return: msi: Dictionary containing MSI values for every pair
    """
    # find all transport reactions
    transport_rxns, model = find_transport_rxns(community)
    # find stuck reactions
    org_info_single, scope_sin, namemap_sin, vis = find_stuck_rxns(model, community, sd_file, 1)
    community = glob.glob('*.xml')
    if not community:
        community = glob.glob('*.sbml')
    community.sort()
    org_info_pair, scope_pair, namemap_pair, vis = find_stuck_rxns(model, community, sd_file, 2)
    # decrypt the stuck reactions
    org_info_single = decrypt_org_info(org_info_single, namemap_sin)
    org_info_pair = decrypt_org_info(org_info_pair, namemap_pair)
    # Filter out the transport reactions from every stuck reaction list
    org_info_single_wo_trans_rxn, org_info_pair_wo_trans_rxn = {}, {}
    for i in org_info_single:
        org_info_single_wo_trans_rxn[i] = list(set(org_info_single[i]) - set(transport_rxns))
    for i in org_info_pair:
        org_info_pair_wo_trans_rxn[i] = list(set(org_info_pair[i]) - set(transport_rxns))
    # find all the relieved reactions in every pairs
    find_relieved_rxn(model, os.path.basename(sd_file).replace('.txt', ''), org_info_single,
                      org_info_pair)
    # calculate MSI for every pair
    msi = {}
    for org1 in model:
        stuck_A = len(org_info_single_wo_trans_rxn[org1.id + '_' + org1.id])
        for org2 in model:
            if org1.id + '_' + org2.id in org_info_pair_wo_trans_rxn.keys():
                stuck_AUB = len(org_info_pair_wo_trans_rxn[org1.id + '_' + org2.id])
                if stuck_A == 0:
                    msi[org1.id + '_' + org2.id] = 0
                else:
                    msi[org1.id + '_' + org2.id] = 1 - (stuck_AUB / stuck_A)
    return msi, model


def calculate_pairwiseMSI(path, sd_file):
    """
    This function calculates pairwise-MSI for all given microbes.

    Creates a csv file containing the MSI values of all pairs.

    Creates an tsv file containing the list of reaction relieved
    in all acceptor microbes in the presence of corresponding donor microbes.

    :param path: path to all xml files
    :param sd_file: path to txt file containing seed metabolites
    """

    warnings.filterwarnings("ignore")
    os.chdir(path)
    file_names = glob.glob('*.xml')
    if not file_names:
        file_names = glob.glob('*.sbml')

    if not file_names:
        print("There are no sbml files. Please check the path")
        return
    print("Filenames", file_names)
    sys.path.append(path)
    community = file_names.copy()
    community.sort()

    msi, model = pMSI(community, sd_file)

    msi_output_file = 'results/MSI_' + os.path.basename(sd_file).replace('.txt', '') + '.csv'
    with open(msi_output_file, 'w') as f:
        header = 'organism,in_the_presence,msi_value\n'
        f.write(header)
        for org1 in model:
            for org2 in model:
                if org1.id + '_' + org2.id in msi.keys():
                    f.write(org1.id + ',' + org2.id + ',' + str(msi[org1.id + '_' + org2.id]) + '\n')
    print('MSI values are written at:\n', msi_output_file)
