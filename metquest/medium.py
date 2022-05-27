import glob
import os
import re
import cobra
import warnings
import itertools
import metquest as mq
from cobra.medium import minimal_medium

def minimal_media_from_cobrapy(path, outputfilename, essential_mets):
    """
    This function creates the seed metabolite file
    :param path:
    :param outputfilename:
    :param essential_mets:
    :return: A txt file containing all the seed metabolites
    """
    os.chdir(path)
    file_names = glob.glob('*.xml')
    if not file_names:
        file_names = glob.glob('*.sbml')

    if not file_names:
        print("There are no .xml files. Please check the path")
    print("Filenames", file_names)

    media, model = [], []
    for i in range(len(file_names)):
        model.append(cobra.io.read_sbml_model(file_names[i]))
        media += list(minimal_medium(model[i], 0.1, minimize_components=True).keys())
    media = list(set(media))
    # exc_rxns = model[0].exchanges
    # met_id = list(exc_rxns[0].metabolites.keys())[0].id
    # pattern = re.compile(r'[_[][a-z]\d*[]]*')
    # exc_marker = pattern.findall(met_id)
    with open(outputfilename, 'w') as f:
        for i in media:
            w = i[3:].replace('_e', '_c')
            f.write(w + '\n')
            w = i[3:].replace('_c', '_e')
            f.write(w + '\n')

    # org_info_single, scope_sin, namemap_sin, _ = find_stuck_rxns(model, file_names, outputfilename, 1)
    # missing_seed = add_missing_seed(model, essential_mets, scope_sin)
    # # media = media + missing_seed
    # for i in range(len(missing_seed)):
    #     missing_seed[i] = missing_seed[i].split(' ')[1]
    # missing_seed = list(set(missing_seed))
    # with open(outputfilename, 'a') as f:
    #     for m in missing_seed:
    #         f.write(m + '\n')


def merge_orgname(i, j):
    return i+' '+j

def add_missing_seed(model, essential_mets, scope):
    """
    Finds and returns essential metabolites that are missing in the seed metabolites
    :param model: List of model identifiers
    :param essential_mets: List of essential metabolites
    :param scope: Dictionary containing list of metabolites that can be produced by the microbe
    :return: Returns the list of missing essential metabolites
    """
    missing_seed = []
    with open(essential_mets, 'r') as f:
        mets=[]
        while True:
            l = f.readline().strip()
            if l == '': break
            mets.append(l)
    for orgname in model:
        name = orgname.id
        if name+'_'+name in scope.keys():
            new_mets = list(map(merge_orgname, [name]*len(mets), mets))
            missing_seed = missing_seed + list(set(new_mets) - set(scope[name+'_'+name]))
    return list(set(missing_seed))


def find_mets_not_produced(model, scope):
    """
    This function extracts the metabolites not produced by the models (mets not in scope)
    :param model: list of models
    :param scope: Dictionary of metabolites produced by every microbe
    :return:
    """
    metabolites = {}
    mets_not_produced = {}
    for i in model:
        metabolites[i.id+'_'+i.id] = []
        for met in i.metabolites:
            metabolites[i.id+'_'+i.id].append(i.id+' '+met.id)
    for org in scope:
        mets_not_produced[org] = list(set(metabolites[org]) - set(scope[org]))
        for i in mets_not_produced[org]:
            if i.find('_e') >= 0:
                mets_not_produced[org].remove(i)
        mets_not_produced[org].sort()
    return mets_not_produced

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