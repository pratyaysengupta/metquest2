import pandas as pd
import re
import glob
import seaborn as sns
import matplotlib.pyplot as plt
import os
from metquest import construct_graph
from metquest import find_stuck_rxns
from metquest import get_models
from pyvis.network import Network
from pyvis.options import Layout


def get_exc_metabolites(filename, seedmet_file):
    with open(filename, 'r') as f:
        temp = f.read().strip()
    data = temp.split('\n')

    with open(seedmet_file, 'r') as g:
        temp = g.read().strip()
    seedmet = temp.split('\n')

    exc_mets = {}
    for i in range(1, len(data), 3):
        line = data[i].strip().split('\t')
        microbes = line[0:2]
        # rxn_collection = line[2:]
        try:
            line = data[i + 2].strip().split('\t')
            full_rxns = line.copy()

            for j in range(len(microbes)):
                microbes[j] = microbes[j] + '.xml'

            metabolites = []
            prog = re.compile('\\b[\\S|\\d]+_e\\d*\\b')
            for j in full_rxns:
                result = prog.findall(j)
                [metabolites.append(a) for a in result]

            exc_mets[microbes[0] + ',' + microbes[1]] = list(set(metabolites))
            exc_mets[microbes[0] + ',' + microbes[1]] = list(set(metabolites) - set(seedmet))
        except IndexError:
            continue
    return exc_mets


def write_exc_metabolites(path_to_models, relrxns_filename, seedmet_file):
    exc_mets = get_exc_metabolites(relrxns_filename, seedmet_file)
    community = glob.glob(path_to_models+'/*.xml')
    extension = '.xml'
    if not community:
        community = glob.glob(path_to_models+'/*.sbml')
        extension = '.sbml'
    community.sort()
    model = get_models(community)

    org_info_single, scope_sin, namemap_sin, vis = find_stuck_rxns(model, community, seedmet_file, 1)
    refined_exc_mets = {}
    for i in exc_mets:
        acceptor = i.split(',')[0].replace(extension, '')
        scope_sin[acceptor + '_' + acceptor] = [x.replace(acceptor+' ', '').replace('_c0', '_e0')
                                                for x in scope_sin[acceptor + '_' + acceptor]]
        refined_exc_mets[i] = [x for x in exc_mets[i] if x not in scope_sin[acceptor + '_' + acceptor]]

    df_exc_mets = pd.DataFrame.from_dict(refined_exc_mets, orient="index")
    df_exc_mets.to_csv(relrxns_filename.replace('.tsv', '') + '_refined_exc_mets.csv')


def get_excmet_stats(filename):
    with open(filename, 'r') as f:
        data = f.read().strip()
    table = data.split('\n')
    line = []
    for i in range(1, len(table)):
        line.append(table[i].split(','))
    mets = {}
    for i in range(len(line)):
        for j in range(2, len(line[i])):
            if line[i][j]:
                mets[line[i][j]] = 0

    for i in line:
        for j in mets:
            if j in i:
                mets[j] += 1
    return mets


def write_excmet_stats(filename):
    mets = get_excmet_stats(filename)
    exc_mets = pd.DataFrame(mets.items(), columns=['metabolites', 'exchange frequency'])
    exc_mets.to_csv(filename.replace('.csv', '_stats.csv'))


def plot_excmet_count(filename):
    mets = get_excmet_stats(filename)
    exc_mets = pd.DataFrame(mets.items(), columns=['metabolites', 'exchange frequency'])
    sns.set(rc={'figure.figsize': (16, 8.35)})
    ax = sns.barplot(exc_mets['metabolites'], exc_mets['exchange frequency'],
                     palette='Set1',
                     order=exc_mets.sort_values('exchange frequency', ascending=False).metabolites,
                     )
    plt.xticks(rotation=45, horizontalalignment='right')
    ax.tick_params(axis='both', which='major', labelsize=9)
    plt.tight_layout()
    title = os.path.basename(filename)
    plt.title(title[:-22].replace('relieved_rxns', 'Metabolites_exchanged'))
    fig = ax.get_figure()
    fig.savefig(filename.replace('.csv', '_stats'))
    fig.clf()


def draw_graph(acceptor, donor, relievedrxn_filename, outfile):
    with open(relievedrxn_filename, 'r') as f:
        temp = f.read().strip()
    data = temp.split('\n')

    exc_mets = {}
    rxn_collection = 0
    microbes = 0
    for i in range(1, len(data), 3):
        line = data[i].strip().split('\t')
        if line[0]+'.xml' == os.path.basename(acceptor).replace('-', '_') and line[1]+'.xml' == os.path.basename(donor).replace('-', '_'):
            microbes = line[0:2]
            rxn_collection = line[2:]
        else:
            continue

    G, namemap_comm = construct_graph.create_graph([acceptor], 1)
    df = pd.DataFrame.from_dict(namemap_comm, orient="index")
    df.to_csv(outfile+'_namemap.csv')
    new_namemap_comm = {v: k for k, v in namemap_comm.items()}

    initial_nodes = []
    for j in rxn_collection:
        initial_nodes.append(new_namemap_comm[j])
    subgraph_nodes = []
    for n in initial_nodes:
        subgraph_nodes.append(n)
        for k in G[0].predecessors(n):
            subgraph_nodes.append(k)
        for j in G[0].successors(n):
            subgraph_nodes.append(j)
    exc_mets_temp = [m for m in subgraph_nodes if '_e' in m if 'Org_' not in m if '_c' not in m if '_p' not in m]
    exc_mets[microbes[0]+','+microbes[1]] = list(set(exc_mets_temp))

    extracted_subgraph = G[0].subgraph(subgraph_nodes)

    net = Network('1000px', '1000px', notebook=True, directed=True)
    net.from_nx(extracted_subgraph)
    net.repulsion(node_distance=500, spring_length=300)
    net.show_buttons(filter_=['nodes'])
    # net.show_buttons(filter_=['edges'])
    print(outfile+'.html')
    net.save_graph(outfile+'.html')
