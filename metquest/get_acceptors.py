import pandas as pd
import os


def get_acceptors(filename):
    """
    This function gets all microbes in which a reaction is relieved by its partner.
    Here, acceptors are the microbes in which reactions are relieved by accepting metabolites from their donors.
    :param filename: filename of the relieved reactions file obtained from calculate_pairwiseMSI() function.
    :return: returns a dictionary containing relieved reactions as key and a list of acceptor microbes
            as corresponding values.
            {reaction:[acceptor1, acceptor2, ...]}
    """
    with open(filename, 'r') as f:
        temp = f.read().strip()
    data = temp.split('\n')
    rxns_collection = []
    for i in range(1, len(data), 3):
        line = data[i].strip().split(',')
        rxns_collection += line[2:]
    rxns_collection = list(set(rxns_collection))
    acceptor_dict = {}
    for r in rxns_collection:
        acceptor_dict[r] = []
        for i in range(1, len(data), 3):
            line = data[i].strip().split(',')
            if r in line[2:]:
                acceptor_dict[r].append(line[0])
        acceptor_dict[r] = list(set(acceptor_dict[r]))
        acceptor_dict[r].sort()
    return acceptor_dict


def write_acceptors(filename):
    """
    This function writes all microbes in which a reaction is relieved by its partner as csv file.
    Here, acceptors are the microbes in which reactions are relieved by accepting metabolites from their donors.
    :param filename: filename of the relieved reactions file obtained from calculate_pairwiseMSI() function.
    :return: returns nothing. The file will be available in directory same as input file (filename)
    """
    acceptors = get_acceptors(filename)
    df1 = pd.DataFrame.from_dict(acceptors, orient="index")
    df1.T.to_csv(os.path.dirname(filename) + "/acceptors_of_" +
                 os.path.basename(filename).replace('.tsv', '.csv'))

