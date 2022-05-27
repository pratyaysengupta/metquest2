import pandas as pd
import os


def get_donors(filename):
    """
    This function gets all the microbes that relieved reactions in their partner, reaction wise.
    Here, donors are the microbes relieving reactions in their partners by donating metabolites.
    :param filename: filename of the relieved reactions file obtained from calculate_pairwiseMSI() function.
    :return: returns a dictionary containing relieved reactions as key and a list of relieving microbes as corresponding
            values.
            {relieved reactions:[donor1, donor2, ...]}
    """
    with open(filename, 'r') as f:
        temp = f.read().strip()
    data = temp.split('\n')
    rxns_collection = []
    for i in range(1, len(data), 3):
        line = data[i].strip().split('\t')
        rxns_collection += line[2:]
    rxns_collection = list(set(rxns_collection))
    donor_dict = {}
    for r in rxns_collection:
        donor_dict[r] = []
        for i in range(1, len(data), 3):
            line = data[i].strip().split('\t')
            if r in line[2:]:
                donor_dict[r].append(line[1])
        donor_dict[r] = list(set(donor_dict[r]))
        donor_dict[r].sort()
    return donor_dict


def write_donors(filename):
    """
    This function writes all the microbes that relieved reactions in their partner, reaction wise.
    Here, donors are the microbes relieving reactions in their partners by donating metabolites.
    :param filename: filename of the relieved reactions file obtained from calculate_pairwiseMSI() function.
    :return: returns nothing. The file will be available in directory same as input file (filename)
    """
    donors = get_donors(filename)
    df1 = pd.DataFrame.from_dict(donors, orient="index")
    df1.T.to_csv(os.path.dirname(filename) + "/donors_of_" +
                 os.path.basename(filename).replace('.tsv', '.csv'))