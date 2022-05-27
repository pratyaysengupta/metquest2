import cobra
import re


def find_transport_rxns(file_names):
    """
    This function finds all the transport and extracellular reactions in all models
    :param file_names: list of GSMM file names
    :return:
        trans_rxn: list of transport and extracellular reactions in all models
    """
    model = []
    for i in range(len(file_names)):
        model.append(cobra.io.read_sbml_model(file_names[i]))
    exc_rxns = model[0].exchanges
    met_id = list(exc_rxns[0].metabolites.keys())[0].id
    pattern = re.compile(r'[_[][a-z]\d*[]]*')
    exc_marker = pattern.findall(met_id)
    trans_rxn = []
    for i in model:
        for rxn in i.reactions:
            for met in rxn.metabolites:
                if met.id.find(exc_marker[-1]) >= 0:
                    trans_rxn.append(rxn.id)
    return trans_rxn, model
