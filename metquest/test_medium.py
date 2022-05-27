from MSI import minimal_media_from_cobrapy

path = '/home/dk/Documents/metGEM_gut/'
outputfilename = '/home/dk/Documents/metGEM_gut/gut_minmedia.txt'
essential_mets = '/home/dk/Documents/metGEM_gut/essential metabolites.txt'

minimal_media_from_cobrapy(path, outputfilename, essential_mets)
