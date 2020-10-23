import ast

def ReadConfig(configfile):
    '''
    Read in a config file. Reads in the parameters and stores them in params object, which will get passed on. Numeric parameters are converted to numbers.
    '''

    params = dict()

    with open(configfile, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                entry = line.strip('\n').split(' = ')
                try:
                    params[entry[0]] = ast.literal_eval(line.strip('\n').replace(entry[0]+' = ',''))
                except ValueError:
                    params[entry[0]] = entry[1]

    return params

# print ReadConfig('HVT_DATA/F2_diatom.txt')


