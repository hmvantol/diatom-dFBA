import numpy as np
import pandas as pd

def ReadVmaxKm(input_file1, input_file2, params, media):
    '''
    This function reads the Vmax and km that were calculated for each model.
    '''

    Vmax = pd.read_csv(input_file1, names=range(0, len(params['species'])), index_col=0)
    km = pd.read_csv(input_file2, names=range(0, len(params['species'])), index_col=0)

    for k in range(0, len(params['species'])):
        for j in Vmax.index:
            if np.isnan(Vmax.loc[j, k]):
                Vmax.loc[j, k] = params['Vmax']
        for j in km.index:
            if np.isnan(km.loc[j, k]):
                km.loc[j, k] = params['km']
        for j in media.keys():
            if j not in Vmax.index:
                Vmax.loc[j, k] = params['Vmax']
            if j not in km.index:
                km.loc[j, k] = params['km']

    return Vmax, km

# print ReadVmaxKm('HVT_DATA/Media/Vmax.csv','HVT_DATA/Media/Km.csv')
