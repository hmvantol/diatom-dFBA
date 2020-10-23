import pandas as pd

def ReadMediumMetsConcen(input_file):
    """
    This function reads list of mets, metNames, EX rates, and molar mass in media_csv as to get the initial
    concentration of the media
    """

    f = pd.read_csv(input_file, header=None, index_col=0)
    f.columns = ['metNames', 'concentrations', 'molar_mass']
    # media = f.iloc[:,:-1].to_dict('index')
    media = f.to_dict('index')

    return media

# print ReadMediumMetsConcen('HVT_DATA/Media/Aquil_thaps.csv')
