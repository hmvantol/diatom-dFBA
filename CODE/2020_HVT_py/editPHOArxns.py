import pandas as pd
from numpy import trapz

def editPHOArxns(model,params,pigment_comp):

    abs_spec_pigm = pd.read_csv(params['input_dir']+'/'+'abs_spec.csv',index_col=0) # m2/mg pigment

    abs_spec_cell20 = pd.DataFrame(columns=abs_spec_pigm.columns)

    # Convert to cm2/cell and integrate over 20 nm bins
    for bm, col in abs_spec_pigm.iteritems():
        pigm_spec_cell = col.values*pigment_comp[bm]*1e3*1e5 # cm2/cell
        for nm in range(400,700,20):
            abs_spec_cell20.loc[str((nm+nm+20)/2),bm] = trapz(pigm_spec_cell[nm-400:nm+20-400],dx=1)

    # Make equations
    for nm, row in abs_spec_cell20.iterrows():
        model.reactions.get_by_id('PHOA'+nm+'_h').reaction = 'photon'+nm+'_e + ' + ' + '.join([' '.join([str(x/sum(row.values)),y]) for x,y in zip(row.values,abs_spec_cell20.columns)]) + ' --> ' + ' + '.join([' '.join([str(x/sum(row.values)),y[:-1]+'exc_h']) for x,y in zip(row.values,abs_spec_cell20.columns)])