import sys, os
import pandas as pd
from CreateModel import CreateModel
from ReadConfig import ReadConfig
from ReadMediumMetsConcen import ReadMediumMetsConcen
from RunCulture import RunCulture

def run_experiments(configfile):
    """
    Run this script to do a simulation.
    """

    params = ReadConfig(configfile)

    # Create output folder
    output_path = params['output_dir']
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    target_bofs = list()
    for k in range(0,len(params['species'])):
        target_bofs.append(pd.read_csv(params['input_dir']+'/'+params['biomass_comp'][k],index_col=0))

    cells_mL = pd.read_csv(params['input_dir']+'/'+params['cellsmL_csv'],squeeze=True,header=None).to_list()

    media = ReadMediumMetsConcen(params['input_dir']+'/'+params['csv_dir']+'/'+params['media_csv'])

    models = CreateModel(params['input_dir']+'/'+params['model_dir']+'/', target_bofs, params['species'], params)

    [biomass, chnops, chnops_new, ctns, initial_flux, allfluxes, allstatus] = RunCulture(params,media,models,cells_mL,target_bofs)

    # Save results
    writer = pd.ExcelWriter(params['output_dir']+'/'+params['output_filename']+'.xlsx', engine='xlsxwriter')
    ctns.to_excel(writer, sheet_name='ctns')
    allstatus.to_excel(writer, sheet_name='status')
    for k in range(0,len(params['species'])):
        biomass[k].to_excel(writer, sheet_name='bio_comps'+str(k))
    for k in range(0,len(params['species'])):
        chnops[k].to_excel(writer, sheet_name='CHNOPS'+str(k))
    for k in range(0,len(params['species'])):
        chnops_new[k].to_excel(writer, sheet_name='CHNOPS_new'+str(k))
    for k in range(0,len(params['species'])):
        allfluxes[k].to_excel(writer, sheet_name='allfluxes'+str(k))
    for k in range(0,len(params['species'])):
        initial_flux[k].to_excel(writer, sheet_name='initial_flux'+str(k))
    pd.DataFrame.from_dict(params,orient='index').to_excel(writer, sheet_name='params')
    writer.save()


from ipdb import launch_ipdb_on_exception

with launch_ipdb_on_exception():
    run_experiments(sys.argv[1])

