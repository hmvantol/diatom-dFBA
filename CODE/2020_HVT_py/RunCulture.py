import pandas as pd
from ReadVmaxKm import ReadVmaxKm
from GrowCell import GrowCell

def RunCulture(params, media, models, cells_mL, target_bofs):
    """
    This is the core of the dFBA model.
    """
    
    # Initialize
    time_steps = range(0, params['totaltime'], params['delta_t'])

    # Get initial bio comps for each species
    for k in range(0, len(params['species'])):

        biomass = list()
        biomass.append(pd.DataFrame())

        for met in target_bofs[k].index:
            biomass[k][met] = [0.] * (len(time_steps) + 1)
        
        # Small correction
        params['biomass_order'] = sum(
            {b for a, b in target_bofs[k]['ME'].to_dict().iteritems() if a.startswith('biomass')}) * 1e-12 * cells_mL[
                                      0] * 1e3

        rxn = models[k].reactions.query('bof')[0]
        for bm, s in rxn.metabolites.iteritems():
            if s < 0:
                biomass[k].loc[0, bm.id] = -1 * s * params['biomass_order']

        biomass[k].loc[0, 'total'] = sum(biomass[k].iloc[0])

        for bm, s in rxn.metabolites.iteritems():
            if s < 0:
                bio_rxn = [models[k].reactions.query(r.id)[0] for r in bm.reactions if r.id.startswith('biomass')][0]
                for met, stoich in bio_rxn.metabolites.iteritems():
                    if stoich < 0 and met.id in biomass[k].columns:
                        if met.id in ['atp_c', 'gtp_c'] and bm.id == 'biomass_prot_c':
                            continue
                        else:
                            if bm.id == 'biomass_prot_c' and met.id.endswith('trna_c'):
                                met.formula = met.formula.replace('R', '')
                                trna = models[k].metabolites.get_by_id('trna' + met.id.replace('trna_c', '') + '_c')
                                trna.formula = trna.formula.replace('R', '')
                                water = models[k].metabolites.h2o_c
                                biomass[k].loc[0, met.id] = -1 * stoich * biomass[k].loc[0, bm.id] * 1e-3 * (
                                        met.formula_weight - trna.formula_weight - water.formula_weight) * (
                                                                    1. / params['conversion_factor'])
                                met.formula = met.formula + 'R'
                                trna.formula = trna.formula + 'R'
                            elif bm.id == 'biomass_rna_c' and met.id in ['atp_c', 'utp_c', 'gtp_c', 'ctp_c']:
                                nmp = models[k].metabolites.get_by_id(met.id.replace('t', 'm'))
                                biomass[k].loc[0, met.id] = -1 * stoich * biomass[k].loc[
                                    0, bm.id] * 1e-3 * nmp.formula_weight * (1. / params['conversion_factor'])
                            elif bm.id == 'biomass_eps_c':
                                proton = models[k].metabolites.h_c
                                udp = models[k].metabolites.udp_c
                                gdp = models[k].metabolites.gdp_c
                                if met.id.startswith('udp'):
                                    biomass[k].loc[0, met.id] = -1 * stoich * biomass[k].loc[0, bm.id] * 1e-3 * (
                                            met.formula_weight - udp.formula_weight - proton.formula_weight) * (
                                                                        1. / params['conversion_factor'])
                                elif met.id.startswith('gdp'):
                                    biomass[k].loc[0, met.id] = -1 * stoich * biomass[k].loc[0, bm.id] * 1e-3 * (
                                            met.formula_weight - gdp.formula_weight - proton.formula_weight) * (
                                                                        1. / params['conversion_factor'])
                            elif met.id not in ['atp_c', 'gtp_c'] and not met.id.startswith('udp') and not met.id.startswith('gdp'):
                                biomass[k].loc[0, met.id] = -1 * stoich * biomass[k].loc[
                                    0, bm.id] * 1e-3 * met.formula_weight * (1. / params['conversion_factor'])
        
        chnops = list()
        chnops.append(pd.DataFrame())
        for element in ['C', 'H', 'N', 'O', 'P', 'S', 'Si']:
            chnops[k][element] = [0.] * (len(time_steps) + 1)
            chnops[k].loc[0, element] = models[k].metabolites.biomass_c.elements[element] * biomass[k].loc[0, 'total'] / \
                                        params['conversion_factor']  # mmol/L

        chnops_new = list()
        chnops_new.append(pd.DataFrame())
        for element in ['C', 'H', 'N', 'O', 'P', 'S', 'Si']:
            chnops_new[k][element] = [0.] * (len(time_steps) + 1)
            chnops_new[k].loc[0, element] = models[k].metabolites.biomass_c.elements[element] * biomass[k].loc[0, 'total'] / \
                                        params['conversion_factor']  # mmol/L

    # Get inital metabolite concentrations in the media
    ctns = pd.DataFrame()
    for j in media.keys():
        ctns[j] = [0.] * (len(time_steps) + 1)
        ctns[j][0] = media[j]['concentrations']

    # Get Vmax and km
    [Vmax, km] = ReadVmaxKm(params['input_dir'] + '/' + params['csv_dir'] + '/' + params['Vmax_csv'],
                            params['input_dir'] + '/' + params['csv_dir'] + '/' + params['Km_csv'], params, media)

    # Make dataframe of fluxes
    allfluxes = list()
    for k in range(0, len(params['species'])):
        allfluxes.append(pd.DataFrame())
        for r in models[k].reactions:
            allfluxes[k][r.id] = [0.] * len(time_steps)

    initial_flux = list()
    for k in range(0, len(params['species'])):
        initial_flux.append(pd.DataFrame())
        for r in models[k].reactions:
            initial_flux[k][r.id] = [0.] * len(time_steps)

    # Status dataframe
    allstatus = pd.DataFrame()
    for k in params['species']:
        allstatus[k] = [None] * len(time_steps)

    # Run FBA for each time point
    GrowCell(biomass, chnops, chnops_new, target_bofs, ctns, initial_flux, allfluxes, allstatus, params, models, Vmax, km, media,
             cells_mL)

    return biomass, chnops, chnops_new, ctns, initial_flux, allfluxes, allstatus
