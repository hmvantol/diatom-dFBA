import math, cobra
import pandas as pd
import numpy as np
from cobra import Reaction
from cobra.core.metabolite import elements_and_molecular_weights
from Solve import Solve
from BuildBiomassReactions import BuildBiomassReactions
from editPHOArxns import editPHOArxns

def GrowCell(biomass, chnops, chnops_new, target_bofs, ctns, initial_flux, allfluxes, allstatus, params, models, Vmax, km, media, cells_mL):
    '''
    FBA at each time point, then update biomass, update environment.
    '''

    # Initial time point
    t = 0
    i = 0

    # Times when it is dark
    dark_h = []
    count = 0
    for x in range(0,params['totaltime']+params['delta_t'],params['delta_t']):
        count += params['delta_t']
        if count > int(params['L_D'].split(':')[0]): #day first, dawn target
            # if count <= int(params['L_D'].split(':')[0]): #night first, dusk target
            # if count <= int(params['L_D'].split(':')[0])/2 or count > (int(params['L_D'].split(':')[0])/2)*3: # mid-night target
            # if count > int(params['L_D'].split(':')[0])/2 and count <= (int(params['L_D'].split(':')[0])/2)*3: # mid-day target
            dark_h.append(x)
        if count == 24:
            count = 0
    print dark_h

    while t < params['totaltime']:

        if params['verbose'] == 1:
            print '==== Time (0:' + str(params['totaltime']) + ' hr): ' + str(t) + ' hr ===='

        # Set light:dark cycle
        if t in dark_h:
            # Dark
            for j in media.keys():
                if j.startswith('photon'):
                    ctns.loc[i,j] = 0
        else:
            # Light
            for j in media.keys():
                if j.startswith('photon'):
                    ctns.loc[i,j] = 1e6

        for k in range(0,len(params['species'])):

            if 'Thalassiosira' in models[k].id:

                bio_t = biomass[k].iloc[i]

                # Decide what is the growth stage for obj setting
                if t < 1.25*24:
                    phase = 'LE1'
                elif t < 3*24:
                    phase = 'LE2'
                elif t < 7*24:
                    phase = 'ES'
                else:
                    phase = 'MS'

                if t < 1*24:
                    next_m = 1
                elif t < 3*24:
                    next_m = 3
                elif t < 7*24:
                    next_m = 7
                else:
                    next_m = 10

                bof_data = target_bofs[k][phase].to_dict()
                mgchla = bio_t['cholphya_h']*1e3
                gDW = bio_t['total'] #g/L

                if t < 1*24:
                    currPhase = 'ME'
                elif t < 2*24:
                    currPhase = 'LE1'
                elif t < 4*24:
                    currPhase = 'LE2'
                elif t < 8.5*24:
                    currPhase = 'ES'
                else:
                    currPhase = 'MS'

                o2_vals = pd.read_csv(params['input_dir']+'/'+params['oxygen'],index_col=0)
                o2_data = o2_vals[currPhase].to_dict()

                curr_pigment_comp = bio_t[['cholphya_h','cholphyc1_h','cholphyc2_h','fxanth_h','diadinx_h','diatox_h','caro_h']]/(cells_mL[t/24]*1e3) #g/cell
                editPHOArxns(models[k],params,curr_pigment_comp)


                delta_dict = dict()
                bio_dict = dict()
                day_time = [0,3,6,9] # dawn target
                # day_time = [18,21,0,3] # mid-day target
                # day_time = [6,9,12,15] # mid-night target
                # day_time = [12,15,18,20] # dusk target
                for bm,target in bof_data.iteritems():

                    current_mass = bio_t[bm]

                    f_chryso = [1.061728395,1.765432099,2.518518519,3.419753086]
                    f_tag = [1.299945153,1.507099273,1.860080728,2.58428030]
                    f_prot = [1.094315429,1.407776326,1.783453556,2.007579428]
                    f_pigm = [1.061287224,1.205763689,1.456292027,1.554658982]

                    # target_gL = target*1e-12*1e3*cells_mL[(t+24)/24]

                    if t in dark_h:
                        target_gL = target*1e-12*1e3*cells_mL[next_m]
                    else:
                        for x in range(0,4):
                            if t-24*(t/24) == day_time[x]:
                                if bm == 'biomass_chryso_c':
                                    target_gL = f_chryso[x]*target*1e-12*1e3*cells_mL[(t+24)/24]
                                elif bm == 'biomass_tag_c' or bm.startswith('tag'):
                                    target_gL = f_tag[x]*target*1e-12*1e3*cells_mL[(t+24)/24]
                                else:
                                    target_gL = target*1e-12*1e3*cells_mL[next_m]


                    delta = target_gL-current_mass

                    if delta > 0:
                        delta_dict[bm] = delta
                        bio_dict[bm] = 0.
                    elif delta < 0:
                        delta_dict[bm] = 0.
                        bio_dict[bm] = -1*delta
                    else:
                        delta_dict[bm] = 0.
                        bio_dict[bm] = 0.

                    if bm.startswith('biomass') and abs(delta) <= params['error_tol']:
                        delta_dict[bm] = 0.
                        bio_dict[bm] = 0.

                if t in dark_h:
                    delta_dict['biomass_pigm_h'] = 0.
                    delta_dict['biomass_chryso_c'] = 0.
                    delta_dict['biomass_tag_c'] = 0.
                    for bm in bof_data:
                        bio_dict[bm] = bio_dict[bm]/(int(params['L_D'].split(':')[0])/params['delta_t'])
                        delta_dict[bm] = delta_dict[bm]/(int(params['L_D'].split(':')[0])/params['delta_t'])

                obj_rxn = models[k].reactions.query('bof')[0]
                if sum([val for key,val in delta_dict.iteritems() if key.startswith('biomass')]) != 0:
                    models[k].objective = 'DM_biomass_c'
                    BuildBiomassReactions(models[k], params, delta_dict)
                else:
                    models[k].objective = 'ATPM_c'
                    print models[k].reactions.ATPM_c
                    BuildBiomassReactions(models[k], params, delta_dict)

            else: # Not photsynthetic
                mgchla = 1
                gDW = 1
                obj_rxn = models[k].reactions.query('bof')[0]
                models[k].objective = 'DM_biomass_c'


            if 'Thalassiosira' in models[k].id and t not in dark_h:
                aslice = 1./params['slices'] # Divide culture into bins for light attenuation
            else:
                aslice = 1.

            for r in models[k].reactions:
                models[k].reactions.get_by_id(r.id).lower_bound = r.lower_bound*aslice
                models[k].reactions.get_by_id(r.id).upper_bound = r.upper_bound*aslice

            CellDensity = bio_t['total'] # biomass at t(i)
            Xctns = ctns.iloc[i]*params['conversion_factor'] # Ctns at t(i)
            # Max uptake by Michaelis-Menten
            r_metX = (Vmax[k]*(mgchla/gDW)*params['conversion_factor']*Xctns)/(Xctns+km[k]*params['conversion_factor'])
            # Max uptake per cell per unit time (assuming the media and species are well-mixed, both have equal access to media metabolites)
            Xctns_scaled = Xctns/(CellDensity*params['delta_t'])

            # Available internal metabolites
            for x in bio_dict.keys():
                if x not in ['biomass_chitin_c','biomass_polyp_c','biomass_chryso_c','gly_c','atp_c','utp_c','gtp_c','ctp_c','cholphya_h'] and not x.endswith('trna_c') and not x.startswith('tag') and '__L' not in x:
                    del bio_dict[x]
            bio_int = pd.Series(bio_dict)
            for m in bio_int.index:
                met = models[k].metabolites.get_by_id(m)
                met.formula = met.formula.replace('R','')
            gmol = pd.Series(data=[models[k].metabolites.get_by_id(m).formula_weight for m in bio_int.index], index=bio_int.index)
            Xbio = (bio_int/gmol)*1e3*params['conversion_factor']
            Xbio['biomass_chryso_c'] = bio_int['biomass_chryso_c']*6.167*params['conversion_factor']
            Xbio['biomass_chitin_c'] = bio_int['biomass_chitin_c']*4.921*params['conversion_factor']
            Xbio['biomass_polyp_c'] = bio_int['biomass_polyp_c']*12.663*params['conversion_factor']
            gmol['biomass_chryso_c'] = models[k].metabolites.chryso_c.formula_weight
            gmol['biomass_chitin_c'] = models[k].metabolites.chitin_c.formula_weight
            gmol['biomass_polyp_c'] = models[k].metabolites.polyp_c.formula_weight
            Xbio_scaled = Xbio/(CellDensity*params['delta_t'])

            # Photon flux absorption
            if 'Medium light' in models[k].reactions.PHOA410_h.name:
                absorb = {'EX_photon410_e': 4730.075289,'EX_photon430_e': 5817.128965,'EX_photon450_e': 5348.203973,'EX_photon470_e': 4050.000013,'EX_photon490_e': 3464.694801,'EX_photon510_e': 2649.794528,'EX_photon530_e': 1876.490736,'EX_photon550_e': 1334.544022,'EX_photon570_e': 873.4095179,'EX_photon590_e': 740.7816246,'EX_photon610_e': 888.7175101,'EX_photon630_e': 1082.718272,'EX_photon650_e': 1178.924274,'EX_photon670_e': 3322.974688,'EX_photon690_e': 1840.91646}
            elif 'High light' in models[k].reactions.PHOA410_h.name:
                absorb = {'EX_photon410_e': 4299.712311,'EX_photon430_e': 4827.773138,'EX_photon450_e': 4682.224642,'EX_photon470_e': 3843.094021,'EX_photon490_e': 3418.26565,'EX_photon510_e': 2978.620455,'EX_photon530_e': 2481.738844,'EX_photon550_e': 2039.337028,'EX_photon570_e': 1451.311983,'EX_photon590_e': 1255.298421,'EX_photon610_e': 1207.59605,'EX_photon630_e': 1475.056726,'EX_photon650_e': 1549.345022,'EX_photon670_e': 2969.670677,'EX_photon690_e': 1924.104339}
            elif 'Low light' in models[k].reactions.PHOA410_h.name:
                absorb = {'EX_photon410_e': 7732.432316,'EX_photon430_e': 7530.191719,'EX_photon450_e': 6777.205304,'EX_photon470_e': 5420.788645,'EX_photon490_e': 4462.950129,'EX_photon510_e': 3778.897393,'EX_photon530_e': 2988.068035,'EX_photon550_e': 2100.53568,'EX_photon570_e': 1300.702583,'EX_photon590_e': 1075.245665,'EX_photon610_e': 1123.387173,'EX_photon630_e': 1460.196353,'EX_photon650_e': 1580.653401,'EX_photon670_e': 4018.115597,'EX_photon690_e': 2155.737164}
            if params['light_source'] == 'cool_white':
                relArea = {'EX_photon410_e': 0.020538848,'EX_photon430_e': 0.044455157,'EX_photon450_e': 0.032857193,'EX_photon470_e': 0.036698467,'EX_photon490_e': 0.075174979,'EX_photon510_e': 0.040060786,'EX_photon530_e': 0.047049592,'EX_photon550_e': 0.233319199,'EX_photon570_e': 0.070267588,'EX_photon590_e': 0.103796078,'EX_photon610_e': 0.1568836,'EX_photon630_e': 0.0624902,'EX_photon650_e': 0.030806635,'EX_photon670_e': 0.024272336,'EX_photon690_e': 0.021329342}
            elif params['light_source'] == 'halogen':
                relArea = {'EX_photon410_e':0.001262457,'EX_photon430_e':0.004585348,'EX_photon450_e':0.008020995,'EX_photon470_e':0.015208428,'EX_photon490_e':0.032867028,'EX_photon510_e':0.049849784,'EX_photon530_e':0.065038153,'EX_photon550_e':0.079970685,'EX_photon570_e':0.093411994,'EX_photon590_e':0.111565142,'EX_photon610_e':0.115420353,'EX_photon630_e':0.117526779,'EX_photon650_e':0.1119338,'EX_photon670_e':0.103590713,'EX_photon690_e':0.089748342}
            SA = params['surface_area']
            if 'Thalassiosira' in models[k].id:
                pfa = {f:relArea[f]*(60*60/1e7)*(mgchla/gDW)*absorb[f]*params['I']*aslice*params['conversion_factor'] for f in relArea.keys()}

            # Edit the bounds in each slice
            for y in range(0,int(1/aslice)):

                if y == 0:
                    flux1 = pd.Series(0,index=[r.id for r in models[k].reactions])
                    sumfluxes = pd.Series(0,index=[r.id for r in models[k].reactions])
                    list_status = list()
                    photon_lb = {'EX_photon410_e': 0,'EX_photon430_e': 0,'EX_photon450_e': 0,'EX_photon470_e': 0, 'EX_photon490_e': 0,'EX_photon510_e': 0,'EX_photon530_e': 0,'EX_photon550_e': 0,'EX_photon570_e': 0,'EX_photon590_e': 0,'EX_photon610_e': 0,'EX_photon630_e': 0,'EX_photon650_e': 0,'EX_photon670_e': 0,'EX_photon690_e': 0}

                pfd = {f:relArea[f]*(60*60/1e3)*params['I']*(SA/CellDensity)*params['conversion_factor']+photon_lb[f] for f in relArea.keys()}
                I_list = list()

                for j in media.keys():
                    if 'EX_' + j in models[k].reactions:
                        if j.startswith('photon'):
                            if t not in dark_h:
                                if pfd['EX_'+j] < 0:
                                    models[k].reactions.get_by_id('EX_'+j).lower_bound = 0.
                                    models[k].reactions.get_by_id('EX_'+j).upper_bound = 0.
                                    I = 0
                                    I_list.append(I)
                                else:
                                    models[k].reactions.get_by_id('EX_'+j).lower_bound = -1*min(pfa['EX_'+j],pfd['EX_'+j])
                                    models[k].reactions.get_by_id('EX_'+j).upper_bound = -0.9999*min(pfa['EX_'+j],pfd['EX_'+j])
                                    if pfa['EX_'+j] < pfd['EX_'+j]:
                                        I = (pfa['EX_'+j]/(relArea['EX_'+j]*absorb['EX_'+j]))*(1e7/(60*60))*(gDW/mgchla)*(1./(aslice*params['conversion_factor']))
                                        I_list.append(I)
                                    elif pfd['EX_'+j] < pfa['EX_'+j]:
                                        I = (pfd['EX_'+j]/relArea['EX_'+j])*(1e3/(60*60))*(CellDensity/SA)*(1./params['conversion_factor'])
                                        I_list.append(I)
                            else:
                                models[k].reactions.get_by_id('EX_'+j).lower_bound = 0.
                                models[k].reactions.get_by_id('EX_'+j).upper_bound = 0.
                                I = 0
                                I_list.append(I)

                        else:
                            models[k].reactions.get_by_id('EX_'+j).lower_bound = -1*min(Xctns_scaled.loc[j],r_metX.loc[j])*aslice # Take conservative one (closer to zero)

                for j in Xbio_scaled.index:
                    if 'sink_'+j.replace('biomass_','') in models[k].reactions:
                        models[k].reactions.get_by_id('sink_'+j.replace('biomass_','')).lower_bound = -1*Xbio_scaled.loc[j]*aslice
                    if j == 'cholphya_h':
                        models[k].reactions.DM_cholphya_h.lower_bound = Xbio_scaled.loc[j]*aslice


                if 'Thalassiosira' in models[k].id:
                    mean_I = int(sum(I_list)/len(I_list))
                    print 'I =', mean_I
                    if mean_I > 0:
                        if 'High light' in models[k].id:
                            D1 = 1.36247E-08*mean_I+4.94719E-06
                        elif 'Medium light' in models[k].id:
                            D1 = 1.16783E-08*mean_I+4.24045E-06
                        elif 'Low light' in models[k].id:
                            D1 = 6.4751E-09*mean_I+2.35114E-06
                        else:
                            D1 = 5.23188E-09*mean_I+1.89972E-06
                        models[k].reactions.PSII_u.reaction = '2.0 h2o_u + {1} h_h + 4.0 p680_exc_u + {2} pq_u + {0} ps2d1_u --> 4.0 h_u + o2_u + 4.0 p680_u + {2} pqh2_u + {0} ps2d1_exc_u'.format(str(D1),str((2 - D1)*2),str(2 - D1))
                        models[k].reactions.PSII_u.lower_bound = o2_data['gross_lb']*(mgchla/gDW)*params['conversion_factor']*aslice
                        models[k].reactions.PSII_u.upper_bound = o2_data['gross_ub']*(mgchla/gDW)*params['conversion_factor']*aslice
                        models[k].reactions.EX_o2_e.lower_bound = 0.
                        models[k].reactions.EX_o2_e.upper_bound = 10000.*params['conversion_factor']*aslice
                        models[k].reactions.EX_co2_e.upper_bound = 0.
                        models[k].reactions.EX_hco3_e.upper_bound = 0.
                    else:
                        D1 = 0.
                        models[k].reactions.PSII_u.reaction = '2.0 h2o_u + {1} h_h + 4.0 p680_exc_u + {2} pq_u + {0} ps2d1_u --> 4.0 h_u + o2_u + 4.0 p680_u + {2} pqh2_u + {0} ps2d1_exc_u'.format(str(D1),str((2 - D1)*2),str(2 - D1))
                        models[k].reactions.PSII_u.lower_bound = 0.
                        models[k].reactions.PSII_u.upper_bound = 10000.*params['conversion_factor']*aslice
                        # models[k].reactions.EX_o2_e.lower_bound = params['dark_respiration'][0]*(mgchla/gDW)*params['conversion_factor']*aslice
                        # models[k].reactions.EX_o2_e.upper_bound = params['dark_respiration'][1]*(mgchla/gDW)*params['conversion_factor']*aslice
                        models[k].reactions.EX_o2_e.lower_bound = -10000.*params['conversion_factor']*aslice
                        models[k].reactions.EX_o2_e.upper_bound = 0.
                        models[k].reactions.EX_co2_e.upper_bound = 10000.*params['conversion_factor']*aslice
                        models[k].reactions.EX_hco3_e.upper_bound = 10000.*params['conversion_factor']*aslice

                    models[k].reactions.ATPM_c.upper_bound = params['Ic']*params['abs_coeff']*60*60*(mgchla/gDW)*1e-6*params['phi_m']*12.*(3./14)*1e3*params['conversion_factor']*aslice


                else:
                    mean_I = 0.

                # Run FBA on slice
                [fluxes,status] = Solve(models[k], params, t, dark_h, mgchla, gDW)

                if 'Thalassiosira' in models[k].id:
                    photon_lb = {f:photon_lb[f]+models[k].reactions.get_by_id(f).lower_bound for f in ['EX_photon410_e','EX_photon430_e','EX_photon450_e','EX_photon470_e', 'EX_photon490_e','EX_photon510_e','EX_photon530_e','EX_photon550_e','EX_photon570_e','EX_photon590_e','EX_photon610_e','EX_photon630_e','EX_photon650_e','EX_photon670_e','EX_photon690_e']}

                if y == 0: #first slice
                    flux1 = fluxes
                sumfluxes += fluxes
                list_status.append(status)


            # Store results
            initial_flux[k].iloc[i] = flux1
            allfluxes[k].iloc[i] = sumfluxes
            allstatus[params['species'][k]].iloc[i] = list_status

            for r in models[k].reactions:
                models[k].reactions.get_by_id(r.id).upper_bound = r.upper_bound/aslice
                models[k].reactions.get_by_id(r.id).lower_bound = r.lower_bound/aslice

        ## Update environment
        X_next = ctns.iloc[i]
        for k in range(0,len(params['species'])):
            # Update biomass
            mu = allfluxes[k].loc[i,'DM_biomass_c'] #mu1
            total_next = bio_t['total']*np.exp(mu*params['delta_t'])
            slope = (total_next - bio_t['total'])/params['delta_t']

            flex_bm = ['biomass_pigm_h','biomass_prot_c','biomass_faa_c','biomass_mem_lipids_c','biomass_tag_c','biomass_eps_c','biomass_rna_c']

            # Growth
            for bm in bio_t.index:
                if bm in [m.id for m in obj_rxn.metabolites.keys()] and 'bof' in obj_rxn.id:
                    stoich = [-1*x for met,x in obj_rxn.metabolites.iteritems() if met.id == bm][0]
                    bio_rxn = [models[k].reactions.query(r.id)[0] for r in models[k].metabolites.get_by_id(bm).reactions if r.id.startswith('biomass')][0]
                    mu_bio = allfluxes[k].loc[i,bio_rxn.id]
                    if bm in flex_bm:
                        for met,s in bio_rxn.metabolites.iteritems():
                            if s < 0 and met.id in bio_t.index:
                                if met.id in ['atp_c','gtp_c'] and bm == 'biomass_prot_c':
                                    continue
                                else:
                                    if bm == 'biomass_prot_c':
                                        trna = models[k].metabolites.get_by_id('trna'+met.id.replace('trna_c','')+'_c')
                                        water = models[k].metabolites.h2o_c
                                        f = -1*s*1e-3*(met.formula_weight-trna.formula_weight-water.formula_weight)*(1./params['conversion_factor'])
                                    elif bm == 'biomass_eps_c':
                                        udp = models[k].metabolites.udp_c
                                        gdp = models[k].metabolites.gdp_c
                                        proton = models[k].metabolites.h_c
                                        if met.id.startswith('gdp'):
                                            f = -1*s*1e-3*(met.formula_weight-gdp.formula_weight-proton.formula_weight)*(1./params['conversion_factor'])
                                        elif met.id.startswith('udp'):
                                            f = -1*s*1e-3*(met.formula_weight-udp.formula_weight-proton.formula_weight)*(1./params['conversion_factor'])
                                    elif bm == 'biomass_rna_c':
                                        nmp = models[k].metabolites.get_by_id(met.id.replace('t','m'))
                                        f = -1*s*1e-3*nmp.formula_weight*(1./params['conversion_factor'])
                                    else:
                                        f = -1*s*1e-3*met.formula_weight*(1./params['conversion_factor'])
                                    growth = bio_t[met.id]+slope*f*stoich*params['delta_t']
                                    biomass[k].loc[i+1,met.id] = growth
                    else:
                        growth = bio_t[bm]+slope*stoich*params['delta_t']
                        biomass[k].loc[i+1,bm] = growth


            for bm in bio_t.index:
                if biomass[k].loc[i+1,bm] == 0 and bm not in flex_bm:
                    biomass[k].loc[i+1,bm] = bio_t[bm]

            # Mobilization
            dx_sink = [r.id for r in models[k].reactions.query('sink_') if r.id != 'sink_Asn_X_Ser_Thr_c']
            v_sink = allfluxes[k].loc[i,dx_sink]/params['conversion_factor']

            for bm in bio_t.index:
                if 'sink_'+bm.replace('biomass_','') in v_sink.index:
                    sink_rxn = 'sink_'+bm.replace('biomass_','')
                    if bm == 'biomass_chryso_c':
                        mobil = v_sink[sink_rxn]*params['delta_t']*bio_t['total']/6.167
                    elif bm == 'biomass_chitin_c':
                        mobil = v_sink[sink_rxn]*params['delta_t']*bio_t['total']/4.921
                    elif bm == 'biomass_polyp_c':
                        mobil = v_sink[sink_rxn]*params['delta_t']*bio_t['total']/12.663
                    else:
                        mobil = v_sink[sink_rxn]*params['delta_t']*(gmol[bm]/1e3)*bio_t['total']
                    biomass[k].loc[i+1,bm] += mobil

            biomass[k].loc[i+1,'biomass_pigm_h'] = sum(biomass[k].loc[i+1,['cholphya_h','cholphyc1_h','cholphyc2_h','fxanth_h','diadinx_h','diatox_h','caro_h']])
            biomass[k].loc[i+1,'biomass_prot_c'] = sum(biomass[k].loc[i+1,[met for met in bio_t.index if met.endswith('trna_c')]])
            biomass[k].loc[i+1,'biomass_faa_c'] = sum(biomass[k].loc[i+1,[met for met in bio_t.index if '__L_' in met or met == 'gly_c']])
            biomass[k].loc[i+1,'biomass_mem_lipids_c'] = sum(biomass[k].loc[i+1,[met for met in bio_t.index if met.startswith('pc') or met.startswith('pe') or met.startswith('pg') or met.startswith('12dgr') or met.startswith('sqdg') or met.startswith('mgdg') or met.startswith('dgdg')]])
            biomass[k].loc[i+1,'biomass_tag_c'] = sum(biomass[k].loc[i+1,[met for met in bio_t.index if met.startswith('tag')]])
            biomass[k].loc[i+1,'biomass_eps_c'] = sum(biomass[k].loc[i+1,[met for met in bio_t.index if met.startswith('udp') or met.startswith('gdp')]])
            biomass[k].loc[i+1,'biomass_rna_c'] = sum(biomass[k].loc[i+1,['atp_c','utp_c','gtp_c','ctp_c']])
            biomass[k].loc[i+1,'total'] = sum(biomass[k].loc[i+1,[met for met in bio_t.index if met.startswith('biomass')]])

            chnops[k].loc[i+1,:] = chnops[k].loc[i,:]
            chnops_new[k].loc[i+1,:] = chnops_new[k].loc[i,:]
            bio_mets = models[k].metabolites.query('biomass')
            for met in bio_mets:
                if met.id != 'biomass_c':
                    for element in ['C','H','N','O','P','S','Si']:
                        if element in met.elements.keys():
                            bio_acc = met.elements[element]*(biomass[k].loc[i+1,met.id]-bio_t[met.id])/params['conversion_factor'] #mmol/L
                            chnops[k].loc[i+1,element] += bio_acc
                else:
                    for element in ['C','H','N','O','P','S','Si']:
                        if element in met.elements.keys():
                            bio_acc = met.elements[element]*(biomass[k].loc[i+1,'total']-bio_t['total'])/params['conversion_factor'] #mmol/L
                            chnops_new[k].loc[i+1,element] = bio_acc

            # Update all metabolites
            # FBA contribution: ((V*bio(t))/mu)*exp(mu*delta_t-1)
            dx_ex = [r.id for r in models[k].reactions if r.id.startswith('EX_')]
            v = allfluxes[k].loc[i,dx_ex]/params['conversion_factor']

            if mu != 0:
                fba_contrib = v*(bio_t['total']/mu)*(math.exp(mu*params['delta_t'])-1)
            else:
                fba_contrib = v*(bio_t['total']/params['error_tol'])*(math.exp(params['error_tol']*params['delta_t'])-1)
            fba_contrib.index = [x.replace('EX_','') for x in fba_contrib.index]
            X_next = X_next.add(fba_contrib, fill_value=0)


        # Avoid negative concentrations
        for x in range(0,len(X_next)):
            if X_next[x] < 0:
                X_next[x] = 0.

        ctns.iloc[i+1] = X_next
        print 'NO3:', round(ctns.loc[i+1,'no3_e']*1e3,3), 'uM'
        print 'NO2:', round(ctns.loc[i+1,'no2_e']*1e3,3), 'uM'
        print 'NH4:', round(ctns.loc[i+1,'nh4_e']*1e3,3), 'uM'

        # Equalize gas to atmosphere
        if params['system'] == 'mixed':
            ctns.loc[i+1,'co2_e'] = 0.0107215
            ctns.loc[i+1,'hco3_e'] = 1.72733
            ctns.loc[i+1,'o2_e'] = 0.225
        else:
            if t % 24 == 0: #mix once a day
                ctns.loc[i+1,'co2_e'] = 0.0107215
                ctns.loc[i+1,'hco3_e'] = 1.72733
                ctns.loc[i+1,'o2_e'] = 0.225

        if params['verbose'] == 1:
            print 'Biomass (g/L):', biomass[k].loc[i+1,'total']


        # Increase time
        i += 1
        t += params['delta_t']
