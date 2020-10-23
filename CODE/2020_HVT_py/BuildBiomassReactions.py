from cobra import Reaction
from cobra.core.metabolite import elements_and_molecular_weights

def BuildBiomassReactions(model, params, bof_data):

    elements_and_molecular_weights['R']=0.0
    elements_and_molecular_weights['Z']=0.0

    # Calculate reaction ratios
    bof_dict = {k:v for k,v in bof_data.items() if k.startswith('biomass')}
    obj_rxn = model.reactions.query('bof')[0]

    if sum(bof_dict.values()) > 0:
        bof_ratios = dict()
        for bofcmp,pg in bof_dict.iteritems():
            met_obj = model.metabolites.get_by_id(bofcmp)
            bof_ratios[met_obj] = -1*pg/sum(bof_dict.values())

        # Build BOF
        for m,s in obj_rxn.metabolites.iteritems():
            if m in bof_ratios.keys():
                stoich = s*-1
                tmp_dict = {m:stoich}
                obj_rxn.add_metabolites(tmp_dict)
        for bm_met,bm_s in bof_ratios.iteritems():
            obj_rxn.add_metabolites({bm_met:bm_s})
        print obj_rxn

    for m,s in obj_rxn.metabolites.iteritems():
        if s < 0:
            biomol = []
            if m.id == 'biomass_pigm_h':
                biomol = ['cholphya_h','cholphyc1_h','cholphyc2_h','fxanth_h','diadinx_h','diatox_h','caro_h']
            elif m.id == 'biomass_prot_c':
                biomol = [k for k in bof_data.keys() if k.endswith('trna_c')]
            elif m.id == 'biomass_faa_c':
                biomol = ['ala__L_c','arg__L_c','asp__L_c','cys__L_c','glu__L_c','gly_c','his__L_c','ile__L_c','leu__L_c','lys__L_c','met__L_c','phe__L_c','pro__L_c','ser__L_c','thr__L_c','tyr__L_c','val__L_c']
            elif m.id == 'biomass_mem_lipids_c':
                biomol = [k for k in bof_data.keys() if k.startswith('pc') or k.startswith('pe') or k.startswith('pg') or k.startswith('12dgr') or k.startswith('sqdg') or k.startswith('mgdg') or k.startswith('dgdg')]
            elif m.id == 'biomass_tag_c':
                biomol = [k for k in bof_data.keys() if k.startswith('tag')]
            elif m.id == 'biomass_eps_c':
                biomol = [k for k in bof_data.keys() if k.startswith('udp') or k.startswith('gdp')]
            elif m.id == 'biomass_rna_c':
                biomol = ['atp_c','utp_c','gtp_c','ctp_c']
            bio_dict = {k:v for k,v  in bof_data.items() if k in biomol}

            if sum(bio_dict.values()) > 0:
                bio_ratios = dict() #mmol/g
                if m.id == 'biomass_prot_c':
                    water = model.metabolites.h2o_c
                    for met,pg in bio_dict.iteritems():
                        met_obj = model.metabolites.get_by_id(met)
                        met_obj.formula = met_obj.formula.replace('R','')
                        trna = model.metabolites.get_by_id('trna'+met.replace('trna_c','')+'_c')
                        trna.formula = trna.formula.replace('R','')
                        bio_ratios[met_obj] = -(1e3*pg/sum(bio_dict.values()))/(met_obj.formula_weight-trna.formula_weight-water.formula_weight)
                        bio_ratios[trna] = (1e3*pg/sum(bio_dict.values()))/(met_obj.formula_weight-trna.formula_weight-water.formula_weight)
                        met_obj.formula = met_obj.formula + 'R'
                        trna.formula = trna.formula + 'R'
                    total = sum([-1*v for k, v in bio_ratios.items() if k.id.endswith('trna_c')])
                    bio_ratios[water] = -total*2.
                    atp = model.metabolites.atp_c
                    bio_ratios[atp] = -total*1.
                    gtp = model.metabolites.gtp_c
                    bio_ratios[gtp] = -total*2.
                    proton = model.metabolites.h_c
                    bio_ratios[proton] = total*3.
                    pi = model.metabolites.pi_c
                    bio_ratios[pi] = total*3.
                    adp = model.metabolites.adp_c
                    bio_ratios[adp] = total*1.
                    gdp = model.metabolites.gdp_c
                    bio_ratios[gdp] = total*2.
                elif m.id == 'biomass_eps_c':
                    udp = model.metabolites.udp_c
                    gdp = model.metabolites.gdp_c
                    proton = model.metabolites.h_c
                    for met,pg in bio_dict.iteritems():
                        met_obj = model.metabolites.get_by_id(met)
                        if met.startswith('udp'):
                            bio_ratios[met_obj] = -(1e3*pg/sum(bio_dict.values()))/(met_obj.formula_weight-udp.formula_weight-proton.formula_weight)
                        elif met.startswith('gdp'):
                            bio_ratios[met_obj] = -(1e3*pg/sum(bio_dict.values()))/(met_obj.formula_weight-gdp.formula_weight-proton.formula_weight)
                    bio_ratios[proton] = -1*sum([v for k,v in bio_ratios.items()])
                    bio_ratios[udp] = -1*sum([v for k,v in bio_ratios.items() if k.id.startswith('udp')])
                    bio_ratios[gdp] = -1*sum([v for k,v in bio_ratios.items() if k.id.startswith('gdp')])
                elif m.id == 'biomass_rna_c':
                    for met, pg in bio_dict.iteritems():
                        met_obj = model.metabolites.get_by_id(met)
                        nmp = model.metabolites.get_by_id(met.replace('t','m'))
                        bio_ratios[met_obj] = -(1e3*pg/sum(bio_dict.values()))/nmp.formula_weight
                    total = sum([-1*v for k, v in bio_ratios.items() if k.id in biomol])
                    ppi = model.metabolites.ppi_c
                    bio_ratios[ppi] = total*1.
                else:
                    for met,pg in bio_dict.iteritems():
                        met_obj = model.metabolites.get_by_id(met)
                        bio_ratios[met_obj] = -(1e3*pg/sum(bio_dict.values()))/met_obj.formula_weight

                bio_rxn = [model.reactions.query(r.id)[0] for r in m.reactions if r.id.startswith('biomass')][0]

                for mx,sx in bio_rxn.metabolites.iteritems():
                    if mx in bio_ratios.keys():
                        x = sx*-1
                        tmp_dict = {mx:x}
                        bio_rxn.add_metabolites(tmp_dict)
                for bm_met,bm_s in bio_ratios.iteritems():
                    bio_rxn.add_metabolites({bm_met:bm_s*params['conversion_factor']})

    # Calculate formula and charge of biomass components
    reactions_to_balance = [r for r in model.reactions if r.id.startswith('biomass') and r.id != 'biomass_vit_c']
    reactions_to_balance.append(obj_rxn)

    for bio_rxn in reactions_to_balance:

        met = [m for m in bio_rxn.products if m.id.startswith('biomass')][0]
        met.charge = 0
        met.formula = ''

        imbal = bio_rxn.check_mass_balance()
        if 'charge' in imbal.keys():
            met_charge = imbal['charge']*-1
            del imbal['charge']
        else:
            met_charge = 0

        met_mass = 0
        formula_string = ''

        for e,v in imbal.iteritems():
            if abs(v) > 1e-12:
                mass = elements_and_molecular_weights[e]
                met_mass = met_mass + (mass*-1*v)
                form_str = e + str(-1*v)
                formula_string = formula_string + form_str

        met.formula = formula_string
        met.charge = met_charge
    

    # Check if balanced correctly
    # for r in model.reactions:
    # 	if not r.id.startswith('EX_') and not r.id.startswith('DM_') and not r.id.startswith('sink_'):
    # 		imbal = r.check_mass_balance()
    # 		if len(imbal) > 0:
    # 			print r.id, imbal

    model.repair()

	