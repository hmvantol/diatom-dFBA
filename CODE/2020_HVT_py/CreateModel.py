import cobra
from cobra import Reaction
from BuildBiomassReactions import BuildBiomassReactions


def CreateModel(FBAmodel_path, target_bofs, species, params):
    """
    Load FBA models and define the medium in which these cells were grown.
    """

    # Load models
    models = list()

    for k in range(0,len(species)):

        print "Loading model..."
        models.append(cobra.io.read_sbml_model(FBAmodel_path+species[k]+'.xml'))

        if 'Thaps' in species[k]:
            models[k].id = 'Thalassiosira pseudonana CCMP 1335'

        print models[k].id, "loaded."

        print 'Configuring model...'
        cobra_config = cobra.Configuration()
        cobra_config.lower_bound = -10000.*params['conversion_factor']
        cobra_config.upper_bound = 10000.*params['conversion_factor']
        cobra_config.solver = 'glpk'
        cobra_config.tolerance = params['error_tol']
        print '...done!'

        if 'Thalassiosira' in models[k].id:

            if 'High light' in models[k].reactions.PHOA410_h.name:
                D1 = 1.36247E-08*params['I']+4.94719E-06
            elif 'Medium light' in models[k].reactions.PHOA410_h.name:
                D1 = 1.16783E-08*params['I']+4.24045E-06
            elif 'Low light' in models[k].reactions.PHOA410_h.name:
                D1 = 6.4751E-09*params['I']+2.35114E-06
            else:
                D1 = 5.23188E-09*params['I']+1.89972E-06
            models[k].reactions.PSII_u.reaction = '2.0 h2o_u + {1} h_h + 4.0 p680_exc_u + {2} pq_u + {0} ps2d1_u --> 4.0 h_u + o2_u + 4.0 p680_u + {2} pqh2_u + {0} ps2d1_exc_u'.format(str(D1),str((2 - D1)*2),str(2 - D1))
            models[k].reactions.PSII_u.upper_bound = 10000.

            bof_data = target_bofs[k]['ME'].to_dict()
            mgchla = bof_data['cholphya_h']*1e-9
            gDW = sum([val for bm, val in bof_data.items() if bm.startswith('biomass')])*1e-12

            models[k].reactions.ATPM_c.upper_bound = params['Ic']*params['abs_coeff']*60*60*(mgchla/gDW)*1e-6*params['phi_m']*12.*(3./14)*1e3
            models[k].reactions.CEF_h.upper_bound = params['CEF_ub']
            models[k].remove_reactions(['DM_eps_c'])

            GAM = Reaction('GAM')
            models[k].add_reactions([GAM])
            models[k].reactions.GAM.build_reaction_from_string('GAM_const_c + {0} atp_c + {0} h2o_c --> {0} adp_c + {0} pi_c + {0} h_c'.format(str(params['GAM_val']*params['conversion_factor'])))
            obj_rxn = models[k].reactions.query('bof')[0]
            gam_met = models[k].metabolites.GAM_const_c
            obj_rxn.add_metabolites({gam_met:1})

            reactions_to_convert = models[k].reactions.query('biomass')
            [reactions_to_convert.remove(r) for r in models[k].reactions.query('DM_biomass')]

            for r in reactions_to_convert:
                tmp_dict = dict()
                r.lower_bound = r.lower_bound*params['conversion_factor']
                r.upper_bound = r.upper_bound*params['conversion_factor']
                for m,s in r.metabolites.iteritems():
                    stoich = s*-1
                    tmp_dict[m] = stoich
                    r.add_metabolites({m:stoich})
                for bm_met,bm_s in tmp_dict.iteritems():
                    if not bm_met.id.startswith('biomass'):
                        r.add_metabolites({bm_met:-1*bm_s*params['conversion_factor']})
                    else:
                        r.add_metabolites({bm_met:-1*bm_s})

            BuildBiomassReactions(models[k], params, bof_data)

            # Mobilize biomass
            mets = ['polyp_c','chitin_c','chryso_c','gly_c','atp_c','utp_c','gtp_c','ctp_c','cholphya_h']
            for m in mets:
                rxn = Reaction('sink_'+m)
                models[k].add_reactions([rxn])
                models[k].reactions.get_by_id('sink_'+m).build_reaction_from_string(m + ' <=> ')

            for m in bof_data.keys():
                if m.endswith('trna_c'):
                    rxn = Reaction('sink_'+m)
                    models[k].add_reactions([rxn])
                    if m == 'glytrna_c':
                        models[k].reactions.get_by_id('sink_'+m).build_reaction_from_string(m.replace('trna','') + ' <=> ')
                    else:
                        models[k].reactions.get_by_id('sink_'+m).build_reaction_from_string(m.replace('trna','__L') + ' <=> ')
                if m.startswith('tag') or '__L_' in m:
                    rxn = Reaction('sink_'+m)
                    models[k].add_reactions([rxn])
                    models[k].reactions.get_by_id('sink_'+m).build_reaction_from_string(m + ' <=> ')

            DM_cholphya_h = Reaction('DM_cholphya_h')
            models[k].add_reactions([DM_cholphya_h])
            models[k].reactions.DM_cholphya_h.build_reaction_from_string('cholphya_h --> ')


            # Add preliminary constraints
            for r in models[k].reactions:
                if r.id == 'EX_cncbl3_e' or r.id.startswith('EX_photon'):
                    models[k].reactions.get_by_id(r.id).lower_bound = 0.
                    models[k].reactions.get_by_id(r.id).upper_bound = 0.
                elif r.id.startswith('EX_'):
                    models[k].reactions.get_by_id(r.id).lower_bound = 0.
                    models[k].reactions.get_by_id(r.id).upper_bound = 10000.*params['conversion_factor']
                elif r.id.startswith('sink_'):
                    if r.id == 'sink_Asn_X_Ser_Thr_c':
                        models[k].reactions.get_by_id(r.id).lower_bound = -10000.*params['conversion_factor']
                        models[k].reactions.get_by_id(r.id).upper_bound = 10000.*params['conversion_factor']
                    else:
                        models[k].reactions.get_by_id(r.id).lower_bound = 0.
                        models[k].reactions.get_by_id(r.id).upper_bound = 0.
                elif r.id.startswith('DM_'):
                    models[k].reactions.get_by_id(r.id).lower_bound = 0.
                    models[k].reactions.get_by_id(r.id).upper_bound = 10000.*params['conversion_factor']
                else:
                    models[k].reactions.get_by_id(r.id).upper_bound = r.upper_bound*params['conversion_factor']
                    models[k].reactions.get_by_id(r.id).lower_bound = r.lower_bound*params['conversion_factor']




    return models
	