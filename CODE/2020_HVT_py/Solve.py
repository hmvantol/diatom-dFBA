import pandas as pd
import cobra
from cobra.exceptions import Infeasible
from chempy import balance_stoichiometry
from sympy import Symbol


def Solve(model, params, t, dark_h, mgchla, gDW):
    """
    Add model-specific constraints and objectives and solve
    """

    if 'Thalassiosira' in model.id:

        with model:

            aslice = 1. / params['slices']

            PR1 = model.problem.Constraint(
                model.reactions.RUBISO_h.flux_expression - 0.025 * model.reactions.RUBISC_h.flux_expression,
                lb=-10000. * params['conversion_factor'] * aslice,
                ub=0.,
                name='photoresp_ub')

            PR2 = model.problem.Constraint(
                model.reactions.RUBISO_h.flux_expression - 0.001 * model.reactions.RUBISC_h.flux_expression,
                lb=0.,
                ub=10000. * params['conversion_factor'] * aslice,
                name='photoresp_lb')

            PQub = calculatePQ(model.metabolites.biomass_c, 'NO3')
            PQlb = calculatePQ(model.metabolites.biomass_c, 'NH4')

            Q1 = model.problem.Constraint(
                model.reactions.EX_o2_e.flux_expression + params['f'] * PQlb * (
                        model.reactions.EX_hco3_e.flux_expression + model.reactions.EX_co2_e.flux_expression),
                lb=0.,
                ub=10000. * params['conversion_factor'] * aslice)

            Q2 = model.problem.Constraint(
                model.reactions.EX_o2_e.flux_expression + params['f'] * PQub * (
                        model.reactions.EX_hco3_e.flux_expression + model.reactions.EX_co2_e.flux_expression),
                lb=-10000. * params['conversion_factor'] * aslice,
                ub=0.)

            resp_day = model.problem.Constraint(
                (model.reactions.PSII_u.flux_expression - model.reactions.EX_o2_e.flux_expression) - (
                        model.reactions.RUBISO_h.flux_expression + model.reactions.PTOX_h.flux_expression + model.reactions.MEHLER_h.flux_expression + model.reactions.GOX_m.flux_expression + model.reactions.GOX_x.flux_expression + model.reactions.AOX_m.flux_expression + model.reactions.CYOO_m.flux_expression),
                lb=0.,
                ub=0.)

            resp_night = model.problem.Constraint(
                model.reactions.EX_o2_e.flux_expression + (
                        model.reactions.RUBISO_h.flux_expression + model.reactions.PTOX_h.flux_expression + model.reactions.MEHLER_h.flux_expression + model.reactions.GOX_m.flux_expression + model.reactions.GOX_x.flux_expression + model.reactions.AOX_m.flux_expression + model.reactions.CYOO_m.flux_expression),
                lb=0.,
                ub=0.)

            if t not in dark_h:
                if model.reactions.DM_biomass_c.objective_coefficient == 1:
                    model.add_cons_vars([PR1, PR2, resp_day, Q1, Q2])
                    print 'Photosynthetic Quotient (Q):', PQlb, '-', PQub
                else:
                    model.add_cons_vars([PR1, PR2, resp_day])
            else:
                model.add_cons_vars([resp_night])

            if model.reactions.ATPM_c.objective_coefficient != 1.:
                model.reactions.DM_biomass_c.upper_bound = 10000. * params['conversion_factor'] * aslice
            else:
                model.reactions.DM_biomass_c.upper_bound = 0.

            model.reactions.H2Ot_m.lower_bound = 0.
            model.reactions.ITPA_c.upper_bound = 0.

            try:
                pfba = cobra.flux_analysis.parsimonious.pfba(model, fraction_of_optimum=1.)
                # pfba = model.optimize()
                fluxes = pfba.fluxes
                status = pfba.status
                print status, pfba.objective_value
            except Infeasible:
                try:
                    cobra_config = cobra.Configuration()
                    cobra_config.tolerance = params['error_tol']*10.
                    pfba = cobra.flux_analysis.parsimonious.pfba(model, fraction_of_optimum=1.)
                    fluxes = pfba.fluxes
                    status = pfba.status
                    print status, pfba.objective_value
                except Infeasible:
                    if params['verbose'] == 1:
                        print 'Infeasible'
                    fluxes = pd.Series(data=[0.] * len(model.reactions), index=[r.id for r in model.reactions])
                    status = 'infeasible'

            if params['verbose'] == 1:
                for r in fluxes.index:
                    if r.startswith('EX_') and abs(fluxes.loc[r]) >= params['error_tol']:
                        print r, fluxes.loc[r]
                        if fluxes.loc[r] == model.reactions.get_by_id(r).lower_bound and not r.startswith('EX_photon'):
                            print r, 'is limiting'
                    elif r.startswith('sink_') and abs(fluxes.loc[r]) > params['error_tol']:
                        print r, fluxes.loc[r]
                        if fluxes.loc[r] == model.reactions.get_by_id(r).lower_bound:
                            print r, 'is limiting'
                    elif r.startswith('DM_') and abs(fluxes.loc[r]) >= params['error_tol']:
                        print(r, fluxes.loc[r])

    return fluxes, status


def find_symbols_from_dict(dict_of_expressions):
    return_set = set()
    [return_set.update(value.atoms(Symbol)) for value in dict_of_expressions.values()]
    return return_set


def calculatePQ(biomass_met, N_source):
    biomass_formula = ''
    reacs = {'H2O'}
    for e in ['C', 'H', 'N', 'O', 'P', 'S', 'Si']:
        if e in biomass_met.elements.keys() and int(round(biomass_met.elements[e])) > 0:
            biomass_formula = biomass_formula + e + str(int(round(biomass_met.elements[e])))
            if e == 'C':
                reacs.add('CO2')
            elif e == 'P':
                reacs.add('PO4')
            elif e == 'S':
                reacs.add('SO4')
            elif e == 'Si':
                reacs.add('SiO4H4')
            elif e == 'N':
                reacs.add(N_source)

    reac, prod = balance_stoichiometry(reacs, {biomass_formula, 'O2'})

    PQ = (float(prod.get('O2')) / reac.get('CO2'))

    return PQ
