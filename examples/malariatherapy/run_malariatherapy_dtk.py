
"""
title: run_malariatherapy_dtk.py

description: An example script for running malariatherapy challenge bite style infections where
infection shapes are drawn using scalable transitions as described in the Malaria 2.0 work.

author: Jon Russell

date: 11/29/2018

notes and dependencies: Uses a special config 'from-cfg.json' in the bin directory to use the
updated immune model

Institute for Disease Modeling, Bellevue, WA
"""

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.generic.climate import set_climate_constant

from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import  ModFn, ModBuilder
from simtools.Analysis.AnalyzeManager import AnalyzeManager

from malaria.interventions.malaria_challenge import add_challenge_trial
from malaria.reports.MalariaReport import add_patient_report

from analyze_infection_durations import DurationsAnalyzer
from immunity_transitions_configuration import set_transition_matrix


# Declare COMPS asset properties
config_filename = "./input/from-cfg.json"
experiment_name = "Malariatherapy_2pt0_infections"
force_immunity = True
debug = False


def set_immune_forcing_builder(transition_matrix=None, scale_factor_array=[2, 5, 10, 100]):
    """Creates an experiment builder setting up immunity based on scale factor"""
    builder = ModBuilder.from_combos(
        [ModFn(set_transition_matrix, transition_matrix, scale_factor)
         for scale_factor in scale_factor_array]
    )
    return builder

# Setup -------------------------------------------------------------------------------------------
cb = DTKConfigBuilder.from_files(config_filename)
cb.update_params({'Vector_Species_Names': [],
                  'Simulation_Duration': 365,
                  'Demographics_Filenames': ['Malariatherapy_demographics.json']
                  })
set_climate_constant(cb)

# Add source of infection (challenge bite or forced EIR) ------------------------------------------
add_challenge_trial(cb, start_day=0)

# ---- CUSTOM REPORTS ----
add_patient_report(cb)
if debug:
    print(f"DEBUG: config builder created")

# Define analyzers, in this case, just the durations analyzer
# show=False will save plots show=True will save and display plots
analyzers = [DurationsAnalyzer(show=True)]

exp_builder = ''
if force_immunity:
    transition_matrix = cb.config['parameters']['Parasite_Peak_Density_Probabilities']
    scale_factor_array = [2, 5, 10, 100]
    exp_builder = set_immune_forcing_builder(transition_matrix, scale_factor_array)
if debug:
    print(f"DEBUG: experiment builder created")


if __name__ == "__main__":
    run_sim_args = {'config_builder': cb,
                    'exp_name': experiment_name,
                    'exp_builder': exp_builder
                    }

    if not SetupParser.initialized:
        SetupParser.init('HPC')

    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
    am = AnalyzeManager(exp_manager.experiment)
    for a in analyzers:
        am.add_analyzer(a)
    am.analyze()
