
"""
title: analyze_infection_durations.py

description: An analyzer aimed to plotthe distribution of infection durations (total days with patent parasitemia from
a challenge bite) that are emitted via the Malaria Patient Report output from DTK simulations.

author: Jon Russell

date: 11/29/2018

notes and dependencies: uses malariatherapy.txt as input, relies on peak_finding.py in same dir

Institute for Disease Modeling, Bellevue, WA
"""

from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

class DurationsAnalyzer(BaseAnalyzer):
    def __init__(self, weight=1):
        super().__init__(filenames = ['output/MalariaPatientReport.json'])
        self.channel = 'true_asexual_parasites'
        # self.datagroupnames = ['sample', 'sim_id', 'channel']

    # df People = Columns, Days = rows, cells: true parasites
    def select_simulation_data(self, data, simulation):
        patient_data = data[self.filenames[0]]['patient_array']
        patient_df = pd.DataFrame({'patient %d' % x['id']: x[self.channel][0] for x in patient_data})
        # patient_df.reset_index(drop=True)
        return patient_df

    def finalize(self, all_data):
        bins = np.linspace(0,365,12)
        # Looping through sims
        for sim, data in all_data.items():
            sim_id = sim.id[0:3]
            sf = sim.tags.get("scale_factor")
            seed = sim.tags.get("Run_Number")
            figure_name = f"{sim_id}_{sf}_{seed}"
            fig= plt.figure(figure_name)
            ax = fig.gca()
            simulation_durations = []

            # Looping through people
            for column in data.columns:
                (positive_day_indices,) = np.where(np.array(data[column]) >0)
                try:
                    simulation_durations.append(max(positive_day_indices)) # Last day with parasite
                except:
                    print(f"EXCEPTION: ScaleFactor : {sf}, Seed: {seed}")

            ax.hist(simulation_durations,bins = bins, normed=True,alpha= 0.5)
            ax.set_ylim([0,0.01])
            ax.annotate("Mean = %s" %np.mean(simulation_durations),xy = (250,0.009),size= 'large')
            ax.annotate(f"ScaleFactor = {sf}", xy=(250, 0.008), size='large')
            
            ax.set_xlabel('Infection Duration (days)')
            ax.set_ylabel('Frequency')

            plt.savefig(os.path.join(self.working_dir, f"{figure_name}.png"))
        plt.show() #Make sure I see 4


