import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer


class PrevalenceAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, report_names=["Daily_Report"], sweep_variables=None, working_dir="."):
        super(PrevalenceAnalyzer, self).__init__(working_dir=working_dir,
                                        filenames=["output/MalariaSummaryReport_{name}.json".format(name=name)
                                                      for name in report_names]
                                           )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.data_channel = 'Smeared True PfPR by Parasitemia and Age Bin'
        self.data_channel_type = 'DataByTimeAndPfPRBinsAndAgeBins'
        self.channel_name = 'parasite prevalence'

    def select_simulation_data(self, data, simulation):

        alldata = data[self.filenames[0]]
        age_bins = alldata['Metadata']['Age Bins']
        parasite_bins = alldata['Metadata']['Parasitemia Bins']
        if parasite_bins[-1] > 1e7 :
            parasite_bins[-1] = 1e7
        pr = alldata[self.data_channel_type][self.data_channel][:-1]

        simdata = pd.DataFrame()
        for m, par in enumerate(pr):
            df = pd.DataFrame.from_records(par, columns=age_bins)
            df = df.stack().reset_index()
            df = df.rename(columns={'level_1': 'age',
                                    0: self.channel_name})
            df['density_bin'] = df['level_0'].apply(lambda x: parasite_bins[x])
            del df['level_0']
            df['day'] = m
            simdata = pd.concat([simdata, df])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        df = pd.concat(selected, sort=False).reset_index(drop=True)

        detection = df['density_bin'].unique()[0]

        df = df[df['density_bin'] >= detection]
        df = df.groupby(['age', 'day', 'Coverage', 'Run_Number'])[self.channel_name].agg(np.sum).reset_index()

        num_age_bins = len(df['age'].unique())
        num_coverages = len(df['Coverage'].unique())

        sns.set_style('whitegrid', {'axes.linewidth' : 0.5})
        fig = plt.figure(self.channel_name, figsize=(5*num_age_bins,5))
        fig.subplots_adjust(left=0.05, right=0.98)
        palette = sns.color_palette('husl', num_coverages)

        for ia, (a, adf) in enumerate(df.groupby('age')):
            ax = fig.add_subplot(1,num_age_bins,ia+1)
            gdf = adf.groupby(['day', 'Coverage'])[self.channel_name].agg([np.min, np.max, np.mean]).reset_index()
            for iarm, (arm, armdf) in enumerate(gdf.groupby('Coverage')):
                ax.plot(armdf['day'], armdf['mean'], '-', color=palette[iarm], label=arm)
                ax.fill_between(armdf['day'], armdf['amin'], armdf['amax'], color=palette[iarm], linewidth=0, alpha=0.3)
            ax.set_title('age < %d' % a)
            ax.set_xlabel('days since campaign start')
            ax.set_ylabel('%s with detection limit %.1f' % (self.channel_name, detection))
            ax.set_ylim(0,1)
            if ia == 0 :
                ax.legend(title='SMC coverage')

        plt.show()


if __name__ == "__main__":

    from simtools.Analysis.AnalyzeManager import AnalyzeManager
    from simtools.SetupParser import SetupParser

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    analyzer = PrevalenceAnalyzer(expt_name='single_node_example_with_interventions',
                              sweep_variables=["Run_Number",
                                               "Coverage"
                                               ])

    am = AnalyzeManager('b0986274-e6ff-e811-a2bd-c4346bcb1555', analyzers=analyzer)
    am.analyze()