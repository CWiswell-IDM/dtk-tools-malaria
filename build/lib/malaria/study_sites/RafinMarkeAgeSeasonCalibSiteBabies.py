import logging
import os
import numpy as np
from calibtool.analyzers.Helpers import season_channel_age_density_csv_to_pandas
from calibtool.study_sites.site_setup_functions import \
    config_setup_fn, summary_report_fn, site_input_eir_fn

from calibtool.study_sites.DensityCalibSite import DensityCalibSite

logger = logging.getLogger(__name__)


class RafinMarkeAgeSeasonCalibSiteBabies(DensityCalibSite):
    metadata = {
        'parasitemia_bins': [0.0, 16.0, 70.0, 409.0, np.inf],  # (, 0] (0, 16] ... (409, inf]
        'age_bins': [0, 1, 4, 8, 18, 28, 43, np.inf],  # (, 1] (1, 4] ... (43, inf],
        'seasons': ['DC2', 'DH2', 'W2'],
        'seasons_by_month': {
            'May': 'DH2',
            'September': 'W2',
            'January': 'DC2'
        },
        'village': 'Matsari'
    }

    def get_reference_data(self, reference_type):
        super(RafinMarkeAgeSeasonCalibSiteBabies, self).get_reference_data(reference_type)

        # Load the Parasitology CSV
        dir_path = os.path.dirname(os.path.realpath(__file__))
        reference_csv = os.path.join(dir_path, 'inputs', 'GarkiDB_data', 'GarkiDBparasitology.csv')
        reference_data = season_channel_age_density_csv_to_pandas(reference_csv, self.metadata).reset_index()
        reference_data = (reference_data[reference_data['Age Bin'] == 1.0]).set_index(['Channel', 'Season', 'Age Bin', 'PfPR Bin'])

        return reference_data

    def get_setup_functions(self):
        setup_fns = super(RafinMarkeAgeSeasonCalibSiteBabies, self).get_setup_functions()
        setup_fns.append(config_setup_fn(duration=365 * 2))  # 60 years (with leap years)
        setup_fns.append(summary_report_fn(start=365, interval=365.0 / 12, description='Monthly_Report',
                                           parasitemia_bins=[0.0, 16.0, 70.0, 409.0, 4000000.0],
                                           age_bins=[1.0, 4.0, 8.0, 18.0, 28.0, 43.0, 400000.0]))
        setup_fns.append(site_input_eir_fn(self.name, birth_cohort=True))
        setup_fns.append(lambda cb: cb.update_params(
            {'Demographics_Filenames': ['Calibration\\birth_cohort_demographics_babies.json'],
             'Age_Initialization_Distribution_Type': 'DISTRIBUTION_SIMPLE',
             'Base_Population_Scale_Factor': 10,
             'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',
             "Death_Rate_Dependence": "NONDISEASE_MORTALITY_OFF",
             'Enable_Birth': 1,
             'Enable_Vital_Dynamics': 1,
             'Maternal_Antibodies_Type': 'SIMPLE_WANING',
             }))

        return setup_fns

    def __init__(self):
        super(RafinMarkeAgeSeasonCalibSiteBabies, self).__init__('Rafin_Marke_babies')
