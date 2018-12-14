import copy
import numpy as np


def probability_shifting_module(TM,
                                immune_stimulation_threshold,
                                row_scale_factor,
                                column_scale_factor):
    number_of_bin_edges = max([len(i) for i in TM])
    bin_numbers = np.arange(0, number_of_bin_edges)
    new_TM = copy.copy(TM)
    rows_to_cycle_through = np.arange(immune_stimulation_threshold, max(bin_numbers) + 1)
    # define a row scaling factor that represents the changing force of immunity across current density classes
    for row in rows_to_cycle_through:
        old_column_values = [x for x in TM[row]]
        bin_decreases = [np.log10(column_scale_factor * row_scale_factor * bin)
                         for bin in bin_numbers]
        new_column_values = [old_column_values[i] / bin_decreases[i]
                             for i in range(len(old_column_values))]

        net_zero_column_change = sum(
            [abs(new_column_values[j] - old_column_values[j])
             for j in np.arange(1, len(old_column_values))])
        new_column_values[0] = old_column_values[0] + net_zero_column_change
        new_TM[row] = new_column_values

    return(new_TM)


def set_transition_matrix(cb, TM, scale_factor, immune_stim_threshold=1):
    shifted_TM = probability_shifting_module(TM,
                                             immune_stim_threshold,
                                             row_scale_factor=scale_factor,
                                             column_scale_factor=scale_factor)
    cb.update_params({"Parasite_Peak_Density_Probabilities": shifted_TM,
                      ".Scale_Factor": scale_factor})

    return {"scale_factor": scale_factor}
