import copy

from malaria import infection, immunity, symptoms
from dtk.vector.species import set_params_by_species
from malaria.interventions.malaria_drugs import drug_params

# --------------------------------------------------------------
# Malaria disease + drug parameters
# --------------------------------------------------------------

disease_params = {
    "Malaria_Model": "MALARIA_MECHANISTIC_MODEL", 
    "Malaria_Strain_Model": "FALCIPARUM_RANDOM_STRAIN"
}

disease_params.update(infection.params)
disease_params.update(immunity.params)
disease_params.update(symptoms.params)

params = copy.deepcopy(disease_params)
params["PKPD_Model"] = "CONCENTRATION_VERSUS_TIME"
params["Malaria_Drug_Params"] = drug_params
params["Genome_Markers"] = []

set_params_by_species(params, ["arabiensis", "funestus", "gambiae"], "MALARIA_SIM")

# --------------------------------------------------------------
# Innate immunity only
# --------------------------------------------------------------

innate_only = copy.deepcopy(params)
innate_only.update({
    "Antibody_Capacity_Growth_Rate": 0,
    "Max_MSP1_Antibody_Growthrate": 0,
    "Min_Adapted_Response": 0
})
