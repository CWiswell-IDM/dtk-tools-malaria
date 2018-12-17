This is an example use of Amelia's 8-site intervention-testing platform

Built by Jaline Gerardin (jgerardin@idmod.org). Other versions also used by abertozzivilla.

Prereqs:
Runs only on HPC (COMPS).
Must have dtk-tools, dtk-tools-malaria, malaria-toolbox packages installed.
(run "python setup.py" or "python setup.py develop" for all 3 packages)


Files:
intervention_simulation.py : Picks up serialized simulation, runs short simulation with interventions, and calls
analyzers on simulation results. Run as
"python ./intervention_simulation.py"

sweep_functions.py - set up vector params, input file paths for basic simulation configuration, no-net IPs, functions
for sweeping over ITN, IRS, health_seeking, and ATSB interventions.
pfpr_analyzer.py - example analyzer for MalariaSummaryReport_* type output

site_details.csv - name, lat-long, country, and nodeIDs for the 8 nodes
species_details.json - vector bionomics parameters for each species

bin/ - executable and dlls
input/ - demographics and climate files.


Description:
Provides example uses of:
- different vector species in different nodes (demographics file)
- user-defined Individual Properties (demographics file)
- changing an individual's Individual Property upon triggering event (assign_net_ip() in sweep_functions)
- targeting an intervention to individuals with a specific Individual Property (add_annual_ITNs() in sweep_functions)
- setting asset collection id's (intervention_simulation)
- picking up a serialized simulation (intervention_simulation)
- generating a serialized simulation (intervention_simulation)
- case management (add_health_seeking call in sweep_functions)
- bednets (add_ITN call in sweep_functions)
- indoor residual spray (add_IRS call in sweep_functions)
- ATSBs (add_ATSB call in sweep_functions)
- requesting MalariaSummaryReport custom report (intervention_simulation)
- requesting multiple custom reports with descriptive names (intervention_simulation)
- sweeping over multiple random seeds (intervention_simulation)
- sweeping over multiple burnin simulations and setting appropriate parameters (intervention_simulation)
- sweeping over different interventions (intervention_simulation)
- writing an analyzer with arguments (both analyzers)
- writing an analyzer for MalariaSummaryReport custom reports and constructing basic plots of this data (pfpr_analyzer)

