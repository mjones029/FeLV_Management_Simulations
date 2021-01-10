# FeLV_Management_Simulations
Manuscript code for FeLV management simulations

This repository features code and data used in the manuscript "Paradoxes and synergies: optimizing management of a deadly virus in an endangered carnivore." Preprint link coming soon.

**NOTE:** Much of the network and transmission simulation pipeline is adapted from Gilbertson et al 2021 ([BioRxiv preprint](https://doi.org/10.1101/2021.01.09.426055)). The [GitHub repository](https://github.com/mjones029/FIV_FeLV_Transmission) for that manuscript is publicly available.

The analysis pipeline follows the following steps:
1. Parameter set generation
2. Simulations with and withouth managment interventions, including:
    1. No intervention scenario
    2. Proactive vaccination
    3. Reactive vaccination
    4. Reactive test-and-removal
    5. Reactive closure of wildlife highway underpasses
3. Sensitivity analysis

**NOTE:** The majority of analyses and generation of figures occurs within main simulation scripts.


## 1. Parameter set generation
This step includes parameter set generation for all simulations (no interventions and all management interventions). Note that the following script also includes code for generation parameter sets for sensitivity analyses, though this relies in part on results from later simulation steps. Parameter sets are generated in (all scripts can be found in the "Scripts" folder):
1. **Design_parameter_sets.R:** code for generating all parameter sets. Some sensitivity analysis parameter generation relies on results of other simulation scripts. 
    1. "Parameter_sets" folder: Actual parameter sets used in manuscript simulations.

## 2. Simulations with management interventions
Main scripts govern network simulation, transmission simulation (with or without interventions), and analysis steps. Each main script relies on external script functions, as well as a parameters dataset. All main and function scripts can be found in the "Scripts" folder.
1. **Nointer_sims.R:** main script for simulations without interventions. Relies on the following external functions:
    1. ***simulate_pop_ergm.R:*** network simulation function
    2. ***trans_sim_basic.R:*** basic FeLV transmission modeling function
    3. ***post_process_outbreak_data_basic.R:*** processes outbreak data to account for respawning process
    4. ***props affected_births included_basic.R:*** calculates proportions of population in each infection state per time point
2. **Proactive_vax_sims.R:** main script for simulations with **proactive vaccination.** Relies on the following external functions:
    1. ***simulate_pop_ergm.R:*** network simulation function
    2. ***proactive_vax.R:*** randomly assigns proactive vaccination to simulated network
    3. ***trans_sim_basic.R:*** basic FeLV transmission modeling function
    4. ***post_process_outbreak_data_basic.R:*** processes outbreak data to account for respawning process
    5. ***props affected_births included_basic.R:*** calculates proportions of population in each infection state per time point
    6. ***analyze_infections.R:*** calculates basic summary values for infection outcomes of interest
3. **Reactive_vax_sims.R:** main script for simulations with **reactive vaccination.** Relies on the following external functions:
    1. ***simulate_pop_ergm.R:*** network simulation function
    2. ***proactive_vax.R:*** randomly assigns proactive vaccination to simulated network
    3. ***trans_sim_reactvax.R:*** reactive vaccination-specific FeLV transmission simulation function
    4. ***post_process_outbreak_data_basic.R:*** processes outbreak data to account for respawning process
    5. ***props affected_births included_basic.R:*** calculates proportions of population in each infection state per time point
    6. ***analyze_infections.R:*** calculates basic summary values for infection outcomes of interest
    7. ***analyze_reactive_vaccinations.R:*** calculates basic summary values for reactive vaccination outcomes of interest (e.g., total vaccines used)
4. **Reactive_TandR_sims.R:** main script for simulations with **reactive test-and-removal.** Relies on the following external functions:
    1. ***simulate_pop_ergm.R:*** network simulation function
    2. ***proactive_vax.R:*** randomly assigns proactive vaccination to simulated network
    3. ***trans_sim_tandr.R:*** FeLV transmission simulation function specific to reactive test-and-removal scenarios
    4. ***post_process_outbreak_data_tandr.R:*** processes outbreak data to account for respawning process AND test-and-removal
    5. ***props affected_births included_tandr.R:*** calculates proportions of population in each infection state per time point, accounting for addition of test-and-removal process
    6. ***analyze_infections.R:*** calculates basic summary values for infection outcomes of interest
    7. ***analyze_reactive_tandr.R:*** calculates basic summary values for reactive test-and-removal outcomes of interest (e.g., proportion of captures that are "successful")
5. **Reactive_upc_sims.R:** main script for simulations with **wildlife underpass closures.** Relies on the following external functions:
    1. ***simulate_pop_ergm.R:*** network simulation function
    2. ***proactive_vax.R:*** randomly assigns proactive vaccination to simulated network
    3. ***up_edge_id.R:*** identifies network edges between individuals on opposing sides of the I-75 freeway
    4. ***trans_sim_upc.R:*** FeLV transmission simulation function specific to reactive underpass closure scenarios
    5. ***post_process_outbreak_data_basic.R:*** processes outbreak data to account for respawning process
    6. ***props affected_births included_basic.R:*** calculates proportions of population in each infection state per time point
    7. ***analyze_infections.R:*** calculates basic summary values for infection outcomes of interest

## Sensitivity analysis