# Recreating the Data

Two notebooks have been provided: `jet_paper_plots.ipynb` and `supplemental_jet_paper.ipynb`. 
The first notebook has the code that generates the plots for the main manuscript. The second
notebook has the code that generates the plots for the supplemental info. HOWEVER, both notebooks
require that all of the PIC simulations have been run and that their data is available. All of the info
needed to run the PIC simulations is included.

### Running the PIC simulations
The `sim_configs` directory has all of the configuration files needed to run every simulation shown in both the
main manuscript and the supplemental info.

To run a simulation:
1. Ensure [Smilei](https://smileipic.github.io/Smilei/) is properly installed.
2. Clone ALFP and checkout an older commit
   - This is used to generate the laser files that are plotted and also used in the PIC sims
   - `git clone https://github.com/kcharbo3/Arbitrary-Laser-Fields-for-PIC.git`
   - `git checkout 8d661427b0228ff9bbd34663dbeac0431829e800`
3. Choose which simulation you want to run out of those listed in either the `sim_configs/high_res` or
`sim_configs/low_res` directories.
4. Replace the following files in ALFP:
- `fourier_prop/laser_input/advanced_parameters.py` with `sim_configs/RES/SIM_NAME/advanced_parameters.py`
- `fourier_prop/laser_input/laser_parameters.py` with `sim_configs/RES/SIM_NAME/laser_parameters.py`
- `fourier_prop/laser_input/propagation_parameters.py` with `sim_configs/RES/SIM_NAME/propagation_parameters.py`
- `fourier_prop/read_laser/sim_grid_parameters.py` with `sim_configs/RES/SIM_NAME/sim_grid_parameters.py`
- `fourier_prop/sim_helpers/sim_parameters.py` with `sim_configs/RES/SIM_NAME/sim_parameters.py`
5. Generate the laser files using ALFP.
6. Use the `sim_configs/RES/SIM_NAME/namelist.py` as your PIC namelist and run the simulation.

### Running the notebooks
Once the simulations have been run, the first cell of both notebooks need to be modified. Modify all of the 
path constants with their proper values. Then run the notebook.