# PX3A0 Overhauled

New file structure and code for PX3A0. Code has been broken down by thesis question, where each folder has its own objects and scripts, running independantly to one another but with a lot of reused code.

Documentation is incomplete currently.



## 10x10

Simulation and analysis to demo light distribution over a 10x10 mesh.

`sim.py`: Main simulation script for 10x10 resolution double PDU. Outputs data to `/output_data/`.

`analysis.py`: Analysis script. Outputs basic statistical data and produces several figures to `/figures/`.

`update_and_run.sh`: As simulations run on the server, this copys the latest data set from godzilla and then re-runs the analysis with the most up-to-data dataset.
