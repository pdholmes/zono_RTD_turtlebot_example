# zonoRTD_turtlebot_example
This repo shows how to use RTD with forward reachable sets (FRSs) represented by zonotopes.
We use turtlebot dynamics as our model system.

# Setup Requirements
You must have [CORA_2018](https://tumcps.github.io/CORA/) on your MATLAB path before running this code

# How to Use
Run the script
```matlab
turtlebot_example
```
The script will generate the turtlebot dynamics, compute an FRS, intersect it with an obstacle in workspace, allow the user to choose a trajectory parameter, and plot the slice of the FRS corresponding to the chosen parameter.

# Authors
Patrick Holmes