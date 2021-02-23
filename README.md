# Adaptive Combigrid for GENE-3D

The adaptation algorithm uses SGpp, and the postprocessing needs the GENE-3D python diagnostics. So first install those using the script

```sh
./install_adaptive_combigrid.sh
```

Next up, set up your simulation by providing a template parameter file that has the resolutions parameterized, e.g. with `$nx0` and `$n_procs_x`.
The paths (to the GENE3D installation and to the directory where you want your simulation to run) need to be set right in `adaptation_parameters.sh`, and you need to provide a parameterized batch job submit script for your machine.
See more details in `adaptation_parameters.sh`.
Also, the minimum level vector needs to be set; So far, l=(5,5,5,5,3), that is a resolution of (x=32,y=32,z=32,v=32,w=8) has worked quite well for me.

Run it once on the "minumum" resolution.

```sh
./initialize_adaptive_combigrid.sh
```

Check and verify the result with the GENE diagnostics. If the result looks OK (not unstable), you can spawn more simulations adaptively:

```sh
./run_adaptive_combigrid.sh
```

Have a look at `run_adaptive_scheme_2d.sh` to see how this can be made a batch submission script too (especially useful for dependency chaining on the simulation jobs -- though you currently have to do the dependency chaining yourself). The script will generate the `oldSet.csv` and `allSet.csv` combination scheme files.

And finally, this here extracts the flux profiles and visualizes them:

```sh
./visualize_flux_adaptive_combigrid.sh
```
