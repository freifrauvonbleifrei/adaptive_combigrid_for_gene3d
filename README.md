# Adaptive Combigrid for GENE-3D

The adaptation algorithm uses SGpp, and the postprocessing needs the GENE-3D python diagnostics. So first install those using the script

```sh
./install_adaptive_combigrid.sh
```

Next up, set up your simulation by providing a template parameter file that has the resolutions parameterized, e.g. with TODO. The paths need to be set right in `sim_launcher.py`, and you need to provide a parameterized batch job submit script for your machine.

Run it once on the "minumum" resolution.

```sh
./initialize_adaptive_combigrid.sh
```

So far, l=(5,5,5,5,3), that is (x=32,y=32,z=32,v=32,w=8) has worked quite well for me.

If the result looks OK (not unstable), you can spawn more simulations adaptively:

```sh
./run_adaptive_combigrid.sh
```

And finally, this here extracts the flux profiles and visualizes them:

```sh
./visualize_flux_adaptive_combigrid.sh
```
