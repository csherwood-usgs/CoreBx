# CoreBx
Code to analyze North Core Banks DEMs

This is in ../proj/2019_DorianOBX/WayneWright_flights/CoreBx on my desktop. Same place on my laptop.

April 9 - Switched DEM for September from _crop to _v3.

`Uncertainty_calcs.ipynb` - Overall estimates of uncertainty in N Core Banks calcs.

`Analyze_stable_points.ipynb` - Calc. statistics on elevations pulled from GlobalMapper DEMs. Now out of date, because not yet run with _v3.

`Depth_change_histograms.ipynb` - Kind of hypsometric view of changes in profiles. Last revised in Santa Cruz. Needs to be redone with latest calcs.

`CoreBx_funcs.py` - Utility functions called by the notebooks below.

`CoreBx_island_v2.ipynb` - Interpolate DEMS onto rotated grid for entire island. Produces `ncbanks.nc`

`Process_CoreBx_island_v4.ipynb` - Current working version to process `ncbanks.nc` Older versions include `Process_CoreBx_island_v3.ipynb`

`CoreBx_multi_v3.ipynb` - Current version to rotate N. Core Banks DEMs by region. Produces `region_X.nc` for nine regions. Starting with v3, a fill map is generated from a September surface that Andy made. Older versions are: `CoreBx_multi.ipynb` and `CoreBx_multi_v2.ipynb`. After this is run, then run:

`Process_CoreBx_multi_v3.ipynb` - This makes the volume calcs and plots. Older versions: `Process_CoreBx_multi.ipynb` and `Process_CoreBx_multi_v2.ipynb`

`Centroid_test_profile.ipynb` - This uses `profile_a.csv` and makes the multi-plot profile results used in the eposter.

### Other files

`Check_centroid_calcs.ipynb` - Just a test of centroid calcs.

`NC_SeaGrant_profile.ipynb` - Code used to make profile for NC SeaGrant article. Superceded by 

`pybeach_example_notebook.ipynb` - Evaluating pybeach profle tool.

`test_coord_transformations.ipynb` - Example of round-trip UTM to rotated coords.

`Tif_read_write.ipynb` - Trying to round trip tif files...not there yet.

`Process_CoreBx` - First bit of code to process DEMS [delete]

`CoreBx_grid_v1` - Make a grid. Can probably delete.

Date_string_to_plot_example.ipynb
