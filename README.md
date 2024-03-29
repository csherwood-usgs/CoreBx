# CoreBx
Code to analyze North Core Banks DEMs

This has been moved to ..crs/src ../proj/CoreBx.

December 11 2021 - Working on rebuilding laptop after new disk encryption installed. Renamed default branch from `master` to `main` and created branch `refac`, hoping to clean cobwebs out.

March, 2022 - Merged `refac` back into `main`.  

July, 2022 - Finalizing calcs

### Data

### July analysis
`CoreBx_island.ipynb` - Convert DSM tifs to rotated island coordinates.  

    input:  
* `small_ilsnad_box.yml` - coordinates for extent and location
* various tiff files  


`Core
#### Original Analysis  


`Uncertainty_calcs.ipynb` - Overall estimates of uncertainty in N Core Banks calcs.  

`Analyze_rotated_stable_points_refac.ipynb` - Calc. statistics on elevations for stable points by pulling them from the gridded, rotated DEMs. Replaces earlier versions.  

`crop_std_err` - Use gdal_translate to trim bounding box of std err DEMs to smaller than rotated island bounding box.  

`Depth_change_histograms.ipynb` - Kind of hypsometric view of changes in profiles. Last revised in Santa Cruz. Needs to be redone with latest calcs.

`CoreBx_funcs.py` - Utility functions called by the notebooks below.

`CoreBx_island_island_refac.ipynb` - Interpolate DEMS onto rotated grid for entire island. Produces 3D array of DEMs called `ncorebx_refac.nc`. This version replaces all of the previous `Process` scripts. Switched to rioxarray. No uses yaml file to define rotated coord. system. Tried to correct axis labelling for along- and cross-shore axes.  

`CoreBx_island_island_stderr_refac.ipynb` - Does same thing as above, but for rotated, cropped std err DEMs

`Process_CoreBx_island_refac` - Main analysis program, performed on .nc files with stack of DEMs and stderr maps.

`Centroid_test_profile.ipynb` - This uses `profile_a.csv` and makes the multi-plot profile results used in the eposter.

### Other files

`Check_centroid_calcs.ipynb` - Just a test of centroid calcs.  

`Elevation_transistion_histogram` - Mostly useless joint probability plot of elevation transitions.  

`NC_SeaGrant_profile.ipynb` - Code used to make profile for NC SeaGrant article. Superceded by

`pybeach_example_notebook.ipynb` - Evaluating pybeach profle tool.

`test_coord_transformations.ipynb` - Example of round-trip UTM to rotated coords.  

`test_gdal.ipynb` - Examples of running some GDAL commands with both Python bindings and as subprocesses.

`test_find_dune_toe.ipynb` - This tests a function I wrote for finding dune toe and also runs the pybeach dune toes.  

`Tif_read_write.ipynb` - Trying to round trip tif files...not there yet.

`Process_CoreBx` - First bit of code to process DEMS [delete]

`CoreBx_grid_v1` - Make a grid. Can probably delete.

Date_string_to_plot_example.ipynb
