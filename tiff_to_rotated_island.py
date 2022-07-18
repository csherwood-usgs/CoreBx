# This srcipt does interpolation of one tiff file to island coordinates, supplementing the ones done in CoreBx_island.ipynb 
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import rioxarray
import os.path
from scipy import interpolate
from CoreBx_funcs import *

def tiff_to_rotated_island_nc(yml_fname, tiff_fname, nc_fname):
    r = yaml2dict(yml_fname)
    print(r)

    # Make a grid 
    xu,yu,xrot,yrot,xcoords,ycoords = make_grid(**r)

    ny,nx = np.shape(xu)
    print('Size of grid:',ny,nx)

    dslist=[]

    iswarned = False
    print('Reading from', tiff_fname)

    da = rioxarray.open_rasterio( tiff_fname, masked=True )

    print( np.shape(np.flipud(da['y'].values)), np.shape(da['x'].values), np.shape( np.flipud(da.values)) )
    x = da['x'].values
    y = np.flipud(da['y'].values)

    # Not sure how da.values got a singleton dimension, but squeeze gets rid of it.
    # However, make sure to squeeze before flipping
    z = np.flipud(np.squeeze(da.values))
    print(np.shape(x),np.shape(y),np.shape(z))

    f = interpolate.RegularGridInterpolator( (y, x), z, method='nearest')   

    # Array for interpolated elevations
    zi=np.NaN*np.ones((ny,nx))

    # this is the fast iteration, which only works when all of the source points fall inside the target box
    try:
        zi=f((yu,xu))

        dar = xr.DataArray(zi,dims=['Cross-shore','Alongshore'],coords={'Cross-shore': ycoords, 'Alongshore' :xcoords })

        dar = dar.chunk()
        dslist.append(dar)

        dsa = xr.concat(dslist, dim='map')

        # TODO - Add some metadata to this netcdf file. Add time to the maps. Are the dimensions in the right order?

        dsa.to_netcdf(nc_fname)
        print('Writing to ',nc_fname)

    # this is a slow iteration through all of the points, but allows us to skip ones that are outside
    except:
        print('Could not do fast interpolation.')

    return None



if __name__ == '__main__':
    yml_fname = 'small_island_box.yml'
    drv = 'D:'
    tiff_path = drv+r'\\crs\\proj\\2019_DorianOBX\\Best_files\\dems\\'
    # tiff_fname = tiff_path+'NCB_Sep_EBK_2022-07-13_cog_enclose.tif'
    # tiff_fname = tiff_path+'NCB_Oct_EBK_SfM_lidar_mosiac_cog.tif'
    tiff_fname = tiff_path+'NCB_Nov_EBK_all_2022-07-18_clip_cog_pad.tif'

    # output folder for rotated DEMs
    nc_path =drv+r'\\crs\\proj\\2019_DorianOBX\\Dorian_paper_analyses\\rotated_dems\\'
    # nc_fname = nc_path+'ncorebx_small_NCB_Sep_EBK_2022-07-13_rotated.nc'
    # nc_fname = nc_path+'ncorebx_small_NCB_Oct_EBK_SfM_lidar_mosiac_cog_rotated.nc'
    nc_fname = nc_path+'ncorebx_small_NCB_Nov_EBK_mosaic_rotated.nc'

    if(os.path.isfile(tiff_fname)):
        tiff_to_rotated_island_nc(yml_fname, tiff_fname, nc_fname)
    else:
        print(tiff_fname,'not found.')
