import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy import interpolate, signal

def box2UTMh(x, y, x0, y0, theta):
    '''
    2D rotation and translation of x, y
    Input:
        x, y - row vectors of original coordinates (must be same size)
        x0, y0 - Offset (location of x, y = (0,0) in new coordinate system)
        theta - Angle of rotation (degrees, CCW from x-axis == Cartesian coorinates)
    Returns:
        xr, yr - rotated, offset coordinates
    '''
    thetar = np.radians(theta)
    c, s = np.cos(thetar), np.sin(thetar)

    # homogenous rotation matrix
    Rh = np.array(((c, -s,  0.),\
                   (s,  c,  0.),\
                   (0., 0., 1.)))
    # homogenous translation matrix
    Th = np.array(((1., 0., x0),\
                   (0., 1., y0),\
                   (0., 0., 1.)))

    # homogenous input x,y
    xyh = np.vstack((x,y,np.ones_like(x)))

    # perform rotation and translation
    xyrh=np.matmul(np.matmul(Th,Rh),xyh)
    xr = xyrh[0,:]
    yr = xyrh[1,:]
    return xr, yr

def map_stats(mp,sfile):
    '''
    Calculate some basic statistics for 3D map arrays
    '''
    mean = np.nanmean(mp,axis=(1,2))
    mad = np.nanmean(np.abs(mp),axis=(1,2))
    dmin = np.nanmin(mp,axis=(1,2))
    dmax = np.nanmax(mp,axis=(1,2))
    rms = np.sqrt(np.nanmean(mp**2.,axis=(1,2)))
    s = np.shape(mp)
    num = []
    numn = []
    for i in range(s[0]):
       num.append(mp[i,:,:].size)
       numn.append(np.count_nonzero(np.isnan(mp[i,:,:])))
    print("Shape: ",s,file=sfile)
    print("mean",mean,file=sfile)
    print("mad",mad,file=sfile)
    print("min",dmin,file=sfile)
    print("max",dmax,file=sfile)
    print("rms",rms,file=sfile)
    print("nans",numn,file=sfile)
    print("size",num,file=sfile)
    return mean, mad

def map_stats2d(mp,sfile):
    '''
    Calculate some basic statistics for 2D map arrays
    '''
    mean = np.nanmean(mp,axis=(0,1))
    mad = np.nanmean(np.abs(mp),axis=(0,1))
    dmin = np.nanmin(mp,axis=(0,1))
    dmax = np.nanmax(mp,axis=(0,1))
    rms = np.sqrt(np.nanmean(mp**2.,axis=(0,1)))
    s = np.shape(mp)
    num = (mp[:,:].size)
    numn = (np.count_nonzero(np.isnan(mp[:,:])))
    print("Shape: ",s,file=sfile)
    print("mean",mean,file=sfile)
    print("mad",mad,file=sfile)
    print("min",dmin,file=sfile)
    print("max",dmax,file=sfile)
    print("rms",rms,file=sfile)
    print("nans",numn,file=sfile)
    print("size",num,file=sfile)
    return mean, mad
