import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy import interpolate, signal

def pvol(dist,profs,pfill,title_str,pnames,imethod='extend',datum=0.4,iverbose=True,iplot=True,iprint=True):
    """
    """
    # Colors from colorbrewer...but one more than needed so we can skip the first one (too light)
    cols=['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']
    dx = dist[1]-dist[0]
    nmaps, lp = np.shape(profs)
    if(iverbose):
        print('dx: ',dx)
        print("nmaps, length profiles: ",nmaps,lp)

    if(iverbose and iplot):
        fig = plt.figure(figsize=(5,3))
        plt.plot(dist,pfill,':r')
        for i in range(0,nmaps):
            plt.plot(dist,profs[i,:],'-',c=cols[i+1])

    # find first good value
    ix = np.zeros((nmaps))
    for i in range(0,nmaps):
        ix[i] = np.argwhere(np.isfinite(profs[i,:]))[0]
    if(iverbose):
        print('indices to first good value: ',ix)

    # make a copy of the unchanged profiles for plotting
    profr = profs.copy()

    # replace NaNs after this point with the fill values
    for i in range((nmaps)):
        # replace NaNs after first good value with the fill values
        idx = np.isnan(profs[i,:])
        profs[i,idx]=pfill[idx]
        # replace any other NaNs with zero
        profs[i,np.isnan(profs[i,:])]=0.
        # restore values before first good value with NaNs
        profs[i,:int(ix[i])]=np.nan

    if(imethod is 'extend' ):
        title_str = title_str+'_extended'
        if iverbose:
            print('extend')
        npts = int(5/dx)
        # fit a straight line to first 5 points
        for i in range((nmaps)):
            p = np.polyfit( dist[int(ix[i]+1):int(ix[i]+1+npts)],\
                profs[i,int(ix[i]+1):int(ix[i]+1+npts)],1)
            if(p[0]>0.):
                # if slope is positive, replace NaNs with line
                profs[i,0:int(ix[i])]=np.polyval(p,dist[0:int(ix[i])])
            else:
                # if slope is not positive, replace NaNs with zeros
                profs[i,0:int(ix[i])]=0.
        # for i in range((nmaps)):
        #     print(np.sum(np.isnan(profs[i])))
    elif(imethod is 'clip'):
        title_str = title_str+'_clip'
        if iverbose:
            print('clipped')
        imx = int(np.max(ix))
        profs[:,0:imx]=0.
        # for i in range((nmaps)):
        #     print(np.sum(np.isnan(profs[i])))

    # Calculate volumes
    profd = profs.copy()-datum
    profd[np.where(profd<=0.)]=0.
    v = np.sum(profd,1)*dx
    if iverbose:
        print("Profile volumes: ", v)

    # Calculate centroids
    cxcy = np.zeros((nmaps,2))
    profc = profs.copy()
    profc[np.where(profc<=datum)]=np.nan
    for i in range(0,nmaps):
        cxcy[i,0],cxcy[i,1] = centroid(dist,profc[i,:])
    if iverbose:
        print("Centroids: \n",cxcy)

    if iplot:
        fig=plt.figure(figsize=(9,6))
        plt.plot(dist,np.ones_like(dist)*datum,'--',c='dimgray',linewidth=2)
        for i in range(0,4):
            lab = '{0} {1: .0f} m$^3$/m'.format(pnames[i],v[i])
            plt.plot(dist,profr[i,:],'-',linewidth=3,c=cols[i+1],label=lab)
            plt.plot(dist,profs[i,:],':',linewidth=3,c=cols[i+1])
        for i in range(0,4):
            plt.plot(cxcy[i,0],cxcy[i,1],'ok',ms=12)
            plt.plot(cxcy[i,0],cxcy[i,1],'o',c=cols[i])
        plt.legend()
        plt.ylim((-1., 6.))
        plt.xlim((lp*dx,0)) # this plots xaxis backwards
        plt.ylabel('Elevation (m NAVD88)')
        plt.xlabel('Across-shore Distance (m)')
        plt.title(title_str)
        if iprint:
            pfn = 'p_'+title_str+'.svg'
            plt.savefig(pfn)

    return v, cxcy


def centroid(x,z):
    cz = np.nanmean(z)
    cx = np.nansum(z*x)/np.nansum(z)
    return(cx,cz)

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
