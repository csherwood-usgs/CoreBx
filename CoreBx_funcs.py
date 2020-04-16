import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy import interpolate, signal

def pvol(dist,profs,pfill,title_str,pnames,imethod='extend',datum=0.4,\
    maxdist=400.,ztoe=2.1,zowp=1.25,nsmooth=51,iverbose=True,iplot=True,iprint=True):
    """
    Calculate cross-sectional volumes for barrier island profiles above datum.
    Assumes distance increases from offshore landward, but plots with ocean to right.

    This is not designed to analyze datum below zero. To do that, fill values that are zeros here should
    be reconsidered...maybe turned into datum.

    Input (lp is length of profiles, nmaps is number of profiles):
        dist(lP) - cross-shore distance (m), starting from arbitrary offshore location
        profs(nmaps, lp) - multiple profiles elevations (m above some datum)
        pfill(lp) - single profile used to fill gaps in other profiles (pre-storm profile)
        title_str - string used for title in plots
        pname - strings with names of profiles
        imethod -"extend" or "clip" TODO: check clip values
        datum - value used to calculate volumes
        iverbose - "True" produces extra output
        iplot - "True" produces plot
        iprint - "True" saves plot
    """
    # Colors from colorbrewer...but one more than needed so we can skip the first one (too light)
    cols=['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']

    dx = dist[1]-dist[0]
    nmaps, lp = np.shape(profs)
    if(iverbose):
        print('dx: ',dx)
        print("nmaps, length profiles: ",nmaps,lp)
        print("Shape of dist: ",np.shape(dist))
        print("Shape of profs: ",np.shape(profs))
        print("Shape of pfill: ",np.shape(pfill))
    if(iverbose and iplot):
        fig = plt.figure(figsize=(5,3))
        plt.plot(dist,pfill,':r')
        for i in range(0,nmaps):
            plt.plot(dist,profs[i,:],'-',c=cols[i+1])

    # make a copy of the unchanged profiles for plotting
    profr = profs.copy()

    # replace NaNs with fill values from September
    for i in range((nmaps)):
        # replace NaNs the fill values
        idx = np.isnan(profs[i,:])
        profs[i,idx]=pfill[idx]

    # find first good value
    ix = np.zeros((nmaps))
    for i in range(0,nmaps):
        try:
            ix[i] = int(np.argwhere(np.isfinite(profs[i,:]))[0])
        except:
            # fails because entire profile is NaN
            ix[i] = 0
    if(iverbose):
        print('indices to first good value, value: ')
        print(ix)

    if(imethod is 'extend' ):
        title_str = title_str+'_extended'
        if iverbose:
            print('extend')
        npts = int(5/dx)
        # fit a straight line to first 5 points
        for i in range((nmaps)):
            try:
                # Not sure why one of these breaks down in
                p = np.polyfit( dist[int(ix[i]+1):int(ix[i]+1+npts)],\
                    profs[i,int(ix[i]+1):int(ix[i]+1+npts)],1)
                if iverbose:
                    print("Slope is: {:.4f}".format(p[0]))
                # if slope is less than 1:50, replace
                if(p[0]>0.02):
                    # if slope is positive, replace NaNs with line
                    profs[i,0:int(ix[i])]=np.polyval(p,dist[0:int(ix[i])])
                else:
                    # if slope is not positive, replace NaNs with zeros
                    profs[i,0:int(ix[i])]=0.

                    # print("warning: replacing slope of {:.4f} with {:.4f}".format(p[0],0.02))
                    # p[0]=0.02
                    # profs[i,0:int(ix[i])]=np.polyval(p,dist[0:int(ix[i])])
            except:
                if iverbose:
                    print('cant calculate slope')
                    print('dist, profs',dist[int(ix[i]+1):int(ix[i]+1+npts)],\
                        profs[i,int(ix[i]+1):int(ix[i]+1+npts)])
                # fill with zeros
                profs[i,0:int(ix[i])]=0.

    elif(imethod is 'clip'):
        title_str = title_str+'_clip'
        if iverbose:
            print('clipped')
        imx = int(np.max(ix))
        profs[:,0:imx]=0.

    for i in range(0,nmaps):
        # replace any other NaNs with zero
        profs[i,np.isnan(profs[i,:])]=0.
        # # restore values before first good value with NaNs
        # profs[i,:int(ix[i])]=np.nan

    # find highest point in first maxdist m
    zmax = np.ones((nmaps))*np.nan
    zmap0 = np.ones((nmaps))*np.nan
    dmax = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        try:
            imxh = np.nanargmax(profs[i,int(ix[i]):int(ix[i]+maxdist)])
            if i == 0:
                ix0 = int(ix[0]+imxh)
            zmax[i] = profs[i,int(ix[i]+imxh)]
            zmap0[i] = profs[i,ix0] # z at location os zmax in first map
            dmax[i] = dist[int(ix[i]+imxh)]
        except:
            pass
        if iverbose:
            print("i, pmax, dmax",i, zmax[i], dmax[i])

    # find dune toe as first point >= ztoe
    dtoe = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        try:
            # have to squeeze because where returns (0,n). Want first one, so [0]
            idt = np.squeeze(np.where(profs[i,int(ix[i]):int(ix[i]+maxdist)]>=ztoe))[0]
            dtoe[i] = dist[int(ix[i]+idt)]
            if iverbose:
                print("i, dtoe",i, dtoe[i])
        except:
            if iverbose:
                print("i, dtoe",i, dtoe[i])

    # find the back of the island using datum
    disl = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        try:
            # find last point >= datum
            #iisl = np.squeeze(np.where(profs[i,int(ix[i]):int(ix[i]+maxdist)]>=datum))[-1]
            iisl = np.squeeze(np.where(profs[i,int(ix[i]):-1]>=datum))[-1]
            disl[i] = dist[int(ix[i]+iisl)]
            if iverbose:
                print("i, disl",i, disl[i])
        except:
            if iverbose:
                print("i, disl",i, disl[i])

    # find the back of the overwash platform using zowp
    dowp = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        # smooth the profile
        ps = smooth(np.squeeze(profs[i,:]),nsmooth)
        # find last point >zowp
        # iowp = np.squeeze(np.where(profs[i,int(ix[i]):int(ix[i]+maxdist)]>=zowp))[-1]
        # iowp = np.squeeze(np.where(ps[int(ix[i]):int(ix[i]+maxdist)]>=zowp))[-1]
        try:
            iowp = np.squeeze(np.where(ps[int(ix[i]):-1]>=zowp))[-1]
            dowp[i] = dist[int(ix[i]+iowp)]
            if iverbose:
                print("i, dowp",i, dowp[i])
        except:
            if iverbose:
                print("i, dowp",i, dowp[i])
    # if back of platforn is not found, use half the distance from the crest to the back of the island
    if not np.isfinite(dowp[0]):
        if np.isfinite(disl[0]):
            dowp[0] = dmax[0]+0.5*(disl[0]-dmax[0])
        else:
            dowp[0] = dmax[0]

    # Calculate volumes
    profd = profs.copy()-datum
    profd[np.where(profd<=0.)]=0.
    v = np.sum(profd,1)*dx
    vp = np.sum(profd[:,0:int(ix[i]+dowp[0])],1)*dx
    if iverbose:
        print("Island volumes: ", v)
        print('Platform volumes:', vp)

    # Calculate centroids
    cxcy = np.zeros((nmaps,2))
    profc = profs.copy()
    profc[np.where(profc<=datum)]=np.nan
    for i in range(0,nmaps):
        try:
            cxcy[i,0],cxcy[i,1] = centroid(dist,profc[i,:])
        except:
            cxcy[i,0],cxcy[i,1] = np.nan, np.nan
    if iverbose:
        print("Centroids: \n",cxcy)

    # nice plot if requested
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
        for i in range(0,4):
            plt.plot(dmax[i],zmax[i],'or',ms=12)
            plt.plot(dmax[i],zmax[i],'o',c=cols[i])
        for i in range(0,4):
            plt.plot(dtoe[i],ztoe,'ob',ms=12)
            plt.plot(dtoe[i],ztoe,'o',c=cols[i])
        for i in range(0,4):
            plt.plot(dowp[i],zowp,'og',ms=12)
            plt.plot(dowp[i],zowp,'o',c=cols[i])
        for i in range(0,4):
            plt.plot(disl[i],datum,'oy',ms=12)
            plt.plot(disl[i],datum,'o',c=cols[i])
        plt.legend()
        plt.ylim((-1., 6.))
        plt.xlim((lp*dx,0)) # this plots xaxis backwards
        plt.ylabel('Elevation (m NAVD88)')
        plt.xlabel('Across-shore Distance (m)')
        plt.title(title_str)
        if iprint:
            pfn = 'p_'+title_str+'.png'
            plt.savefig(pfn,format='png',dpi=300)

    return v, vp, cxcy, zmax, dmax, zmap0, dtoe, dowp

def smooth(y, npts):
    '''
    Smooth a 1-d array with a moving average
    https://stackoverflow.com/questions/20618804/how-to-smooth-a-curve-in-the-right-way

    Input:
        y - 1-d array
        npts - number of points to average
    Returns:
        ys - smoothed arrays
    '''
    box = np.ones(npts)/npts
    ys = np.convolve(y, box, mode='same')
    return ys

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
