import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy import interpolate, signal

def stat_summary(x,iprint=False):
    n = len(x)
    nnan = np.sum(np.isnan(x))
    # intitialize with NaNs

    if n > nnan:
        meanx = np.nanmean(x)
        stdx = np.nanstd(x)
        minx = np.nanmin(x)
        d5 = np.nanpercentile(x,5.)
        d25 = np.nanpercentile(x,25.)
        d50 = np.nanpercentile(x,50.)
        d75 = np.nanpercentile(x,75.)
        d95 = np.nanpercentile(x,95.)
        maxx = np.nanmax(x)
    else:
        meanx = np.NaN
        stdx = np.NaN
        minx = np.NaN
        d5   = np.NaN
        d25 = np.NaN
        d50 = np.NaN
        d75 = np.NaN
        d95 = np.NaN
        maxx = np.NaN

    # return it in a dict
    s = {'n':n,'nnan':nnan,'mean':meanx,'std':stdx,'min':minx,'max':maxx,\
         'd5':d5,'d25':d25,'d50':d50,'d75':d75,'d95':d95}
    if iprint:
        for key,value in s.items():
            print('{:6s} = {:.3f}'.format(key,value))

    return s

def analyze_channels(x,diff,dx=1.,vthresh=0.5):
    """
    Calculate channel data from alongshore difference vector

    Input:
        x - vector of alongshore locations
        diff - vector of alongshore elevations (m)
        dx - spacing of points in diff (m)
        vthres - vertical threshold for channel id (m)
        hthresh - horizonal threshold (width) for channel id (m)

        Assumes diff is postive

    """
    # sumple calculation of channel area
    diff[diff <= vthresh]=0.
    chana = np.cumsum(diff)*dx
    dlength = len(diff)*dx
    print('Total channel area m^2/m: {:.2f}'.format(chana[-1]/dlength) )

    nc = 0
    channel_strt = np.array([])
    channel_width = np.array([])
    channel_max_depth = np.array([])
    channel_avg_depth = np.array([])
    channel_area = np.array([])

    run = False
    nc = 0
    for i, z in enumerate(diff):
        if i == 0:
            if z >= vthresh:
                # handle first points
                run=True
                nc = nc+1
                channel_strt = np.append( channel_strt, x[i] )
                channel_width = np.append( channel_width, dx )
                channel_max_depth = np.append( channel_max_depth, z)
                channel_avg_depth = np.append( channel_avg_depth, z)
                channel_area = np.append( channel_area, z*dx )
                channel_sum_depth = z
        else:
            if z >= vthresh and run is False:
                # start new channel
                run = True
                nc = nc+1
                channel_strt = np.append( channel_strt, x[i] )
                channel_width = np.append( channel_width, dx )
                channel_max_depth = np.append( channel_max_depth, z)
                channel_avg_depth = np.append( channel_avg_depth, z)
                channel_area = np.append( channel_area, z*dx )
                channel_sum_depth = z
            elif z >= vthresh and run is True:
                # update existing channel
                run = True
                channel_width[nc-1] = channel_width[nc-1]+dx
                channel_max_depth[nc-1] = np.max( (channel_max_depth[nc-1], z) )
                channel_sum_depth = channel_sum_depth + z
                channel_avg_depth[nc-1] = channel_sum_depth/(channel_width[nc-1]/dx)
                channel_area[nc-1] = channel_avg_depth[nc-1]*channel_width[nc-1]
            elif z <= vthresh:
                # reset
                run = False

    channel_ctr = channel_strt + 0.5*channel_width
    return nc, channel_ctr, channel_area, channel_width, channel_max_depth, channel_avg_depth


def pvol(dist,profs,pfill,dcrest_est,dback,\
    title_str,pnames,imethod='extend',\
    dx = 1.,\
    datum=0.4,\
    maxdist=200.,ztoe=2.4,zowp=1.25,nsmooth=51,
    iverbose=True,iplot=True,iprint=True):
    """
    Calculate cross-sectional volumes for barrier island profiles above datum.
    Assumes distance increases from offshore landward, but plots with ocean to right.

    This is not designed to analyze datum below zero. To do that, fill values that are zeros here should
    be reconsidered...maybe turned into datum.

    Input (lp is length of profiles, nmaps is number of profiles):
        dist(lP) - cross-shore distance (m), starting from arbitrary offshore location, equally spaced at dx
        profs(nmaps, lp) - multiple profiles elevations (m above some datum)
        pfill(lp) - single profile used to fill gaps in other profiles (pre-storm profile)
        dcrest_est - cross-shore location of dune crest (estimated)
        dback - cross-shore location of barrier platform (estimated as 1.25-m contour)
        title_str - string used for title in plots
        pnames - strings with names (dates) of profiles
        imethod -"extend" or "clip" TODO: check clip code...that code is stale
        dx - profile spacing (m)
        datum - elevation used as floor to calculate volumes (m)
        ztoe=2.4 - elevation for estimating dune toe (m)
        maxdist=200.,,zowp=1.25,nsmooth=51,
        iverbose - "True" produces extra output
        iplot - "True" produces plot
        iprint - "True" saves plot

    Returns:
        v, vp, cxcy, zmax, dmax, zcrest, dcrest, zcrest0, dtoe, width_island, width_platform
        width_platform - distance from first point above datum to dback

    """
    # Colors from colorbrewer...but one more than needed so we can skip the first one (too light)
    cols=['#feedde','#fdbe85','#fd8d3c','#e6550d','#a63603']

    nmaps, lp = np.shape(profs)
    if(iverbose):
        print('dx: ',dx)
        print("nmaps, length profiles: ",nmaps,lp)
        print("Shape of dist: ",np.shape(dist))
        print("Shape of profs: ",np.shape(profs))
        print("Shape of pfill: ",np.shape(pfill))
    if(iverbose and iplot):
        fig=plt.figure(figsize=(12,8))
        plt.plot(dist,pfill,':r')
        for i in range(0,nmaps):
            plt.plot(dist,profs[i,:],'-',c=cols[i+1])

    # make a copy of the unchanged profiles for plotting
    profr = profs.copy()

    # find first good value (do this before fitting profile or filling)
    ix = np.zeros((nmaps), dtype=int)
    for i in range(0,nmaps):
        try:
            ix[i] = int(np.argwhere(np.isfinite(profs[i,:]))[0])
            if iverbose:
                print(i,ix[i],profs[i,ix[i]-3:ix[i]+3])
        except:
            # fails because entire profile is NaN
            ix[i] = 0

    # extend the profiles with linear fit or zeros
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
        # truncate the profiles to start at common point (profile w/ least data)
        title_str = title_str+'_clip'
        if iverbose:
            print('clipped')
        imx = int(np.max(ix))
        profs[:,0:imx]=0.

    # determine first point >= datum (do this after fitting profile)
    ixd = np.zeros((nmaps), dtype=int)
    for i in range(0,nmaps):
        try:
            ixd[i] = int(np.argwhere((profs[i,:]>=datum))[0])
            if iverbose:
                print(i,ix[i],profs[i,ixd[i]-3:ixd[i]+3])
        except:
            # fails because entire profile is NaN
            ixd[i] = 0

    # replace NaNs with fill values from September
    for i in range((nmaps)):
        # replace NaNs the fill values
        idx = np.isnan(profs[i,:])
        profs[i,idx]=pfill[idx]

    # replace any other NaNs with zero
    for i in range(0,nmaps):
        profs[i,np.isnan(profs[i,:])]=0.

    # find the back of the island using datum
    iisl = np.zeros((nmaps), dtype=int)
    disl = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        try:
            # find last point >= datum
            #iisl = np.squeeze(np.where(profs[i,int(ix[i]):int(ix[i]+maxdist)]>=datum))[-1]
            iisl[i] = np.squeeze(np.where(profs[i,int(ix[i]):-1]>=datum))[-1]
            disl[i] = dist[int(ix[i]+iisl[i])]
        except:
            pass
        if iverbose:
            print("iisl, disl",iisl[i], disl[i])

    # find the highest point in the profile
    zmax = np.ones((nmaps))*np.nan
    dmax = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        try:
            imxh = int ( np.nanargmax(profs[i,:]) )
            zmax[i] = profs[i,imxh]
            dmax[i] = dist[imxh]
        except:
            pass
        if iverbose:
            print("i, zmax, dmax",i, zmax[i], dmax[i])

    # find highest point within 10 meters of estimated dune crest
    idc = np.ones((nmaps),dtype=int)
    ni = 15
    zcrest0 = np.ones((nmaps))*np.nan
    zcrest = np.ones((nmaps))*np.nan
    dcrest = np.ones((nmaps))*np.nan
    if np.isfinite(dcrest_est) and dcrest_est >= 0:
        idcrest =     int(max(dcrest_est/dx,0.))
        idcrest_min = int(max(idcrest-ni,0))
        idcrest_max = int(min(idcrest+ni,lp))
        if iverbose:
            print('dcrest_est, idcrest: ',dcrest_est, idcrest)
        for i in range((nmaps)):
            try:
                idc[i] = int ( np.nanargmax( profs[i,idcrest_min:idcrest_max]) )
                if i == 0:
                    idc0 = idc[0]
                zcrest[i] = profs[i,idc[i]+idcrest-ni]
                zcrest0[i] = profs[i,idc0+idcrest-ni] # z at location os zmax in first map
                dcrest[i]  = dist[idc[i]+idcrest-ni]
            except:
                pass
            if iverbose:
                print("idc, zcrest, dcrest",idc[i], zcrest[i], dcrest[i])

    # find dune toe as first point >= ztoe
    idt = np.zeros((nmaps), dtype=int)
    dtoe = np.ones((nmaps))*np.nan
    for i in range((nmaps)):
        try:
            # have to squeeze because where returns (0,n). Want first one, so [0]
            idt[i] = np.squeeze(np.where(profs[i,int(ix[i]):int(ix[i]+maxdist)]>=ztoe))[0]
            dtoe[i] = dist[int(ix[i]+idt[i])]
        except:
            pass
    if iverbose:
        print("i, dtoe",idt, dtoe)


    # find the back of the overwash platform using zowp
    # this code is no longer used...back of platform now comes in as dback
    # dowp = np.ones((nmaps))*np.nan
    # for i in range((nmaps)):
    #     # smooth the profile
    #     ps = smooth(np.squeeze(profs[i,:]),nsmooth)
    #     # find last point >zowp
    #     # iowp = np.squeeze(np.where(profs[i,int(ix[i]):int(ix[i]+maxdist)]>=zowp))[-1]
    #     # iowp = np.squeeze(np.where(ps[int(ix[i]):int(ix[i]+maxdist)]>=zowp))[-1]
    #     try:
    #         iowp = np.squeeze(np.where(ps[int(ix[i]):-1]>=zowp))[-1]
    #         dowp[i] = dist[int(ix[i]+iowp)]
    #         if iverbose:
    #             print("i, dowp",i, dowp[i])
    #     except:
    #         if iverbose:
    #             print("i, dowp",i, dowp[i])
    # # if back of platforn is not found, use half the distance from the crest to the back of the island
    # if not np.isfinite(dowp[0]):
    #     if np.isfinite(disl[0]):
    #         dowp[0] = dmax[0]+0.5*(disl[0]-dmax[0])
    #     else:
    #         dowp[0] = dmax[0]


    # calculate total width of Island
    width_island = np.zeros((nmaps))*np.NaN
    width_platform = np.zeros((nmaps))*np.NaN
    for i in range((nmaps)):
        try:
            width_island[i] = disl[i]-dist[ixd[i]]
        except:
            pass
        try:
            # print('dback, ixd[i], dist[ixd[i]]: ',dback, ixd[i], dist[ixd[i]])
            width_platform[i]= dback-dist[ixd[i]]
        except:
            pass

        if iverbose:
            print("width, platform width",width_island[i], width_platform[i])

    # Calculate volumes
    profd = profs.copy()-datum
    profd[np.where(profd<=0.)]=0.
    try:
        v = np.sum(profd,1)*dx
    except:
        v = np.NaN
    try:
        vp = np.sum(profd[:,ixd[i]:int(dback/dx)],1)*dx
    except:
        vp = np.NaN
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
        fig=plt.figure(figsize=(12,8))
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
            plt.plot(dist[ixd[i]],datum,'vr',ms=12)
            plt.plot(dist[ixd[i]],datum,'v',c=cols[i])
        for i in range(0,4):
            plt.plot(dcrest[i],zcrest[i],'^r',ms=12)
            plt.plot(dcrest[i],zcrest[i],'^',c=cols[i])
        plt.plot(dback,zowp,'vr',ms=12)
        for i in range(0,4):
            plt.plot(disl[i],datum,'<y',ms=12)
            plt.plot(disl[i],datum,'<',c=cols[i])
        plt.legend()
        plt.ylim((-1., 6.))
        plt.xlim((lp*dx,0)) # this plots xaxis backwards
        plt.ylabel('Elevation (m NAVD88)')
        plt.xlabel('Across-shore Distance (m)')
        plt.title(title_str)
        if iprint:
            pfn = 'p_'+title_str+'.png'
            plt.savefig(pfn,format='png',dpi=300)

    return v, vp, cxcy, zmax, dmax, zcrest, dcrest, zcrest0, dtoe, width_island, width_platform


def running_mean(y, npts):
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

def running_nanmean(y, npts):
    '''
    Smooth a 1-d array with a moving average
    https://stackoverflow.com/questions/40773275/sliding-standard-deviation-on-a-1d-numpy-array

    Input:
        y - 1-d array
        npts - number of points to average
    Returns:
        ys - smoothed arrays
    '''
    sy = np.ones_like(y)*np.nan
    nrows = y.size - npts + 1
    n = y.strides[0]
    y2D = np.lib.stride_tricks.as_strided(y,shape=(nrows,npts),strides=(n,n))
    nclip = int((npts-1)/2)
    # print(nclip)
    sy[nclip:-nclip] = np.nanmean(y2D,1)
    return sy

def running_stddev(y, npts):
    """
    Smooth a 1-d array w/ moving average of npts
    Return array of smoothed data and moving std. deviation
    https://stackoverflow.com/questions/40773275/sliding-standard-deviation-on-a-1d-numpy-array

    Input:
        y - 1-d array
        npts - number of points to average
    Returns:
        sy -  array of running std. deviation
    """
    sy = np.ones_like(y)*np.nan
    nrows = y.size - npts + 1
    n = y.strides[0]
    y2D = np.lib.stride_tricks.as_strided(y,shape=(nrows,npts),strides=(n,n))
    nclip = int((npts-1)/2)
    # print(nclip)
    sy[nclip:-nclip] = np.nanstd(y2D,1)
    return sy

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

def pcoord(x, y):
    """
    Convert x, y to polar coordinates r, az (geographic convention)
    r,az = pcoord(x, y)
    """
    r  = np.sqrt( x**2 + y**2 )
    az=np.degrees( np.arctan2(x, y) )
    # az[where(az<0.)[0]] += 360.
    az = (az+360.)%360.
    return r, az

def xycoord(r, az):
    """
    Convert r, az [degrees, geographic convention] to rectangular coordinates
    x,y = xycoord(r, az)
    """
    x = r * np.sin(np.radians(az))
    y = r * np.cos(np.radians(az))
    return x, y

def UTM2rot(xutm,yutm,r):
    """
    Convert UTM coordinates to rotated coordinates
    """
    # Convert origin to UTM
    xu,yu = box2UTMh(0.,0.,r['e0'],r['n0'],r['theta'])
    # reverse the calc to find the origin (UTM =0,0) in box coordinates.
    # First, just do the rotation to see where Box = 0,0 falls
    xb0,yb0 = box2UTMh(xu,yu,0.,0.,-r['theta'])
    # Then put in negative values for the offset
    #TODO: why does this return a list of arrays?
    xbl,ybl = box2UTMh(xutm,yutm,-xb0,-yb0,-r['theta'])
    # this fixes it...probably should fix box2UTMh
    xb = np.concatenate(xbl).ravel()
    yb = np.concatenate(ybl).ravel()
    return xb, yb



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
