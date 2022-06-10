# pylint: disable=C0103
import numpy as np
import os
import matplotlib.pyplot as plt
# from scipy import interpolate
from scipy.stats import linregress
from pyproj import Proj, transform
import yaml

def which_computer():
    computername = os.environ['COMPUTERNAME']
    drive = 'C:/'
    if computername == 'IGSAGIEGWSCSH10':
        drive = 'D:/'
    return drive, computername


def beach_slope(dist, prof, zm=0., dx=1., npts=5):
    """Determine whether the slope extends below zm and calculate slope over first npts

    Returns
        zinit - Elevation of first valid point (m), even if not near zm
        slope - Slope of the first npts
    """
    zinit = np.nan
    slope = np.nan

    if(np.all(np.isnan(prof))):

        # bail if the profile is empty
        print('all nans in beach_slope')

    else:

        # find first point with data
        ix = int(np.argwhere(np.isfinite(prof))[0])
        zinit = prof[ix]
        # print(ix, zinit)

        # if first point is <zm, find first point >=zm
        if(prof[ix]<zm):
            ix = np.argmax(prof >= zm)

        # fit a line to the first five points
        # find first five points
        dg = dist[ix:ix+npts]
        zg = prof[ix:ix+npts]
        #print('zg:', zg)

        # discard NaNs
        idx = np.isfinite(dg+zg)
        if(len(idx)<2):
            print('too many NaNs in first npts')

        else:
            try:
                # fit line
                #print('dg:', dg)
                #print('zg:', zg)
                p = np.polyfit(dg[idx], zg[idx], 1)
                slope = p[0]
                # print('slope=',slope)

            except:
                # Not sure if this can ever happen
                print('polyval failed')

    return zinit, slope


def extend_beach(dist, prof, zm=0., slope=0.09, dx=1., npts=7):
    """ Find beginning of profile on seaward side
    If that depth is > zm, extend the profile seaward with specified slope.

    dist - array of cross-shore distances
    prof - array of cross-shore elevations (same size)
    zm - elevation criterion: extend beach if first point is >zm
    slope - slope of extended beach profile
    npts - number of points to fit when extending the profile

    return ix, icross, zinit, del_vol, bp, ext_prof
    """
    ext_prof = prof.copy()
    del_vol = 0.
    bp = 0

    if np.all(np.isnan(prof)):

        # bail if the profile is empty
        print('all nans in extend beach')
        return np.nan, np.nan, np.nan, 0., bp, ext_prof

    else:

        # find first point with data
        ix = int(np.argwhere(np.isfinite(prof))[0])
        zinit = prof[ix]
        ext_prof[0:ix]=np.nan
        #print('ix:, zinit: ',ix, zinit)

        # if second point is less than first, and less than 1m
        # assume first point is bad and advance one
        bp = 0
        while (prof[ix] >= prof[ix+1]) and (prof[ix+1] < 1.5) and (bp < 10):
            prof[ix] = np.nan
            ix += 1
            bp += 1

        # if first point is <zm, NaN out depths < zm
        if(prof[ix]<zm):
            # NaN out deeper than zm
            icross = np.argmax(prof >= zm)

            #print('icross',icross)
            ext_prof[0:icross]=np.nan
            #print('returning with first point')
            return ix, icross, zinit, 0., bp, ext_prof

        # fit a line to the first npts
        else:
            # find first npts
            dg = dist[ix:ix+npts]
            zg = prof[ix:ix+npts]
            # print('zg:', zg)

            # discard NaNs
            idx = np.isfinite(dg+zg)
            ngood = npts - np.sum(np.isnan(dg+zg))
            if ngood < 2:
                print('too many NaNs in first npts',ngood, dg[idx],zg[idx])
                return ix, ix, zinit, 0., bp, ext_prof

            else:
                try:
                    # fit line
                    #print('dg:', dg)
                    #print('zg:', zg)
                    p = np.polyfit(dg[idx], zg[idx], 1)
                    # print('p=',p[0])
                except:
                    # Not sure if this can ever happen
                    print('polyval failed: len, dg, zg:',ngood, dg[idx],zg[idx])
                    return ix, ix, np.nan, 0., bp, ext_prof

                try:
                    # ensure slope is not too shallow
                    if p[0] < slope:
                        p[0] = slope
                        # print('changed slope to ',p[0])

                    # extend seaward with slope
                    x = dist[0:ix]-dist[ix]
                    a = prof[ix]
                    ext_prof[0:ix] = a + slope*x
                    # print('New profile for ix = ',ix)
                    # print(dist[ix-2:ix+2])
                    # print(ext_prof[ix-2:ix+2])

                    # NaN out deeper than zm
                    icross = np.argmax(ext_prof >= zm)

                    #print('icross',icross)
                    ext_prof[0:icross]=np.nan

                    # compute volume added to profile
                    #print(ext_prof[icross-1:ix+1])
                    del_vol = np.sum(ext_prof[icross:ix])*dx
                    #print(del_vol)
                    return ix, icross, zinit, del_vol, bp, ext_prof

                except:
                    # Not sure if this can ever happen
                    print('post-polyval failed: x={},a={}'.format(a,x))
                    return ix, ix, np.nan, 0., bp, ext_prof



    return ix, icross, zinit, del_vol, bp, ext_prof


def find_first_valid(dist, prof, pno):
    """ Find beginning of profile on seaward side

    dist - array of cross-shore distances
    prof - array of cross-shore elevations (same size)
    pno - profile number

    return isy, zshore, bp
    """
    isy = -1
    zshore = np.nan
    bp = -1

    if np.all(np.isnan(prof)):
        # bail if the profile is empty
        print(pno,'all nans in find_island_points')
        return isy, zshore, bp

    else:
        # find first point with data
        isy = int(np.argwhere(np.isfinite(prof))[0])
        bp = 0
        # but skip points with reverse beach slope in first 10 points
        while (prof[isy] >= prof[isy+1]) and (prof[isy+1] < 1.5) and (bp < 10):
            prof[isy] = np.nan
            isy += 1
            bp += 1

        zshore = prof[isy]
    return isy, zshore, bp


def find_dune(dist, prof, isy, idy_guess, pno, zb=.5, ni=20):
    """ Find dune crest
    Should be run on smoothed arrays

    dist - array of cross-shore distances
    prof - array of cross-shore elevations (same size)
    isy - first valid point
    idy_guess - guess at dune crest location
    pno = profile number
    zb = elevation criterion for island back
    ni = number of grid cells around idy_guess to search for dune crest

    return idy, zdune
    """
    idy= idy_guess
    zdune = -99.9

    if np.all(np.isnan(prof)):
        # bail if the profile is empty
        print(pno,'all nans in find_dune')
        return idy, zdune
    else:
        # find highest point within ni grid points of estimated dune crest
        if np.isfinite(idy_guess) and idy_guess >= isy:
            idcrest = int(max(idy_guess, isy))
            idcrest_min = int(max(idcrest-ni, 0))
            idcrest_max = int(min(idcrest+ni, len(prof)-1))
            try:
                idy = int(np.argmax( prof[idcrest_min:idcrest_max]))+idcrest_min
                zdune = prof[idy]
                if(pno==12500):
                    print(idcrest_min, idcrest_max, idcrest, idy)
            except:
                #print(pno,np.sum(np.isnan(prof[idcrest_min:idcrest_max])))
                idy = idy_guess
                zdune = -99.9
        else:
            idy = idy_guess
            zdune = -99.9

        if(pno == 12500):
            plt.plot(dist,prof)
            plt.plot(idy,zdune,'ok')

    return idy, zdune


def find_back(dist, prof, idy, pno, zb=.75):
    """ Find island dback
    Should be run on super-smoothed arrays

    dist - array of cross-shore distances
    prof - array of cross-shore elevations (same size)
    idy - dune crest
    pno = profile number
    zb = elevation criterion for island back

    return iby, zback
    """
    iby = -1
    zback = -99.9

    if np.all(np.isnan(prof)):
        # bail if the profile is empty
        print(pno,'all nans in find_dune_and_back')
        bp = 1
        return iby, zback
    else:
        # find back side
        try:
            iby = np.argwhere(prof[int(idy):-1] >= zb)[-1]
            zback = prof[int(iby)]
        except:

            print(pno,'error idy, iby=',idy, iby)
            iby = idy
            zback = -99.9

    return iby, zback


def find_toe(dist, z, s=0.05, zz=2.4, izero='offshore', debug=False):
    """
    Find the toe of the dune using three algorithms:
      * Most offshore occurrence of z>=ztoe
      * Maximum inflection point
      * Longest chord from reference slope to elevations

    The input profile should not extend much beyond the primary dune line,
    and the highest point should be the primary dune crest.

    Input:
      dist = distance along profile (m; no origin necessary)
      z = elevation profile (m)
      s = slope of reference line (m/m)
      zz = elevation of fixed dune toe (m)
      izero = 'offshore' or 'onshore' - origin of profile: offshore indicates profile distances increase shoreward

    Returns:
      izz - index of toe using fixed elevation method (same as input value; m)
      izc - index of toe using chord method (distance from reference line) (m)
      izip - index of toe using inflection point (m)
      zz - elevation of toe using fixed elevation method (same as input value; m)
      ztoe - elevation of toe using distance from reference line (m)
      zipt - elevation of toe using inflection point (m)

    """

    # how many nans in profile?
    nan_frac = np.sum(np.isnan(z))/len(z)
    if debug:
        print('nan_frac: ', nan_frac)

    # only process if at least 10% data
    if nan_frac < 0.9:
        zsmooth = running_mean(z, 7)

        if izero == 'offshore':  # transect goes from offshore to onshore
            # highest point for offset starting point in offshore segment
            izmax = np.nanargmax(zsmooth)
            zb = zsmooth[0:izmax]
            db = dist[0:izmax]-dist[izmax]

            # chord method: make a sloping reference line
            zr = z[izmax] + s*db
            # find greatest vertical diff between ref line and topography
            try:
                izc = np.nanargmax((zr-zb))
                zc = zb[izc]
                dzc = dist[izc]
            except:
                izc = 0
                zc = np.nan
                dzc = np.nan

            # inflection point method: find minimum in smoothed dz/dx
            dz2dx = np.diff(zb, 2)
            # inflection point at (first) minimum dz2/dx
            try:
                izip = np.nanargmax(dz2dx)+2
                zipt = z[izip]
                dipt = dist[izip]
            except:
                izip = 0
                zipt = np.nan
                dipt = np.nan

            # fixed elevation method: find values greater than zz
            # chose first one headed onshore
            try:
                izz = np.squeeze(np.where(zb >= zz))[0]
                zz = z[izz]
                dzz = dist[izz]
            except:
                izz = 0
                zz = np.nan
                dzz = np.nan

            if debug:
                print('offshore')
                print('izmax=', izmax, 'dist=', dist[izmax], 'elev=', z[izmax])
                print('izc=', izc, 'dist=', dzc, 'elev=', zc)
                print('izip=', izip, 'dist=', dipt, 'elev=', zipt)
                print('izz=', izz, 'dist=', dzz, 'elev=',zz)

                fig,ax = plt.subplots(nrows=1, ncols=3)
                ax[0].plot(dist, z)
                ax[0].plot(dzc, zc, 'o', label='chord')
                ax[0].plot(dipt, zipt, 'o', label='inflec. pt.')
                ax[0].plot(dzz, zz, 'o',label='fixed z.')
                ax[0].legend()

                ax[1].plot(db, zr)
                ax[1].plot(db, zb)
                ax[1].plot(db[izc], zc, 'o')
                # ax[1].invert_xaxis()
                ax[1].set_ylim([-1, 7])

                ax[2].plot(dist[2:izmax], dz2dx)

        else: # transect from onshore to offshore
            # This has not been tested
            # highest point for starting point in offshore segment
            if debug:
                print('onshore')
            izmax = np.argmax(zsmooth)
            zb = zsmooth[0:izmax]
            db = dist[0:izmax]

            # chord method: make a sloping reference line
            zr = z[izmax] - s*np.flip(db)
            # find greatest vertical diff between ref line and topography
            zd = zr-zb
            izc = np.argmax(zd)
            zc = zb[izc]

            # inflection point method: find minimum in dz/dx
            dzb = np.array([0])
            dzb = -np.append(dzb, np.diff(zb)/np.diff(db))
            # inflection point at (first) minimum dz/dx
            izip = np.argmin(dzb)
            zipt = zb[izip]

            # fixed elevation method: find values greater than zz
            iztall = np.squeeze(np.where(zb >= zz))
            # chose last one headed offshore
            izz = iztall[0]

    return izz, izc, izip, zz, zc, zipt


def calcR2(H, T, slope, igflag=0):
    """
    %
    % [R2, S, setup,  Sinc,  SIG,  ir] = calcR2(H, T, slope, igflag);
    %
    % Calculated 2% runup (R2), swash (S), setup (setup), incident swash (Sinc)
    % and infragravity swash (SIG) elevations based on parameterizations from runup paper
    % also Iribarren (ir)
    % August 2010 - Included 15% runup (R16) statistic that, for a Guassian distribution,
    % represents mean+sigma. It is calculated as R16 = setup + swash/4.
    % In a wave tank, Palmsten et al (2010) found this statistic represented initiation of dune erosion.
    %
    %
    % H = significant wave height, reverse shoaled to deep water
    % T = deep-water peak wave period
    % slope = radians
    % igflag = 0 (default)use full equation for all data
    %        = 1  use dissipative-specific calculations when dissipative conditions exist (Iribarren < 0.3)
    %        = 2  use dissipative-specific (IG energy) calculation for all data
    %
    % based on:
    %  Stockdon, H. F., R. A. Holman, P. A. Howd, and J. Sallenger A. H. (2006),
    %    Empirical parameterization of setup, swash, and runup,
    %    Coastal Engineering, 53, 573-588.
    % author: hstockdon@usgs.gov
    # Converted to Python by csherwood@usgs.gov
    """
    g = 9.81

    # make slopes positive!
    slope = np.abs(slope)

    # compute wavelength and Iribarren
    L = (g*T**2) / (2.*np.pi)
    sqHL = np.sqrt(H*L)
    ir = slope/sqHL

    if igflag == 2:                     # use dissipative equations (IG) for ALL data
        R2 = 1.1*(0.039 * sqHL)
        S = 0.046*sqHL
        setup = 0.016*sqHL

    elif igflag == 1 and ir < 0.3:      # if dissipative site use diss equations
        R2 = 1.1*(0.039 * sqHL)
        S = 0.046*sqHL
        setup = 0.016*sqHL

    else:                               # if int/ref site, use full equations
        setup = 0.35*slope*sqHL
        Sinc = 0.75*slope*sqHL
        SIG = 0.06*sqHL
        S = np.sqrt(Sinc**2 + SIG**2)
        R2 = 1.1*(setup + S/2.)
        R16 = 1.1*(setup + S/4.)

    return R2, S, setup, Sinc, SIG, ir, R16


def nanlsfit(x, y):
    """least-squares fit of data with NaNs"""
    ok = ~np.isnan(x) & ~np.isnan(y)
    n = len(ok(bool(ok)))
    xx = x[ok]
    yy = y[ok]
    slope, intercept, r, p, stderr = linregress(xx, yy)
    print("n={}; slope, intercept= {:.4f},{:.4f}; r={:.4f} p={:.4f}, stderr={:.4f} ".format(n, slope, intercept, r, p, stderr))
    return n, slope, intercept, r, p, stderr


def stat_summary(x, iprint=False):
    n = len(x)
    nnan = np.sum(np.isnan(x))
    nvalid = n-nnan
    # intitialize with NaNs

    if n > nnan:
        meanx = np.nanmean(x)
        stdx = np.nanstd(x)
        minx = np.nanmin(x)
        d5 = np.nanpercentile(x, 5.)
        d25 = np.nanpercentile(x, 25.)
        d50 = np.nanpercentile(x, 50.)
        d75 = np.nanpercentile(x, 75.)
        d95 = np.nanpercentile(x, 95.)
        maxx = np.nanmax(x)
    else:
        meanx = np.nan
        stdx = np.nan
        minx = np.nan
        d5 = np.nan
        d25 = np.nan
        d50 = np.nan
        d75 = np.nan
        d95 = np.nan
        maxx = np.nan

    # return it in a dict
    s = {'n':n, 'nnan':nnan, 'nvalid':nvalid, 'mean':meanx, 'std':stdx, 'min':minx, 'max':maxx,
         'd5':d5, 'd25':d25, 'd50':d50, 'd75':d75, 'd95':d95}
    # if iprint:
    #     for key, value in s.items():
    #         print('{:6s} = {:.3f}'.format(key, value)),
    if iprint:
        print("  n, nnan, nvalid: ",s['n'],s['nnan'],s['nvalid'])
        print("  mean, std, min, max   : {:.3f} {:.3f} {:.3f} {:.3f}"
            .format(s['mean'], s['std'], s['min'], s['max']))
        print("  d5, d25, d50, d75, d95: {:.3f} {:.3f} {:.3f} {:.3f} {:.3f}"
            .format(s['d5'], s['d25'], s['d50'], s['d75'], s['d95']))

    return s

def analyze_channels(x, diff, dx=1., vthresh=0.5):
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
    diff[diff <= vthresh] = 0.
    chana = np.cumsum(diff)*dx
    dlength = len(diff)*dx
    print('Total channel area m^2/m: {:.2f}'.format(chana[-1]/dlength))

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
                run = True
                nc = nc+1
                channel_strt = np.append(channel_strt, x[i])
                channel_width = np.append(channel_width, dx)
                channel_max_depth = np.append(channel_max_depth, z)
                channel_avg_depth = np.append(channel_avg_depth, z)
                channel_area = np.append(channel_area, z*dx)
                channel_sum_depth = z
        else:
            if z >= vthresh and run is False:
                # start new channel
                run = True
                nc = nc + 1
                channel_strt = np.append(channel_strt, x[i] )
                channel_width = np.append(channel_width, dx )
                channel_max_depth = np.append(channel_max_depth, z)
                channel_avg_depth = np.append(channel_avg_depth, z)
                channel_area = np.append(channel_area, z*dx)
                channel_sum_depth = z
            elif z >= vthresh and run is True:
                # update existing channel
                run = True
                channel_width[nc-1] = channel_width[nc-1]+dx
                channel_max_depth[nc-1] = np.max( (channel_max_depth[nc-1], z))
                channel_sum_depth = channel_sum_depth + z
                channel_avg_depth[nc-1] = channel_sum_depth/(channel_width[nc-1]/dx)
                channel_area[nc-1] = channel_avg_depth[nc-1]*channel_width[nc-1]
            elif z <= vthresh:
                # reset
                run = False

    channel_ctr = channel_strt + 0.5*channel_width
    return nc, channel_ctr, channel_area, channel_width, channel_max_depth, channel_avg_depth


def find_dune_crest_and_back(ix, odcrest_est, dist, prof, zowp = 1.23):
    """
    """
    # find the back of the island using zowp
    # find last point >= zowp
    iisl = np.squeeze(np.where(prof[int(ix):-1] >= zowp))[-1]
    dback = dist[int(iisl)]

    # find highest point within ni meters of estimated dune crest
    if np.isfinite(idcrest_est) and idcrest_est >= 0:
        idcrest = int(max(dcrest_est, 0.))
        idcrest_min = int(max(idcrest-ni, 0))
        idcrest_max = int(min(idcrest+ni, lp))

        idc = int(np.nanargmax( prof[i,idcrest_min:idcrest_max]))

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


def running_nanmin(y, npts):
    '''
    Smooth a 1-d array with a moving minimum
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
    y2D = np.lib.stride_tricks.as_strided(y, shape=(nrows, npts), strides=(n, n))
    nclip = int((npts-1) / 2)
    # print(nclip)
    sy[nclip:-nclip] = np.nanmin(y2D, 1)
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
    y2D = np.lib.stride_tricks.as_strided(y, shape=(nrows, npts), strides=(n, n))
    nclip = int((npts-1) / 2)
    # print(nclip)
    sy[nclip:-nclip] = np.nanstd(y2D, 1)
    return sy


def centroid(x, z):
    cz = np.nanmean(z)
    cx = np.nansum(z * x) / np.nansum(z)
    return(cx, cz)


def make_grid(name=None, e0=None, n0=None, xlen=None, ylen=None, dxdy=None, theta=None):
    """
    Make a rectangular grid to interpolate elevations onto.
    Takes as argument array of dicts, like:
    r={'name': 'ncorebx_refac', 'e0': 378490., 'n0': 3855740., 'xlen': 36500.0, 'ylen': 1500.0, 'dxdy': 1.0, 'theta': 42.0}
    where:
      e0 - UTM Easting of origin [m]
      n0 - UTM Northing of origin [m]
      xlen - Length of alongshore axis [m]
      ylen - Length of cross-shore axis [m]
      dxdy - grid size (must be isotropic right now) [m]
      theta - rotation CCW from x-axis [deg]
    """
    nx = int((1./dxdy) * xlen)
    ny = int((1./dxdy) * ylen)

    xcoords = np.linspace(0.5*dxdy,xlen-0.5*dxdy,nx)
    ycoords = np.linspace(0.5*dxdy,ylen-0.5*dxdy,ny)

    # these will be the coordinates in rotated space
    xrot, yrot = np.meshgrid(xcoords, ycoords, sparse=False, indexing='xy')

    print('make_grid: Shape of xrot, yrot: ',np.shape(xrot),np.shape(yrot))
    shp = np.shape(xrot)
    xu, yu = box2UTMh(xrot.flatten(), yrot.flatten(), e0, n0, theta)
    xu=np.reshape(xu,shp)
    yu=np.reshape(yu,shp)
    # write the UTM coords of the corners to an ASCII file
    corners = np.asarray(  [[xu[0][0],yu[0][0]],
                           [xu[0][-1],yu[0][-1]],
                           [xu[-1][-1],yu[-1][-1]],
                           [xu[-1][0],yu[-1][0]],
                           [xu[0][0],yu[0][0]]])

    print('corners x, corners y]')
    print(corners)
    fn_corners = name+'.csv'
    print('Saving to '+fn_corners)
    np.savetxt(fn_corners, corners, delimiter=",")
    return xu, yu, xrot, yrot, xcoords, ycoords


def box2UTMh(x, y, x0, y0, theta):
    '''
    2D rotation and translation of x, y
    Input:
        x, y - row vectors of original coordinates (must be same size)
        x0, y0 - Offset (location of x, y = (0,0) in new coordinate system)
        theta - Angle of rotation (degrees, CCW from x-axis == Cartesian coorinates)
    Returns:
        x_r, y_r - rotated, offset coordinates
    '''
    thetar = np.radians(theta)
    c, s = np.cos(thetar), np.sin(thetar)

    # homogenous rotation matrix
    Rh = np.array(((c, -s,  0.),
                   (s,  c,  0.),
                   (0., 0., 1.)))
    # homogenous translation matrix
    Th = np.array(((1., 0., x0),
                   (0., 1., y0),
                   (0., 0., 1.)))

    # homogenous input x,y
    xyh = np.vstack((x,y,np.ones_like(x)))

    # perform rotation and translation
    xyrh = np.matmul(np.matmul(Th,Rh), xyh)
    x_r = xyrh[0,:]
    y_r = xyrh[1,:]
    return x_r, y_r


def pcoord(x, y):
    """
    Convert x, y to polar coordinates r, az (geographic convention)
    r,az = pcoord(x, y)
    """
    r = np.sqrt(x**2 + y**2)
    az = np.degrees(np.arctan2(x, y))
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


def UTM2Island(eutm, nutm, eoff=378489.45785127, noff=3855740.50113774, rot=42.):
    """
    Convert UTM NAD83 Zone 18N easting, northing to N. Core Banks alongshore, cross-shore coordinates
    xisl, yisl = UTM2Island( eutm, nutm )
    """
    [r, az] = pcoord(eutm-eoff, nutm-noff)
    az = az + rot
    [xisl,yisl] = xycoord(r,az)
    return xisl, yisl


def island2UTM(alongshore, across_shore, eoff=378489.45785127, noff=3855740.50113774, rot=42.):
    """Convert island coordinates to UTM
       Inverse of UTM2Island()
       Better to use values from the dict than defaults for translation/rotation values

       Here is code for UTM2island:
          [r, az] = pcoord(eutm-eoff, nutm-noff)
          az = az + rot
          [xisl,yisl] = xycoord(r,az)
    """
    r, az = pcoord(alongshore, across_shore)
    az = az - rot
    eUTM, nUTM = xycoord(r, az)
    eUTM = eUTM + eoff
    nUTM = nUTM + noff
    return eUTM, nUTM


def LatLon2UTM(lat,lon,init_epsg='epsg:26918'):
    """
    Convert lat lon (WGS84) to UTM.
    Defaults to Zone 18N

    TODO: Update to Proj 6 and correct this syntax
    """
    inProj = Proj(init='epsg:4326')
    outProj = Proj(init=init_epsg)
    outx,outy = transform(inProj,outProj,lon,lat)
    return outx, outy


def UTM2LatLon(easting, northing, initepsg='epsg:26918'):
    """
    Convert UTM to lat, lon (WGS84)
    Defaults to Zone 18N

    TODO: Update to Proj 6 and correct this syntax
    """
    outProj = Proj(init='epsg:4326')
    inProj = Proj(init=initepsg)
    lon,lat = transform(inProj,outProj,easting,northing)
    return lon, lat


def UTM2rot(xutm,yutm,r):
    """
    Convert UTM coordinates to rotated coordinates

    Now deprecated by UTM2Island ... delete
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


def yaml2dict(yamlfile):
    """Import contents of a YAML file as a dict

    Args:
        yamlfile (str): YAML file to read

    Returns:
        dict interpreted from YAML file

    Raises:

    """
    dictname = None
    with open(yamlfile, "r") as infile:
        try:
            dictname = yaml.safe_load(infile)
        except yaml.YAMLError as exc:
            print(exc)

    return dictname


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
