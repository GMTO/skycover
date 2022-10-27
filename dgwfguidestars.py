# select the guide stars given a field in DGWF mode

import warnings
warnings.filterwarnings("ignore")

import ceo
import yaml
from yaml import Loader
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
sys.path.append(os.path.expanduser('~/PYTHON/ATP/'))
from astropy.units import Quantity
from astropy.time import Time
import astropy.io.fits as pyfits
import atp
import copy
from agwsprobes import *
import scipy.interpolate

L=25.5
nPx = 201
nLenslet = 48
threshold = 0.5 # subaperture flux threshold
gmt = ceo.GMT_MX()

def modestogmtstate(gmt,modes):
    state = gmt.state

    M1Txyz = modes[[0,1,2,12,13,14,24,25,26,36,37,38,48,49,50,60,61,62,72,73,74]]
    state['M1']['Txyz'] = np.reshape(M1Txyz,(7,3))

    M2Txyz = modes[[6,7,8,18,19,20,30,31,32,42,43,44,54,55,56,66,67,68,77,78,79]]
    state['M2']['Txyz'] = np.reshape(M2Txyz,(7,3))

    M1Rxyz = np.zeros(21)
    M1Rxyz[0:20] = modes[[3,4,5,15,16,17,27,28,29,39,40,41,51,52,53,63,64,65,75,76]]
    state['M1']['Rxyz'] = np.reshape(M1Rxyz,(7,3))
    
    M2Rxyz = np.zeros(21)
    M2Rxyz[0:20] = modes[[9,10,11,21,22,23,33,34,35,45,46,47,57,58,59,69,70,71,80,81]]
    state['M2']['Rxyz'] = np.reshape(M2Rxyz,(7,3))
    gmt^=state

girmode = "tracking"
config = "dgwf"
maxstarsperprobe = 8 # maximum number of stars per probe to analyze

# this contains information about the observation
cfg = yaml.load(open(os.path.expanduser('~/PYTHON/ATP/atp.yaml')),Loader)

# load the data about the observation
data = pyfits.getdata(os.path.expanduser('~/PYTHON/AGWS/SAO/weather_pointing_data.fits'))
temperature = data[:,0] # Temperature draw in deg C
pressure = data[:,1] # Pressure draw in mbar
humidity = data[:,2] # Relative humidity in %
windspeed = data[:,3] # Wind speed in m/s
winddir = data[:,4] # Wind direction in degrees
datetime = data[:,5] # Date vector (Clay is first 500, Baade is second 500)
scira = data[:,6] # RA vector (degrees; Clay is first 500, Baade is second 500)
scidec = data[:,7] #  Declination vector (degrees; Clay is first 500, Baade is second 500)
telaz = data[:,8] # Azimuth angle of telescope (degrees; Clay is first 500, Baade is second 500)
telel = data[:,9] # Elevation angle of telescope (degrees; Clay is first 500, Baade is second 500)
r0 = data[:,10] # r0 (cm)
L0 = data[:,11] # L0 (m)

validator = agwsinit(config)

for fieldno in range(0,1000):
    telzen = (90-telel[fieldno])*np.pi/180.
    airmass = 1./np.cos(telzen)

    cfg['Observation']['time'] = Time(datetime[fieldno] , format='decimalyear', scale='utc')
    cfg['Target']['ra'] = scira[fieldno]
    cfg['Target']['dec'] = scidec[fieldno]

    cfg['Target']['pointing alt/az']['alt'][0] = telel[fieldno]
    cfg['Target']['pointing alt/az']['az'][0] = telaz[fieldno]

    cfg['Atmosphere']['r0'][0] = r0[fieldno]*100.
    cfg['Atmosphere']['L0'][0] = L0[fieldno]

    # calculate the seeing 
    gs_wavelength = Quantity(*cfg['SH']['guide star']['wavelength']).to('m').value
    r0_wavelength = Quantity(*cfg['Atmosphere']['wavelength']).to('m').value
    r0_val = Quantity(*cfg['Atmosphere']['r0']).to('m').value
    r0_val *= atp.r0_scaling(r0_wavelength,gs_wavelength,telzen)
    seeingRad = gs_wavelength/r0_val
    seeingArcsec = seeingRad*ceo.constants.RAD2ARCSEC
    print("seeing (arcsec) : %5.3f" %seeingArcsec)
    print("elevation (degrees) : %4.1f" %telel[fieldno])
    
    obs = atp.Observatory(**cfg['Observatory'],**cfg['Observation'])
    target = atp.Target(obs,**cfg['Target'])

    fielddir = "~/PYTHON/AGWS/Fields"
    starfield = os.path.expanduser(fielddir+f"/field_{fieldno:04d}.csv")

    outputfilename = os.path.expanduser('~/')+f'PYTHON/AGWS/Fields/{config}_asterism_{fieldno:04d}.csv'
    probefunction = [None,None,None,None]
    gsno = [None,None,None,None] # guide star number
    vismaggs = [0,0,0,0]
    xposgs = [0,0,0,0]
    yposgs = [0,0,0,0]
   
    stars  = atp.StarField(obs,target,field=starfield,**cfg['Star Catalog'])

    vismag = 0.46*stars.I+0.54*stars.R
    # remove stars that have an NaN for magnitude
    valid = ~np.isnan(vismag)
    validpos = np.where(valid)[0]
    
    xposarcsec = stars.local[0,:]*180./np.pi*3600
    yposarcsec = stars.local[1,:]*180./np.pi*3600

    plot = 0
    if plot:
        plt.figure(8)
        plt.clf()
        ax = plt.axes()
        ax.set_aspect('equal')
        plt.xlim((-610,610))
        plt.ylim((-610,610))
        plt.plot(xposarcsec[validpos],yposarcsec[validpos],'.',color='red')
        
    gspos = np.transpose(np.array([xposarcsec,yposarcsec]))

    nstars = len(xposarcsec)
    for k in range(nstars):
        if valid[k]:
            (success,loc,idx) = agwscheck(validator,gspos[[k],:])
            valid[k] = success

    validpos = np.where(valid)[0]
    
    if len(validpos) == 0:
        print('There are no guide stars!')
        bestasterism = {'Probe number': [0,1,2,3], 'Probe function': probefunction,'Guide star':gsno,'Visible magnitude':vismaggs,'xpos':xposgs,'ypos':yposgs}
        
        print(outputfilename)
        df = pd.DataFrame(bestasterism,index=None)
        print (df)
        df.to_csv(outputfilename, index = False, header=True)
        continue
    
    vismag = vismag[validpos]
    xposarcsec = xposarcsec[validpos]
    yposarcsec = yposarcsec[validpos]

    if plot:
        if config == "dgwf":
            minradarcsec = 365 # trial and error
            maxradarcsec = 572.4 # trial and error
        if config == "dgnf":
            # the AGWS must be > 6' away to avoid vignetting the field
            minradarcsec = 364.7 
            maxradarcsec = 600
        if config == "m3":
            # radius of guide star must be greater than 357.9 mm
            minradarcsec = 365.07
            maxradarcsec = 600

        plt.plot(xposarcsec,yposarcsec,'.',color='darkgreen')
        #innercircle = plt.Circle([0,0],radius=minradarcsec,color='darkblue',fill=False)
        #ax.add_artist(innercircle)
        #outercircle = plt.Circle([0,0],radius=maxradarcsec,color='darkblue',fill=False)
        #ax.add_artist(outercircle)
        plt.plot([0],[0],'*',color='black')

        plt.xlabel('X-position (arcsec)')
        plt.ylabel('Y-position (arcsec)')

    gspos = np.transpose(np.array([xposarcsec,yposarcsec]))
    nvalidstars = len(xposarcsec)    

    # calculate the TT7 error
    # find the TT7 error for all stars
    print('Calculating TT7 error')
    dist = np.arange(1,102)/10.
    tt7_aniso_dist = [atp.tt7_tt_error(zz,0.,telzen,**cfg) for zz,magnitude in zip(dist,dist*0.)]
    aniso_interp_function = scipy.interpolate.interp1d(dist,tt7_aniso_dist,bounds_error=False)

    radialdistArcmin = np.hypot(xposarcsec,yposarcsec)/60.
    tt7_aniso_rms = aniso_interp_function(radialdistArcmin)
    
    # saturation is attained for magnitude 10 stars
    # values calculated using tt7noise.i for 0.8" seeing
    magvec = np.array([0,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,20])
    # this is the error in the segment tip-tilt estimate using a 24x24 SH WFS
    tterr_for_800mas = np.array([0.00204073,0.00204073,0.00248421,0.00313843,0.00404597,0.0052902,0.00649479,0.00812451,0.0108732,0.0140705,0.018404,0.0248977,0.0369819,0.0556996,0.0910204,0.13293,0.179323,0.205629,0.205629])*1000.
    tterr_mas = tterr_for_800mas*seeingArcsec/0.8
    
    interp_function = scipy.interpolate.interp1d(magvec,tterr_mas,bounds_error=False)
    tt7_noise_rms = interp_function(np.clip(vismag,np.min(magvec),np.max(magvec)))
    tt7_res_rms = tt7_noise_rms + tt7_aniso_rms
    
    # determine which probes reach which stars
    print('Finding which probes reach which stars')
    probesreachstars = agwsreachstars(validator,gspos)

    # calculate the obscuration by each probe for each star
    startloc = np.array([[0,0.159],[-0.159,0],[0,-0.159],[0.159,0]])*3600.
    (success,loc,idx) = agwscheck(validator,startloc)
    minshadow = validator.shadowfrac
    
    shadowfraction = np.zeros(probesreachstars.shape)
    for k1 in range(4):
        for k2 in range(nvalidstars):
            if probesreachstars[k1,k2]:
                pos = startloc.copy()
                pos[k1,:] = gspos[k2,:]
                (success,loc,idx) = agwscheck(validator,pos)
                shadowfraction[k1,k2] = validator.shadowfrac - 0.75*minshadow

    probes_with_stars = probesreachstars.any(axis=1) 
    n_probes_with_stars = np.sum(probes_with_stars)

    if n_probes_with_stars == 4:
        print('Finding four suitable stars')
        
        # find the valid 4 star asterisms that meet the shadowing requirement
        s0 = np.where(probesreachstars[0,:])[0]
        s1 = np.where(probesreachstars[1,:])[0]
        s2 = np.where(probesreachstars[2,:])[0]
        s3 = np.where(probesreachstars[3,:])[0]

        # limit the calculation to the brightest 30 stars for each probe
        s0 = s0[np.argsort(vismag[s0])][0:30]
        s1 = s1[np.argsort(vismag[s1])][0:30]
        s2 = s2[np.argsort(vismag[s2])][0:30]
        s3 = s3[np.argsort(vismag[s3])][0:30]            
        
        pos = np.zeros((4,2))
        validasterisms = []

        # keep track of the shadowing and report it        
        for m0 in s0.tolist():
            frac0 = shadowfraction[0,m0]
            pos[0,:] = gspos[m0,:]

            for m1 in s1.tolist():
                if m1 == m0:
                    continue
                pos[1,:] = gspos[m1,:]

                frac1 = shadowfraction[1,m1]
                for m2 in s2.tolist():
                    if m2 == m1 or m2 == m0:
                        continue
                    frac2 = shadowfraction[2,m2]
                    pos[2,:] = gspos[m2,:]

                    for m3 in s3.tolist():
                        if m3 == m2 or m3 == m1 or m3 == m0:
                            continue
                        frac3 = shadowfraction[3,m3]
                        pos[3,:] = gspos[m3,:]

                        frac = frac0+frac1+frac2+frac3 
                        if frac <= 0.15:
                            # check if asterism is valid
                            (success,loc,idx) = agwscheck(validator,pos)
                            if success:
                                validasterisms.append([m0,m1,m2,m3])
                            
        nAsterism = len(validasterisms)
        if nAsterism > 0:

            # select the 25 asterisms with the brightest faint star
            maxasterisms = 25
            maxmag = np.zeros(nAsterism)
            counter = 0
            for k,asterism in enumerate(validasterisms):
                maxmag[k] = np.max(vismag[asterism])

            sortmag = np.argsort(maxmag)
            brightasterisms = np.array(validasterisms)[sortmag][0:np.clip(nAsterism,0,maxasterisms)]

            # find the best TT7 star in those bright asterisms
            brightstars = np.sort(np.array(list(set(np.ndarray.flatten(brightasterisms)))))
                        
            # sort from lowest to highest error
            tt7_gs_indices = np.argsort(tt7_res_rms[brightstars])
            tt7idx = brightstars[tt7_gs_indices[0]]

            # find the asterisms that include the brightest TT7 star
            selectedasterisms = [] 
            for asterism in brightasterisms:
                if tt7idx in asterism:
                    selectedasterisms.append(asterism)

            # find the best performing active optics asterism
            nAsterism = len(selectedasterisms)
            res = np.zeros(nAsterism)
            for ast,asterism in enumerate(selectedasterisms):
                gsidx = asterism.tolist()
                gsidx.remove(tt7idx)
                gsidx = np.array(gsidx)

                # generate a reconstructor for each one for one segment
                wfs = ceo.GeometricShackHartmann(nLenslet,L/nLenslet,3,coupled=True)

                pos = gspos[gsidx,:]
                zen = [np.hypot(pos[i,0],pos[i,1])/206265. for i in range(3)]
                azi = [np.arctan2(pos[i,0],pos[i,1]) for i in range(3)]
                mag = vismag[gsidx]

                gs = ceo.Source(photometric_band="R+I",zenith=zen,azimuth=azi,magnitude=np.array([10,10,10]),rays_box_size=L,rays_box_sampling=nLenslet*8+1,rays_origin=[0,0,25])

                # run my own calibration, because the interaction matrix is not consistent with the measured values
                gs.reset()
                gmt.reset()
                gmt.propagate(gs)
                wfs.calibrate(gs,threshold)
                gs>>(gmt,wfs)
                nmes = len(wfs.get_measurement())
                nvl = np.sum(wfs.valid_lenslet.f.host(shape=(3,48**2))>0,axis=1)

                magvec = np.array([0,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20])

                # this is the error in the segment tip-tilt estimate using a 24x24 SH WFS
                meserr_for_800mas = np.array([0.00136217,0.00136217,0.00171104,0.00216745,0.00275779,0.00351932,0.00453948,0.00590709,0.00780148,0.0104783,0.0144437,0.0204345,0.0296933,0.044133,0.0674427,0.110046]) # measurement error per subaperture
                meserr_mas = meserr_for_800mas*seeingArcsec/0.8
                sh_interp_function = scipy.interpolate.interp1d(magvec,meserr_mas,bounds_error=False)

                noisediag = np.zeros(nmes)
                idx0 = 0
                n2 = int(nmes/2)
                for c1 in range(3):
                    sh_noise_rms = sh_interp_function(np.clip(mag[c1],np.min(magvec),np.max(magvec)))
                    noisediag[idx0:idx0+nvl[c1]] = sh_noise_rms # noise in arcsec
                    noisediag[idx0+n2:idx0+n2+nvl[c1]] = sh_noise_rms # noise in arcsec
                    idx0 += nvl[c1]

                poke = 5e-6
                Rxyz = [3,4,5,9,10,11]# modes are M1S1Rxyz and M2S1Rxyz
                nmodes = 82
                H = np.zeros((nmes,nmodes))

                ~gmt
                for k in Rxyz:
                    modes = np.zeros(nmodes)
                    modes[k] = poke
                    modestogmtstate(gmt,modes)
                    ~wfs
                    gs>>(gmt,wfs)
                    +gs
                    +wfs
                    mes = wfs.get_measurement()
                    H[:,k] = mes/poke

                H = H[:,Rxyz]
                #penmat = np.eye(6)*1e-7*0.
                #Hinv = np.linalg.solve(np.transpose(H)@H+penmat,np.transpose(H))
                Hinv = np.linalg.solve(np.transpose(H)@H,np.transpose(H))
                Cnn = np.diag(noisediag)

                J = np.trace(Hinv@Cnn@np.transpose(Hinv))
                res[ast] = J

            sidx = np.argsort(res)[0]

            asterism = selectedasterisms[sidx]
            tt7probe = np.where(tt7idx == asterism)[0][0]       
            tt7star = validpos[tt7idx]

            acoprobes = list(range(4))
            acoprobes.remove(tt7probe)
            acostars = validpos[gsidx]

            probefunction = ['aco','aco','aco','aco']
            probefunction[tt7probe] = 'tt7'

            gsno = validpos[asterism]
            vismaggs = vismag[asterism]
            xposgs = xposarcsec[asterism]
            yposgs =  yposarcsec[asterism]

            bestasterism = {'Probe number': [0,1,2,3], 'Probe function': probefunction,'Guide star':gsno,'Visible magnitude':vismaggs,'xpos':xposgs,'ypos':yposgs}

            # report the best asterism
            print(outputfilename)
            df = pd.DataFrame(bestasterism,index=None)
            print (df)
            df.to_csv(outputfilename, index = False, header=True)
            continue
    
    # two or more probes with stars, but no valid 4 star asterism
    print('No valid 4 star asterism')
    
    # is there a valid 3 star asterism?
    nguidestars = gspos.shape[0]
    
    if n_probes_with_stars >= 3:



        
        continue

    print('No valid 3 star asterism')
    
    # is there a valid 2 star asterism?
    if n_probes_with_stars >= 2:
        # check all the two star combinations

        wp = np.where(probes_with_stars)[0]
        
        s0 = np.where(probesreachstars[wp[0],:])[0]
        s1 = np.where(probesreachstars[wp[1],:])[0]
        
        validasterisms = []
        for m0 in s0.tolist():
            for m1 in s1.tolist():
                pos = gspos[[m0,m1],:]

                # check if asterism is valid
                (success,loc,idx) = agwscheck(validator,pos)
                if success:
                    validasterisms.append([m0,m1])

        nAsterism = len(validasterisms)

        if nAsterism >= 1:
            # select the asterism that has the brightest faintest star
            starmags = vismag[np.array(validasterisms)]
            maxmag = np.max(starmags,axis=1)
            sidx = np.argsort(maxmag)[0]            

            asterism = validasterisms[sidx]
            mags = starmags[sidx]
            
            # brightest star in asterism is used for TT7
            tt7idx = np.where(mags == np.min(mags))[0][0]
            gsidx = list(asterism)
            gsidx.remove(tt7idx)
            gsidx = np.array(gsidx)[0]

            pos = gspos[[tt7idx,gsidx],:]
            (success,loc,idx) = agwscheck(validator,pos)

            tt7probe = np.where(idx == 0)[0][0]
            tt7star = validpos[tt7idx]

            acoprobes = np.where(idx == 1)[0][0]
            acostars = validpos[gsidx]

            probefunction[tt7probe] = 'tt7'
            probefunction[acoprobes] = 'aco'

            gsno = [None,None,None,None] # guide star number
            gsno[tt7probe] = tt7star
            gsno[acoprobes] = acostars
            vismaggs[tt7probe] = vismag[tt7idx]
            vismaggs[acoprobes] = vismag[gsidx]
            xposgs[tt7probe] = xposarcsec[tt7idx]
            yposgs[tt7probe] = yposarcsec[tt7idx]
            xposgs[acoprobes] = xposarcsec[gsidx]
            yposgs[acoprobes] = yposarcsec[gsidx]
                            
            bestasterism = {'Probe number': [0,1,2,3], 'Probe function': probefunction,'Guide star':np.array(gsno),'Visible magnitude':vismaggs,'xpos':xposgs,'ypos':yposgs}

            print(outputfilename)
            df = pd.DataFrame(bestasterism,index=None)
            print (df)
            df.to_csv(outputfilename, index = False, header=True)
            continue
                
    print('No valid 2 star asterism')    

    
