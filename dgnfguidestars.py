# select the guide stars given a field

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
config = "m3"
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
    for k in range(len(xposarcsec)):
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
        innercircle = plt.Circle([0,0],radius=minradarcsec,color='darkblue',fill=False)
        ax.add_artist(innercircle)
        outercircle = plt.Circle([0,0],radius=maxradarcsec,color='darkblue',fill=False)
        ax.add_artist(outercircle)
        plt.plot([0],[0],'*',color='black')

        plt.xlabel('X-position (arcsec)')
        plt.ylabel('Y-position (arcsec)')

    gspos = np.transpose(np.array([xposarcsec,yposarcsec]))

    # determine which probes reach which stars
    probesreachstars = agwsreachstars(validator,gspos)

    # select the segment tip-tilt star and the associated probe
    # zenith angle in radians (used to find the anisoplanatic error)

    # find the star with the lowest segment tip-tilt error
    # inputs: distance from center in arcminutes, mV (currently); need to fix
    radialdistArcmin = np.hypot(xposarcsec,yposarcsec)/60.

    # calculate the seeing 
    gs_wavelength = Quantity(*cfg['SH']['guide star']['wavelength']).to('m').value
    r0_wavelength = Quantity(*cfg['Atmosphere']['wavelength']).to('m').value
    r0_val = Quantity(*cfg['Atmosphere']['r0']).to('m').value
    r0_val *= atp.r0_scaling(r0_wavelength,gs_wavelength,telzen)
    seeingRad = gs_wavelength/r0_val
    seeingArcsec = seeingRad*ceo.constants.RAD2ARCSEC
    print("seeing (arcsec) : %5.3f" %seeingArcsec)

    # TT7 tip-tilt error per subaperture (24x24 subapertures)
    # use the interaction matrix to obtain the error in segment tip-tilt
    # 

    # calculate at what magnitude we saturate and define as minimum magnitude

    # calculate the anisoplanatic error by setting magnitude to 0
    tt7_aniso_rms = [atp.tt7_tt_error(zz,0.,telzen,**cfg) for zz,magnitude in zip(radialdistArcmin,vismag)]
    # calculate the approximate TT7 error
    # saturation is attained for magnitude 10 stars
    # values calculated using tt7noise.i for 0.8" seeing
    magvec = np.array([0,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,20])
    # this is the error in the segment tip-tilt estimate using a 24x24 SH WFS
    tterr_for_800mas = np.array([0.00204073,0.00204073,0.00248421,0.00313843,0.00404597,0.0052902,0.00649479,0.00812451,0.0108732,0.0140705,0.018404,0.0248977,0.0369819,0.0556996,0.0910204,0.13293,0.179323,0.205629,0.205629])*1000.
    tterr_mas = tterr_for_800mas*seeingArcsec/0.8

    interp_function = scipy.interpolate.interp1d(magvec,tterr_mas,bounds_error=False)
    tt7_noise_rms = interp_function(np.clip(vismag,np.min(magvec),np.max(magvec)))
    tt7_res_rms = tt7_noise_rms + tt7_aniso_rms

    # sort from lowest to highest error
    tt7_gs_indices = np.argsort(tt7_res_rms)
    tt7_gs_idx = tt7_gs_indices[0]

    # for now, we will use the best TT7 star, but later we could also try the other stars if we do not succeed in finding a suitable active optics asterism

    if fieldno == 810:
        if config == 'dgnf' or config == 'm3':
            # here, we have to use the second best TT7 star
            tt7_gs_idx = tt7_gs_indices[1]
            
    tt7probes = np.where(probesreachstars[:,tt7_gs_idx])[0]
    stt7 = np.array([tt7_gs_idx])[0]
    
    probes_with_stars = probesreachstars.any(axis=1) 
    n_probes_with_stars = np.sum(probes_with_stars)

    results = {'minangle':[],'minradius':[],'maxmag':[],'tt7':[],'aco':[],'tt7probe':[],'acoprobes':[]} 
        
    if n_probes_with_stars == 4:  
        for tt7probe in tt7probes:    
            aco_idx = np.arange(4)
            aco_idx = np.delete(aco_idx,tt7probe)

            s2 = np.where(probesreachstars[aco_idx[2],:])[0]
            if s2.tolist() == []:
                np.delete(aco_idx,aco_idx[2])

            s1 = np.where(probesreachstars[aco_idx[1],:])[0]
            if s1.tolist() == []:
                np.delete(aco_idx,aco_idx[1])

            s0 = np.where(probesreachstars[aco_idx[0],:])[0]
            if s0.tolist() == []:
                np.delete(aco_idx,aco_idx[0])

            # sort according to the visible magnitude
            s0 = s0[np.argsort(vismag[s0])]
            s1 = s1[np.argsort(vismag[s1])]
            s2 = s2[np.argsort(vismag[s2])]

            m0 = np.clip(len(s0),0,maxstarsperprobe)
            m1 = np.clip(len(s1),0,maxstarsperprobe)
            m2 = np.clip(len(s2),0,maxstarsperprobe)

            testpos = np.zeros((4,2))
            testpos[tt7probe,:] = gspos[stt7]

            for k0 in range(m0):
                testpos[aco_idx[0],:] = gspos[s0[k0],:]
                for k1 in range(m1):
                    testpos[aco_idx[1],:] = gspos[s1[k1],:]
                    for k2 in range(m2):
                        testpos[aco_idx[2],:] = gspos[s2[k2],:]

                        (success,loc,idx) = agwscheck(validator,testpos)
                        if success:
                            # evaluate the maximum magnitude
                            maxmag = np.max(vismag[[s0[k0],s1[k1],s2[k2]]])

                            # evaluate the angle that each set of stars makes with respect to the center
                            az0 = np.arctan2(gspos[s0[k0],0],gspos[s0[k0],1])
                            az1 = np.arctan2(gspos[s1[k1],0],gspos[s1[k1],1])
                            az2 = np.arctan2(gspos[s2[k2],0],gspos[s2[k2],1])

                            radius0 = np.hypot(gspos[s0[k0],0],gspos[s0[k0],1])
                            radius1 = np.hypot(gspos[s1[k1],0],gspos[s1[k1],1])
                            radius2 = np.hypot(gspos[s2[k2],0],gspos[s2[k2],1])

                            az = np.sort(np.array([az0,az1,az2]))
                            minangle = np.min([az[1]-az[0],az[2]-az[1],az[0]+2*np.pi-az[2]])

                            minradius = np.min([radius0,radius1,radius2])

                            results['minradius'].append(minradius)
                            results['minangle'].append(minangle*180./np.pi)
                            results['maxmag'].append(maxmag)
                            results['tt7'].append(stt7)
                            results['tt7probe'].append(tt7probe)
                            results['acoprobes'].append(aco_idx)
                            results['aco'].append(np.array([s0[k0],s1[k1],s2[k2]]))
        nAsterism = len(results['maxmag'])
    else:
        nAsterism = 0

    if nAsterism > 0:
        validAst = np.ones(nAsterism,dtype=bool)

        for k1 in range(nAsterism):
            for k2 in range(k1+1,nAsterism):
                if results['maxmag'][k1] > results['maxmag'][k2] and results['minangle'][k1] < results['minangle'][k2] and results['minradius'][k1] < results['minradius'][k2]:
                    validAst[k1] = 0
                if results['maxmag'][k2] > results['maxmag'][k1] and results['minangle'][k2] < results['minangle'][k1] and results['minradius'][k2] < results['minradius'][k1]:
                    validAst[k2] = 0

        wv = np.where(validAst)[0]
        nAsterism = len(wv)

        for key in results.keys():
            results[key] = [results[key][i] for i in wv]
                
        # generate a reconstructor for each one for one segment
        wfs = ceo.GeometricShackHartmann(nLenslet,L/nLenslet,3,coupled=True)
        res = np.zeros(nAsterism)
        for ast in range(nAsterism):
            gsidx = results['aco'][ast]
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

        tt7star = validpos[results['tt7'][sidx]].tolist()
        acostars = validpos[results['aco'][sidx]].tolist()

        tt7probe = results['tt7probe'][sidx]
        acoprobes = results['acoprobes'][sidx]

        probefunction[tt7probe] = 'tt7'
        gsno[tt7probe] = int(tt7star)
        vismaggs[tt7probe] = vismag[results['tt7'][sidx]]
        xposgs[tt7probe] = xposarcsec[results['tt7'][sidx]]
        yposgs[tt7probe] = yposarcsec[results['tt7'][sidx]]

        for idx,k in enumerate(acoprobes.tolist()):
            probefunction[k] = 'aco'
            gsno[k] = int(acostars[idx])
            vismaggs[k] = vismag[results['aco'][sidx][idx]]
            xposgs[k] = xposarcsec[results['aco'][sidx][idx]]
            yposgs[k] = yposarcsec[results['aco'][sidx][idx]]

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
        testpos = np.zeros((3,2))
        testpos[0,:] = gspos[stt7]
        
        # loop over all of the other combinations
        for k0 in range(nguidestars):
            if k0 == stt7:
                continue

            testpos[1,:] = gspos[k0]

            for k1 in range(k0+1,nguidestars):
                if k1 == stt7:
                    continue

                testpos[2,:] = gspos[k1]

                (success,loc,idx) = agwscheck(validator,testpos)

                if success:
                    tt7probe = np.where(idx == 0)[0][0]
                    acoprobes = np.array([np.where(idx == 1)[0][0],np.where(idx == 2)[0][0]])
                    results['tt7'].append(stt7)
                    results['tt7probe'].append(tt7probe)
                    results['acoprobes'].append(acoprobes)
                    results['aco'].append([k0,k1])
                                              
        nAsterism = len(results['tt7'])

        if nAsterism == 0:
            continue

        # compare the possible asterisms
        # generate a reconstructor for each one for one segment
        wfs = ceo.GeometricShackHartmann(nLenslet,L/nLenslet,2,coupled=True)
        res = np.zeros(nAsterism)
        for ast in range(nAsterism):
            gsidx = results['aco'][ast]
            pos = gspos[gsidx,:]
            zen = [np.hypot(pos[i,0],pos[i,1])/206265. for i in range(2)]
            azi = [np.arctan2(pos[i,0],pos[i,1]) for i in range(2)]
            mag = vismag[gsidx]

            gs = ceo.Source(photometric_band="R+I",zenith=zen,azimuth=azi,magnitude=np.array([10,10,10]),rays_box_size=L,rays_box_sampling=nLenslet*8+1,rays_origin=[0,0,25])

            # run my own calibration, because the interaction matrix is not consistent with the measured values
            gs.reset()
            gmt.reset()
            gmt.propagate(gs)
            wfs.calibrate(gs,threshold)
            gs>>(gmt,wfs)
            nmes = len(wfs.get_measurement())
            nvl = np.sum(wfs.valid_lenslet.f.host(shape=(2,48**2))>0,axis=1)

            magvec = np.array([0,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20])

            # this is the error in the segment tip-tilt estimate using a 24x24 SH WFS
            meserr_for_800mas = np.array([0.00136217,0.00136217,0.00171104,0.00216745,0.00275779,0.00351932,0.00453948,0.00590709,0.00780148,0.0104783,0.0144437,0.0204345,0.0296933,0.044133,0.0674427,0.110046]) # measurement error per subaperture
            meserr_mas = meserr_for_800mas*seeingArcsec/0.8
            sh_interp_function = scipy.interpolate.interp1d(magvec,meserr_mas,bounds_error=False)

            noisediag = np.zeros(nmes)
            idx0 = 0
            n2 = int(nmes/2)
            for c1 in range(2):
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

        tt7star = validpos[results['tt7'][sidx]].tolist()
        acostars = validpos[results['aco'][sidx]].tolist()

        tt7probe = results['tt7probe'][sidx]
        acoprobes = results['acoprobes'][sidx]

        probefunction[tt7probe] = 'tt7'
        gsno[tt7probe] = int(tt7star)
        vismaggs[tt7probe] = vismag[results['tt7'][sidx]]
        xposgs[tt7probe] = xposarcsec[results['tt7'][sidx]]
        yposgs[tt7probe] = yposarcsec[results['tt7'][sidx]]

        for idx,k in enumerate(acoprobes.tolist()):
            probefunction[k] = 'aco'
            gsno[k] = int(acostars[idx])
            vismaggs[k] = vismag[results['aco'][sidx][idx]]
            xposgs[k] = xposarcsec[results['aco'][sidx][idx]]
            yposgs[k] = yposarcsec[results['aco'][sidx][idx]]

        bestasterism = {'Probe number': [0,1,2,3], 'Probe function': probefunction,'Guide star':np.array(gsno),'Visible magnitude':vismaggs,'xpos':xposgs,'ypos':yposgs}

        print(outputfilename)
        df = pd.DataFrame(bestasterism,index=None)
        print (df)
        df.to_csv(outputfilename, index = False, header=True)
        continue

    print('No valid 3 star asterism')
    
    # is there a valid 2 star asterism?
    nguidestars = gspos.shape[0]
    if n_probes_with_stars >= 2:
        testpos = np.zeros((2,2))
        testpos[0,:] = gspos[stt7]
        
        # loop over all of the other combinations
        for k0 in range(nguidestars):
            if k0 == stt7:
                continue

            testpos[1,:] = gspos[k0]
            (success,loc,idx) = agwscheck(validator,testpos)

            if success:
                tt7probe = np.where(idx == 0)[0][0]
                acoprobes = np.array(np.where(idx == 1)[0])
                results['tt7'].append(stt7)
                results['tt7probe'].append(tt7probe)
                results['acoprobes'].append(acoprobes)
                results['aco'].append(k0)
                                              
        nAsterism = len(results['tt7'])
        
        # make the brightest star the active optics star
        mag = vismag[results['aco']]
        sidx = np.argsort(mag)[0]

        tt7star = validpos[results['tt7'][sidx]].tolist()
        acostars = validpos[results['aco'][sidx]].tolist()

        tt7probe = results['tt7probe'][sidx]
        acoprobes = results['acoprobes'][sidx][0]

        probefunction[tt7probe] = 'tt7'
        probefunction[acoprobes] = 'aco'
        
        gsno[tt7probe] = int(tt7star)
        vismaggs[tt7probe] = vismag[results['tt7'][sidx]]
        xposgs[tt7probe] = xposarcsec[results['tt7'][sidx]]
        yposgs[tt7probe] = yposarcsec[results['tt7'][sidx]]

        gsno[acoprobes] = int(acostars)
        vismaggs[acoprobes] = vismag[results['aco'][sidx]]
        xposgs[acoprobes] = xposarcsec[results['aco'][sidx]]
        yposgs[acoprobes] = yposarcsec[results['aco'][sidx]]
                
        bestasterism = {'Probe number': [0,1,2,3], 'Probe function': probefunction,'Guide star':np.array(gsno),'Visible magnitude':vismaggs,'xpos':xposgs,'ypos':yposgs}

        print(outputfilename)
        df = pd.DataFrame(bestasterism,index=None)
        print (df)
        df.to_csv(outputfilename, index = False, header=True)
        continue
                
    print('No valid 2 star asterism')    

    # make the TT7 star the active optics probe
    tt7probes = np.where(probesreachstars[:,tt7_gs_idx])[0]
    tt7probe = tt7probes[0] # it does not matter which probe we use 
    stt7 = np.array([tt7_gs_idx])[0]
    tt7star = validpos[stt7].tolist()
    
    probefunction[tt7probe] = 'aco'
    gsno[tt7probe] = int(tt7star)
    vismaggs[tt7probe] = vismag[stt7]
    xposgs[tt7probe] = xposarcsec[stt7]
    yposgs[tt7probe] = yposarcsec[stt7]
        
    bestasterism = {'Probe number': [0,1,2,3], 'Probe function': probefunction,'Guide star':np.array(gsno),'Visible magnitude':vismaggs,'xpos':xposgs,'ypos':yposgs}
        
    print(outputfilename)
    df = pd.DataFrame(bestasterism,index=None)
    print (df)
    df.to_csv(outputfilename, index = False, header=True)
