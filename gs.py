import warnings
import argparse
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

def main():
    parser = argparse.ArgumentParser(
        description='-----AGWS sky coverage------')
        
    parser.add_argument('config',
                    choices=('fp', 'dgwf', 'dgnf', 'gmacs'),
                    help='optical configuration')
    parser.add_argument('aomode',
                    choices=('ns', 'glao'),
                    help='AO mode')
    parser.add_argument('-field_start', dest='field_start', default=0, type=int,
                        help='start with field no')
    parser.add_argument('-field_end', dest='field_end', default=1, type=int,
                        help='end with field no')
    parser.add_argument('-writecsv', help='create csv for the star fields',
                        action='store_true')
    parser.add_argument('-writepng', help='make png plots for the star fields',
                        action='store_true')
    args = parser.parse_args()
    
    config = args.config
    if config == 'fp':
        agws_config = 'm3'
    else:
        agws_config = config
    if args.aomode == 'ns':
        maxstarsperprobe = 8 # maximum number of stars per probe to analyze
    else:
        maxstarsperprobe = 8 #12
        
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

    validator = agwsinit(agws_config)

    nfail = 0
    if args.writecsv:
        for fieldno in range(args.field_start, args.field_end+1):
            outputfilename = f'output/{config}_{(args.aomode).upper()}_asterism_{fieldno:04d}.csv'
            
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

            #these are the input catalog files for the 1000 standard fields
            fielddir = "~/PYTHON/AGWS/Fields"
            starfield = os.path.expanduser(fielddir+f"/field_{fieldno:04d}.csv")

            # make use of ATP code by Rod for managing star fields
            stars  = atp.StarField(obs,target,field=starfield,**cfg['Star Catalog'])

            vismag = 0.46*stars.I+0.54*stars.R
            # remove stars that have an NaN for magnitude
            valid = ~np.isnan(vismag)
            validpos = np.where(valid)[0]

            xposarcsec = stars.local[0,:]*180./np.pi*3600
            yposarcsec = stars.local[1,:]*180./np.pi*3600

            print("------------------------fieldno = ", fieldno, " N stars = ", len(vismag), len(valid))

            if args.writepng:
                ax = plt.axes()
                ax.set_aspect('equal')
                plt.plot(xposarcsec,yposarcsec,'.',color='darkgreen')
                
                #obscuration poly
                if config == 'fp':
                    plt.plot(np.array([-334, -495, -274, 110, 145, -334]),
                         np.array([129, 339, 661, 590, 327, 129]), '-b')
                elif config == 'gmacs':
                    plt.plot(np.array([-10.11/2, 10.11/2, 10.11/2, -10.11/2, -10.11/2])*60,
                         np.array([9.11/2, 9.11/2, -9.11/2, -9.11/2, 9.11/2])*60, '-b')

                #inner circle
                if config == 'fp':
                    minradarcsec = 357.9
                    obscuration = plt.Circle([0,0],radius=minradarcsec,color='darkblue',fill=False)
                    ax.add_artist(obscuration)

                #outer circle
                if config == 'fp':
                    maxradarcsec = 600
                elif config == 'gmacs':
                    maxradarcsec = 420
                outercircle = plt.Circle([0,0],radius=maxradarcsec,color='darkblue',fill=False)
                ax.add_artist(outercircle)
                
                #plt.plot([0],[0],'*',color='black')

                plt.xlim([-maxradarcsec*1.1, maxradarcsec*1.1])
                plt.ylim([-maxradarcsec*1.1, maxradarcsec*1.1])
                plt.xlabel('X-position (arcsec)')
                plt.ylabel('Y-position (arcsec)')
                plt.title(f'field no {fieldno}')
                plt.savefig(outputfilename.replace('csv','png'))
                plt.clf()
            
            ## run agwscheck() to remove invalid stars
            gspos = np.transpose(np.array([xposarcsec,yposarcsec]))
            for k in range(len(xposarcsec)):
                if valid[k]:
                    (success,loc,idx) = agwscheck(validator,gspos[[k],:])
                    valid[k] = success
            validpos = np.where(valid)[0]
            print("after agwscheck on individual stars: ", len(validpos))

            if len(validpos) == 0:
                print('There are no guide stars!')
                continue

            vismag = vismag[validpos]
            xposarcsec = xposarcsec[validpos]
            yposarcsec = yposarcsec[validpos]

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
            
            probes_with_stars = probesreachstars.any(axis=1) #for example: [true, true, true, true]
            n_probes_with_stars = np.sum(probes_with_stars) # number of probes that can reach any star, ideally, =4

            if args.aomode == 'ns':
                ### we next determine the best TT7 star, i.e., which star will be measured with the TT7 probe
                #we do this by predicting the TT7 error for each star
                #This error has two parts, the noise_rms and aniso_rms.

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
                tterr_for_800mas = np.array([0.00204073,0.00204073,0.00248421,0.00313843,0.00404597,0.0052902,0.00649479,
                                             0.00812451,0.0108732,0.0140705,0.018404,0.0248977,0.0369819,0.0556996,0.0910204,
                                             0.13293,0.179323,0.205629,0.205629])*1000.
                tterr_mas = tterr_for_800mas*seeingArcsec/0.8

                interp_function = scipy.interpolate.interp1d(magvec,tterr_mas,bounds_error=False)
                tt7_noise_rms = interp_function(np.clip(vismag,np.min(magvec),np.max(magvec)))
                tt7_res_rms = tt7_noise_rms + tt7_aniso_rms

                #### for now, we will use the best TT7 star, 
                ## but later we could also try the other stars if we do not succeed in 
                ## finding a suitable active optics asterism
                # still, there might be multiple probes that can point to this TT7 star. We loop over that below.

                # sort from lowest to highest error
                tt7_gs_indices = np.argsort(tt7_res_rms)
                tt7_gs_idx = tt7_gs_indices[0] #index of gs used for tt7

                tt7probes = np.where(probesreachstars[:,tt7_gs_idx])[0] #which probes can reach the selected tt7 gs
                stt7 = np.array([tt7_gs_idx])[0] #index of gs used for tt7

                results = {'minangle (deg)':[],'minradius (arcsec)':[],'maxmag (=faintest)':[],
                           'tt7 index':[],'aco indices':[],'tt7probe':[],'acoprobes':[]} 

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

                                        results['minradius (arcsec)'].append(minradius)
                                        results['minangle (deg)'].append(minangle*180./np.pi)
                                        results['maxmag (=faintest)'].append(maxmag)
                                        results['tt7 index'].append(stt7)
                                        results['tt7probe'].append(tt7probe)
                                        results['acoprobes'].append(aco_idx)
                                        results['aco indices'].append(np.array([s0[k0],s1[k1],s2[k2]]))

            elif args.aomode == 'glao':
                results = {'minangle (deg)':[],'minradius (arcsec)':[],'maxmag (=faintest)':[],
                           'gs indices':[], 'star mag':[], 'xpos arcsec':[], 'ypos arcsec':[]}

                s0 = np.where(probesreachstars[0,:])[0]
                s1 = np.where(probesreachstars[1,:])[0]
                s2 = np.where(probesreachstars[2,:])[0]
                s3 = np.where(probesreachstars[3,:])[0]

                # sort by magnitude and clip at maximum number
                s0 = s0[np.argsort(vismag[s0])[0:maxstarsperprobe]]
                s1 = s1[np.argsort(vismag[s1])[0:maxstarsperprobe]]
                s2 = s2[np.argsort(vismag[s2])[0:maxstarsperprobe]]
                s3 = s3[np.argsort(vismag[s3])[0:maxstarsperprobe]]
                
                # see which combinations are allowable
                for k0 in s0:
                    for k1 in s1:
                        for k2 in s2:
                            for k3 in s3:
                                testpos = gspos[[k0,k1,k2,k3]]
                                (success,loc,idx) = agwscheck(validator,testpos)
                                if success and np.all(idx == np.array([0,1,2,3])):
                                    
                                    # evaluate the maximum magnitude
                                    maxmag = np.max(vismag[[k0, k1, k2, k3 ]])

                                    # evaluate the angle that each set of stars makes with respect to the center
                                    az0 = np.arctan2(gspos[k0,0],gspos[k0,1])
                                    az1 = np.arctan2(gspos[k1,0],gspos[k1,1])
                                    az2 = np.arctan2(gspos[k2,0],gspos[k2,1])
                                    az3 = np.arctan2(gspos[k3,0],gspos[k3,1])
                                    
                                    az = np.sort(np.array([az0,az1,az2,az3]))
                                    minangle = np.min([az[1]-az[0],az[2]-az[1], az[3]-az[2],az[0]+2*np.pi-az[3]])
                                        
                                    radius0 = np.hypot(gspos[k0,0],gspos[k0,1])
                                    radius1 = np.hypot(gspos[k1,0],gspos[k1,1])
                                    radius2 = np.hypot(gspos[k2,0],gspos[k2,1])
                                    radius3 = np.hypot(gspos[k3,0],gspos[k3,1])

                                    minradius = np.min([radius0,radius1,radius2, radius3])
                                    
                                    results['minradius (arcsec)'].append(minradius)
                                    results['minangle (deg)'].append(minangle*180./np.pi)
                                    results['maxmag (=faintest)'].append(maxmag)
                                    results['gs indices'].append(np.array([k0,k1,k2,k3]))
                                    results['star mag'].append(np.round(vismag[[k0, k1, k2, k3 ]],6))
                                    results['xpos arcsec'].append(np.round(xposarcsec[[k0, k1, k2, k3 ]],2))
                                    results['ypos arcsec'].append(np.round(yposarcsec[[k0, k1, k2, k3 ]],2))

            df = pd.DataFrame(results)
            if len(df) == 0:
                nfail += 1
                print('------------xxxxxxxxx-----------------------------------------------------', nfail)
                df['selected']=0
            else:
                df['selected']=1
                # there are various solutions to the problem

                # (1) We select the asterims whose fourth brightest star is the brightest.
                # (2) We weight all of the stars according to their noise. In this case, we want to maximimize the sum of 1 divided by the noise terms.
                # (3) We weight all of the stars equally. In this case, we minimize the sum of the noise terms

                ## approach (1)
                sortedmags = np.sort(np.array(results['star mag']),axis=1)
                i = 3
                while sum(df['selected'])>1 and i>=0:
                    idx3 = np.where((sortedmags[:,i] != np.min(sortedmags[df['selected']==1,i])))[0]
                    df['selected'][idx3]=0
                    #print(sortedmags[df['selected']==1,:])
                    i-=1

                #we may have >1 asterisms made of same set of stars
                if sum(df['selected'])>1:
                    idx = df['selected'] == 1
                    df['selected'][np.where(np.array(idx))[0][1:]]=0 #keep the 1st asterism and throw out others
            df.to_csv(outputfilename)
    
    # get stats of the field results
    stats = {'field No.':[],'N asterisms':[],'brightest faint star':[],
                       'minradius (arcsec)':[], 'star mag':[], 'xpos arcsec':[], 'ypos arcsec':[]}
    for fieldno in range(args.field_start, args.field_end+1):
        outputfilename = f'output/{config}_{(args.aomode).upper()}_asterism_{fieldno:04d}.csv'
        df = pd.read_csv(outputfilename)
        
        stats['field No.'].append(fieldno)
        stats['N asterisms'].append(len(df))
        if len(df)>0:
            idx = df['selected'] == 1
            stats['brightest faint star'].append(round(min(df['maxmag (=faintest)']),2))
            stats['minradius (arcsec)'].append(round(min(df['minradius (arcsec)']),2))
            stats['star mag'].append(np.array(df['star mag'][idx]))
            stats['xpos arcsec'].append(np.array(df['xpos arcsec'][idx]))
            stats['ypos arcsec'].append(np.array(df['ypos arcsec'][idx]))
        else:
            stats['brightest faint star'].append(0)
            stats['minradius (arcsec)'].append(0)
            stats['star mag'].append(0)
            stats['xpos arcsec'].append(0)
            stats['ypos arcsec'].append(0)
        stats_df = pd.DataFrame(stats,index=None)
        print('total number of fields = ', fieldno - args.field_start + 1)
        print('fields with valid asterisms = ', sum(stats_df['N asterisms']>0))

    stats_df.to_csv(f'output/{args.config}_gs_{args.aomode}.csv')
        
if __name__ == "__main__":
    main()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
