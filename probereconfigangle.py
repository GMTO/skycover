# calculate how long we can hold the stars for and how many reconfigurations are needed
# initially, use the saved asterism and see if I can rotate the stars

import warnings
warnings.filterwarnings("ignore")
import sys
import os
sys.path.append(os.path.expanduser('~/PYTHON/ATP/'))
import ceo
import yaml
import pandas as pd
import numpy as np
import atp
import matplotlib.pyplot as plt
from agwsprobes import *
import copy

from astropy.units import Quantity
from astropy.time import Time
import astropy.io.fits as pyfits

plot = 1
if plot:
    resm3 = pyfits.getdata('m3fieldrotation.fits')
    resdgnf = pyfits.getdata('dgnffieldrotation.fits')

    # plot the results
    plt.figure(2,figsize=(9,5.25))
    plt.clf()
    plt.plot(np.arange(1000)/10.,np.sort(resm3),color='steelblue',label='FP')
    plt.plot(np.arange(1000)/10.,np.sort(resdgnf),color='indianred',label='DGNF')
    plt.xlabel('Field percentile')        
    plt.ylabel('Maximum field rotation (degrees)')
    plt.grid()
    plt.legend()
    plt.savefig('fieldrotation.png')
    sys.exit()


def rotate(gspos,rotation_angle):
    # rotation angle in degrees
    rad = rotation_angle*np.pi/180. # angle in radians
    rotmat = np.array([[np.cos(rad),np.sin(rad)],[-np.sin(rad),np.cos(rad)]])
    return gspos@rotmat

def range_of_motion(validator,gspos,angle0):
    # gspos is the gspos at angle of 0 degrees
    # angle0 is used to determine the initial probe locations
    
    (success,loc,idx) = agwscheck(validator,rotate(gspos,angle0))
    if success == False:
        print('Something went wrong!')
        print(fieldno,success,loc,idx)
        sys.exit()
    else:
        idx0 = idx.tolist()

    # idx0 are the indices correponding to angle0

    minangle = copy.copy(angle0)
    while success and idx.tolist() == idx0:
        minangle -= 1
        (success,loc,idx) = agwscheck(validator,rotate(gspos,minangle))        

    (success,loc,idx) = agwscheck(validator,rotate(gspos,angle0))
    maxangle = copy.copy(angle0)
    while success and idx.tolist() == idx0:
        maxangle += 1
        (success,loc,idx) = agwscheck(validator,rotate(gspos,maxangle))        
           
    return [minangle+1,maxangle-1]
    

    
    # now rotate the field and see how long the asterism is valid for
#    for theta_deg in range(1,61):
#        theta_rad = theta_deg*np.pi/180
#        rotmat = np.array([[np.cos(theta_rad),np.sin(theta_rad)],[-np.sin(theta_rad),np.cos(theta_rad)]])
#        (success,loc,idx) = agwscheck(validator,gspos@rotmat)
        
#        if not success or idx.tolist() != idx0:
#            print(fieldno,theta_deg-1)
#            results[fieldno] = theta_deg-1
#            break


    
config = "m3"
girmode = "fixed"
validator = agwsinit(config)

results = np.zeros(1000)+60
for fieldno in range(1000):
    
    gsfilename = os.path.expanduser('~/')+f'PYTHON/AGWS/Fields/{config}_asterism_{fieldno:04d}.csv'
    gsdata = pd.read_csv(gsfilename)

    probefunction = gsdata['Probe function'].values
    vismag = gsdata['Visible magnitude'].values
    xpos = gsdata['xpos'].values
    ypos = gsdata['ypos'].values
    
    acoprobe = np.where(probefunction == 'aco')[0]
    tt7probe = np.where(probefunction == 'tt7')[0]

    validprobes = np.where(vismag != 0)[0]

    nstars = len(validprobes)
    gspos = np.transpose(np.array([xpos,ypos]))
    gspos = gspos[validprobes]

    angle0 = 0
    minmaxangle = range_of_motion(validator,gspos,angle0)

    results[fieldno] = minmaxangle[1]

