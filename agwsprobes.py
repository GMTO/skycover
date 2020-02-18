## assign a probe to each guide star if possible and return "True"
## otherwise, return "False"

import agwsvalid
import numpy as np
import sys
import itertools

def agwsprobes(obsmode,gsloc):
    # gsloc is currently a 4x2 matrix of guide star positions
    # will be upgraded to accept fewer guide stars
    # obsmode is the observing mode, 'dgnf', 'm3', 'dgwf' 

    obsmode = obsmode.lower()
    
    if obsmode not in ['dgnf','m3','dgwf']:
        print('*** obsmode must be one of dgnf, m3 or dgwf ***')
        return (False,gsloc)

    dims = gsloc.shape
    if len(dims) != 2:
        print('*** gsloc must be an Nx2 array of star coordinates ***')
        return (False,gsloc)

    if dims[1] != 2:
        print('*** gsloc must be an Nx2 array of star coordinates ***')
        return (False,gsloc)

    nstars = dims[0]
    
    if nstars > 4:
        print('*** gsloc must have 4 or fewer stars ***')
        return (False,gsloc)
    
    validator = agwsvalid.validator(obsmode,silent=True)

    if obsmode == 'dgnf':
        startloc = np.array([[0,0.166],[-0.166,0],[0,-0.166],[0.166,0]])
    if obsmode == 'm3':
        startloc = 0.166*np.array([[np.sin(0.2),np.cos(0.2)],[-1,0.],[0,-1],[1,0]])
    if obsmode == 'dgwf':
        startloc = np.array([[0,0.159],[-0.159,0],[0,-0.159],[0.159,0]])
            
    if nstars < 4:
        gslocs = np.concatenate((gsloc,startloc),axis=0) 
        combinations = np.array(list(itertools.combinations([0,1,2,3,4],4)))

        if nstars == 1:
            ncomb = 4
        if nstars == 2:
            ncomb = 6
        if nstars == 3:
            ncomb = 4

        combinations = combinations[0:ncomb,:] # only include the ones with the specified gsloc
        permutations = np.array(list(itertools.permutations([0,1,2,3])))
        for comb in combinations:
            gsloc = gslocs[comb,:]
            for perm in permutations:
                loc = gsloc[perm,:]
                t = np.ndarray.flatten(loc)
                check = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
                if check:
                    gsloc = loc
                    return (True,gsloc)

        return (False,gsloc)
    
    permutations = np.array(list(itertools.permutations([0,1,2,3])))
    for perm in permutations:
        loc = gsloc[perm,:]
        t = np.ndarray.flatten(loc)
        check = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
        if check:
            gsloc = loc
            return (True,gsloc)

    return (False,gsloc)

## let us try an example with four guide stars and let the software decide on the order

dist = 0.159

gsloc = np.array([[0,dist],[-dist,0],[0,-dist],[dist,0]])
gsloc = gsloc[[3,2,0,1],:]

(success,loc) = agwsprobes('dgnf',gsloc)
print(success)

gs2loc = gsloc[0:2,:]
(success,loc) = agwsprobes('m3',gs2loc)
print(success)

gs1loc = np.array([[0.1,0.1]])
(success,loc) = agwsprobes('dgwf',gs1loc)
print(success)

