## assign a probe to each guide star if possible using the agwsvalid tools
## Marcos van Dam, Flat Wavefronts, February 202

import agwsvalid
import numpy as np
import sys
import itertools
import math

def nCk(n,k):
    return np.int(math.factorial(n)/(math.factorial(n-k)*math.factorial(k)))

def agwsinit(obsmode,silent=True):

    obsmode = obsmode.lower()
    
    if obsmode not in ['dgnf','m3','dgwf','gmacs']:
        print('*** obsmode must be one of dgnf, m3, dgwf or gmacs ***')
        return (False,gsloc)

    validator = agwsvalid.validator(obsmode,silent=silent)
    validator.obsmode = obsmode
    
    return validator

def agwsreachstars(validator,gspos):

    nstars = (gspos.shape)[0]
    results = np.zeros((4,nstars),dtype=bool)

    if validator.obsmode == 'dgnf':
        startloc = np.array([[0,0.166],[-0.166,0],[0,-0.166],[0.166,0]])*3600.
    elif validator.obsmode == 'm3':
        startloc = 0.166*np.array([[np.sin(0.2),np.cos(0.2)],[-1,0.],[0,-1],[1,0]])*3600
    elif validator.obsmode == 'dgwf':
        startloc = np.array([[0,0.159],[-0.159,0],[0,-0.159],[0.159,0]])*3600.
    elif validator.obsmode == 'gmacs':
        ##in arcsec, and slightly within the boundary; small angle preferred due to vignetting & spot size
        startloc = np.array([[0, 9.2/2], [10.2/2,0], [0, -9.2/2], [-10.2/2, 0] ])*60

    for ns in range(nstars):
        for k in range(4):
            loc = startloc.copy()
            loc[k,:] = gspos[ns,:]
            t = np.ndarray.flatten(loc / 3600.)
            results[k,ns] = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])

    return results

def agwscheck(validator,gsloc,nstars_out=4):
    # gsloc is currently a 4x2 matrix of guide star positions, in arcsec!
    # will be upgraded to accept fewer guide stars
    # obsmode is the observing mode, 'dgnf', 'm3', 'dgwf', 'gmacs'
    
    dims = gsloc.shape
    if len(dims) != 2:
        print('*** gsloc must be an Nx2 array of star coordinates ***')
        return (False,gsloc,[])

    if dims[1] != 2:
        print('*** gsloc must be an Nx2 array of star coordinates ***')
        return (False,gsloc,[])

    nstars = dims[0]
        
    if validator.obsmode == 'dgnf':
        startloc = np.array([[0,0.166],[-0.166,0],[0,-0.166],[0.166,0]])*3600.
    elif validator.obsmode == 'm3':
        startloc = 0.166*np.array([[np.sin(0.2),np.cos(0.2)],[-1,0.],[0,-1],[1,0]])*3600
    elif validator.obsmode == 'dgwf':
        startloc = np.array([[0,0.159],[-0.159,0],[0,-0.159],[0.159,0]])*3600.
    elif validator.obsmode == 'gmacs':
        ##in arcsec, and slightly within the boundary; small angle preferred due to vignetting & spot size
        startloc = np.array([[0, 9.2/2], [10.2/2,0], [0, -9.2/2], [-10.2/2, 0] ])*60
        
    if nstars_out < 4:
        vals = []
        idx = []

        # here, we input 4 or more stars but return fewer than 4 valid stars
        gslocs = np.concatenate((gsloc,startloc),axis=0)

        comb1 = list(itertools.combinations(range(0,nstars),nstars_out))
        comb2 = list(itertools.combinations(range(nstars,nstars+4),4-nstars_out))

        combinations = list(itertools.product(*(comb1,comb2)))
        
        flatten = lambda l: [item for sublist in l for item in sublist]
        for k in range(len(combinations)):
            combinations[k] = flatten(combinations[k])
            
        for comb in combinations:
            for perm in list(itertools.permutations(comb)):
                #import pdb;pdb.set_trace()

                loc = gslocs[perm,:]
                t = np.ndarray.flatten(loc / 3600.)
                check = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
                if check:
                    vals.append(loc)
                    idx.append(np.array(perm))
                    break

        if vals == []:
            return (False,gsloc,idx)
        else:
            return (True,vals,idx)

    if nstars < 4:
        gslocs = np.concatenate((gsloc,startloc),axis=0) 

        ncomb = nCk(4,nstars)
        print(nstars, ncomb)
        combinations = np.array(list(itertools.combinations(range(nstars+4),4)))
        combinations = combinations[0:ncomb,:] # only include the ones with the specified gsloc
        
        for comb in combinations:
            for perm in list(itertools.permutations(comb)):
                loc = gslocs[perm,:]
                t = np.ndarray.flatten(loc/3600.)
                check = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
                print(perm, check)
                if check:
                    gsloc = loc
                    return (True,gsloc,np.array(perm))

        return (False,gsloc,[])

    if nstars == 4:
        permutations = np.array(list(itertools.permutations([0,1,2,3])))
        for perm in permutations:
            loc = gsloc[perm,:]
            t = np.ndarray.flatten(loc / 3600.)

            check = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
            if check:
                gsloc = loc
                return (True,gsloc,np.array(perm))

        return (False,gsloc,[])

    if nstars > 4:
        vals = []
        idx = []
        combinations = np.array(list(itertools.combinations(range(nstars),4)))
        for comb in combinations:
            permutations = np.array(list(itertools.permutations(comb)))
            for perm in permutations:
                loc = gsloc[perm,:]
                t = np.ndarray.flatten(loc / 3600.)
                check = validator.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
                if check:
                    vals.append(loc)
                    idx.append(np.array(perm))
                    break

        if vals == []:
            return (False,gsloc,idx)
        else:
            return (True,vals[0],idx[0])

def test():
    config  = 'dgnf'        
    validator = agwsinit(config)
    
    gsloc = np.array([[0,0.166],[-0.166,0],[0,-0.166],[0.166,0]])*3600
    (success,gspos,idx) = agwscheck(validator,gsloc)
    print(success,idx)
    
    
    gsloc = np.array([[0,-0.166]])*3600
    (success,gspos,idx) = agwscheck(validator,gsloc)
    print(success,idx)
    
    gsloc = np.array([[0,0.166],[0.14,0],[0,-0.166],[0.166,0]])*3600
    (success,gspos,idx) = agwscheck(validator,gsloc)
    print(success,idx)
    
    (success,gspos,idx) = agwscheck(validator,gsloc,nstars_out=3)
    print(success,idx)
