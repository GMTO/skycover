from agwsprobes import *
import sys

dgnf = agwsinit('dgnf')
m3   = agwsinit('m3')
dgwf = agwsinit('dgwf')
gmacs = agwsinit('gmacs')

# check single star
#input to agwscheck needs to be in arcsec, but we like to think in arcmin
testpos = 60* np.transpose(np.array(([0,7.2])).reshape(2,1))
(success,loc,idx) = agwscheck(m3,testpos)
assert not success

(success,loc,idx) = agwscheck(dgnf,testpos)
assert success

(success,loc,idx) = agwscheck(dgwf,testpos)
assert success

(success,loc,idx) = agwscheck(gmacs,testpos)
assert not success

testpos = 60* np.transpose(np.array(([0,6.9])).reshape(2,1))
(success,loc,idx) = agwscheck(gmacs,testpos)
assert success

#sys.exit()

#input to agwscheck needs to be in arcsec, but we like to think in arcmin
testpos = 60* np.array(([0,6.9],[-6.9,0], [0,-6.9], [6.9, 0]))

(success,loc,idx) = agwscheck(dgnf,testpos)
assert success

(success,loc,idx) = agwscheck(m3,testpos)
assert not success

(success,loc,idx) = agwscheck(dgwf,testpos)
assert not success

(success,loc,idx) = agwscheck(gmacs,testpos)
assert success
