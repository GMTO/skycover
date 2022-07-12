import numpy as np
import sys
import agwsvalid

dgnf = agwsvalid.validator('dgnf')
m3   = agwsvalid.validator('m3')
dgwf = agwsvalid.validator('dgwf')
gmacs = agwsvalid.validator('gmacs')

#input to agwscheck needs to be in degree, but we like to think in arcmin
t0 = np.array([0 , 1.0, -1.0, 0, 0, -1.0, 1.0, 0])/60. #1 arcmin circle, in deg

#the ordering the points matter.
#because agwsvalid always uses probe0 to try to reach the first point, probe1 for 2nd point, and so on
# probe 0 rotation center is at 12 o'clock
# probe 1 rotation center is at 9  o'clock
# probe 0 rotation center is at 6 o'clock
# probe 1 rotation center is at 3  o'clock

## --------------------------dgnf
t = t0 * 7.2
print('dgnf 7.2 arcmin circle, expect success')
success = dgnf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

t = t0 * 6.4
print('dgnf 6.4 arcmin circle, expect failure')
success = dgnf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

t = t0 * 10.1
print('dgnf 10.1 arcmin circle, expect failure')
success = dgnf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

## --------------------------m3
t = t0 * 7.2
print('m3 7.2 arcmin, evenly-spaced, expect failure')
success = m3.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

#sys.exit()

t0 = np.array([0.707, 0.707, -1.0, 0, 0, -1.0, 1.0, 0])/60. #1 arcmin circle, in deg
t = t0 * 7.2
print('m3 7.2 arcmin, top star shift to right, expect success')
success = m3.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

t = t0 * 9.9
print('m3 9.9 arcmin, top star shift to right, expect success')
success = m3.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

t = t0 * 10.1
print('m3 10.1 arcmin, top star shift to right, expect failure')
success = m3.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

t = t0 * 6.6
print('m3 6.6 arcmin, top star shift to right, expect sucess')
success = m3.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

t = t0 * 6.4
print('m3 6.4 arcmin, top star shift to right, expect failure')
success = m3.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

## --------------------------dgwf

t0 = np.array([0 , 1.0, -1.0, 0, 0, -1.0, 1.0, 0])/60.

t = t0 * 7.2
print('dgwf 7.2 arcmin, evenly-spaced, expect failure')
success = dgwf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

t = t0 * 8.0
print('dgwf 8.0 arcmin, evenly-spaced, expect success')
success = dgwf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

#sys.exit()

t = t0 * 10.1
print('dgwf 10.1 arcmin, evenly-spaced,, expect failure')
success = dgwf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

## --------------------------gmacs

t = t0 * 7.1
print('gmacs 7.1 arcmin, evenly-spaced, expect failure')
success = gmacs.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success

t = t0 * 6.9
print('gmacs 6.9 arcmin, evenly-spaced, expect success')
success = gmacs.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

t = np.array([5.1 , 4.6, -5.1, 4.6, -5.1, -4.6, 5.1, -4.6])/60.
print('gmacs, just outside the corners of the obscuration')
success = gmacs.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert success

t = np.array([5.0 , 4.5, -5.0, 4.5, -5.0, -4.5, 5.0, -4.5])/60.
print('gmacs, just inside the corners of the obscuration')
success = gmacs.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
assert not success
