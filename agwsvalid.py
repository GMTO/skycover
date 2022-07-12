import subprocess as sub
import numpy as np

class validator:

    # Set up the pipes to the agwsvalid executable
    # mode should be 'dgnf', 'm3', or 'dgwf'
    def __init__(self, mode, silent=False):
        if silent: 
            self.pipes = sub.Popen(['agwsvalid',  '--'+mode, '--boolonly'], stdin=sub.PIPE, stdout=sub.PIPE)
        else:
            self.pipes = sub.Popen(['agwsvalid',  '--'+mode, '--bool'], stdin=sub.PIPE, stdout=sub.PIPE)

    def close(self):
        self.pipes.stdin.close()
        self.pipes.stdout.close()
        
    def check(self, probeax, probeay, probebx, probeby, probecx, probecy, probedx, probedy,**kwargs):
        probestr = "%f %f %f %f %f %f %f %f" % (probeax, probeay, probebx, probeby, probecx, probecy, probedx, probedy) + '\n'
        try:
            self.pipes.stdin.write(bytes(probestr, 'UTF-8'))
            self.pipes.stdin.flush()
            r = self.pipes.stdout.readline().decode('UTF-8').strip().split()
        except:
            return False
        
        self.shadowfrac = float(r[1])

        if (str(r[0]) == "1"):
            return True
        else:
            return False

def test():
        
    dgwf = validator('dgwf',silent=True)
    gsloc = np.array([[0,0.159],[-0.159,0],[0,-0.159],[0.159,0]])*0.9
    t = np.ndarray.flatten(gsloc)
    check = dgwf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
    print(check)
    print(dgwf.shadowfrac)
    
    dgnf = validator('dgnf',silent=True)
    gsloc = np.array([[0,0.159],[-0.159,0],[0,-0.159],[0.159,0]])*0.9
    t = np.ndarray.flatten(gsloc)
    check = dgnf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
    print(check)

    dgnf = validator('gmacs',silent=True)
    gsloc = np.array([[0, 9.2/2], [10.7/2,0], [0, -9.2/2], [-10.7/2, 0] ])/60 #in degree
    t = np.ndarray.flatten(gsloc)
    check = dgnf.check(t[0],t[1],t[2],t[3],t[4],t[5],t[6],t[7])
    print(check)
