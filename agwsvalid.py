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

    # Write probe positions 
    def check(self, probeax, probeay, probebx, probeby, probecx, probecy, probedx, probedy):
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

