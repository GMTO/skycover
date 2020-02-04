import subprocess as sub
import numpy as np


class validator:

    # Set up the pipes to the agwsvalid executable
    # mode should be 'dgnf', 'm3', or 'dgwf'
    def __init__(self, mode):
        self.pipes = sub.Popen(['agwsvalid',  '--'+mode, '--bool'], stdin=sub.PIPE, stdout=sub.PIPE)


    # Write probe positions 
    def check(self, probeax, probeay, probebx, probeby, probecx, probecy, probedx, probedy):
        probestr = "%f %f %f %f %f %f %f %f" % (probeax, probeay, probebx, probeby, probecx, probecy, probedx, probedy) + '\n'
        self.pipes.stdin.write(bytes(probestr, 'UTF-8'))
        self.pipes.stdin.flush()
        r = self.pipes.stdout.readline().decode('UTF-8').strip()
        if (str(r) == "1"):
            return True
        else:
            return False

