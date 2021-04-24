import numpy as np


class quiet(object):

    bands = ['Q', 'W']

    detectors = {'Q': np.arrange(0,18), 'W':np.arrange(0,90)}

    diodes = [0, 1, 2, 3]
    

    @staticmethod
    def instrument_filename(version):
        return 'quiet_instrument_v' + str(version) + '.h5'

    def centFreq(detName):

        if detName[0] == 'Q':
            return NotImplementedError
        if detName[0] == 'W':
            
