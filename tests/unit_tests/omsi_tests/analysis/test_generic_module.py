"""
Basic testing for core functionality of the omsi_file package.
This includes functionality across many of the modules.
"""
import unittest
from omsi.workflow.executor.greedy_executor import greedy_executor

SCRIPT =  \
"""
from omsi.analysis.generic import analysis_generic
import numpy as np

msidata = np.arange(10*10*20).reshape(10,10,20)

# Define a simple function to compute the total intensity image
def total_intensity(msidata, axis=2):
    import numpy as np
    return np.sum(msidata, axis=axis)

# Define a simple function to normalize an MSI data cube by per-spectrum normalization factors
def normalize_intensities(msidata, normfactors):
    import numpy as np
    return msidata / normfactors[:,:,np.newaxis]

# Define compute of total intensity image
a2 = analysis_generic.from_function(analysis_function=total_intensity,
                                    output_names=['total_intensities'])
a2['msidata'] = msidata

# Define the normalization of the peak cube
a3 = analysis_generic.from_function(normalize_intensities)
a3['msidata'] = msidata
a3['normfactors'] = a2['total_intensities']
"""

class test_generic(unittest.TestCase):

    def setUp(self):
        from omsi.shared.log import log_helper
        log_helper.set_log_level('WARNING')

    def tearDown(self):
        pass

    def test_from_scripts(self):
        try:
            driver = greedy_executor.from_scripts(SCRIPT)
        except:
            self.fail("Creation  of the workflow fomr script failed")
        try:
            driver.execute()
        except:
            self.fail("Execution of the workflow failed.")



if __name__ == '__main__':
    unittest.main()





