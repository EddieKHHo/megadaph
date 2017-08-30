#!/usr/bin/env python

import os, sys, unittest
from glob import glob

if __name__ == '__main__':                                                      
    # Add submodules to path
    sys.path.append('../modules/')

    test_modules = glob("test_*.py")

    # Strip .py from name and convert to string
    for i in range(len(test_modules)):
        test_modules[i] = test_modules[i][:-3]
    
    suite = unittest.TestSuite()

    for t in test_modules:
        suite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))

    unittest.TextTestRunner().run(suite)
