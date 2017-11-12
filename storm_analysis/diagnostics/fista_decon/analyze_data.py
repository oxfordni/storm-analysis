#!/usr/bin/env python
"""
Analyze test data using FISTA Spiner

Hazen 11/17
"""
import glob
import os
import time

import storm_analysis.spliner.spline_analysis as sp_ana

dirs = sorted(glob.glob("test*"))

for a_dir in dirs:
    print()
    print("Analyzing:", a_dir)
    print()
    
    mlist = a_dir + "/test_mlist.bin"

    # Remove stale results, if any.
    if os.path.exists(mlist):
        os.remove(mlist)

    # Run analysis.
    start_time = time.time()
    sp_ana.analyze(a_dir + "/test.dax", mlist, "fdecon.xml")
    stop_time = time.time()

    # Save timing results.
    print("Analysis completed in {0:.2f} seconds".format(stop_time - start_time))

    with open(a_dir + "/timing.txt", "w") as fp:
        fp.write(str(stop_time - start_time) + "\n")

