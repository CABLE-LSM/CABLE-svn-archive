#!/usr/bin/env python

"""
Create an integration branch in the shared space.

That's all folks.
"""

__author__ = "Martin De Kauwe"
__version__ = "1.0 (09.03.2019)"
__email__ = "mdekauwe@gmail.com"

import os
import sys
import subprocess
import datetime

root = "https://trac.nci.org.au/svn/cable"
msg = "\"setup integration branch\""

cmd = "svn copy %s/trunk %s/branches/Share/integration -m %s" % \
        (root, root, msg)
error = subprocess.call(cmd, shell=True)
if error is 1:
    raise("Error copying repo")
