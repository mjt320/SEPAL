# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 19:50:31 2021

@author: mthripp1
"""

import time
import numpy as np

st=time.time()
a=np.linspace(1,100,1000)
b=np.linspace(1,100,1000)
for n in range(1000):
    np.convolve(a,b,mode='full')
print(time.time()-st)