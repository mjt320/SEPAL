# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 17:06:01 2021

@author: mthripp1
"""

import numpy as np
from scipy.optimize import root, minimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

t1=np.linspace(1.,10.,10)-0.5
print(t1)
max_t = np.max(t1)


delta_t_interp_req = 0.25
n_interp = np.round(max_t/delta_t_interp_req + 0.5).astype(int)
delta_t_interp_act = max_t/(n_interp-0.5)
print(delta_t_interp_act)

t2=np.linspace(0.5*delta_t_interp_act, max_t, num=n_interp)
print(t2)

#y1=10+0.05*t1**2-0.005*t1**3
y1=[100, 105, 95, 115, 145, 180, 200, 190, 187, 180]
ax=plt.plot(t1,y1,'o')

f=interp1d(t1,y1,kind='quadratic', bounds_error=False, fill_value=y1[0])

y2=f(t2)
plt.plot(t2,y2)