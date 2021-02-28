import time
import numpy as np
from scipy.optimize import root, minimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from models import interpolate_time_series, pkp_to_c, pkp_to_c_2, pkp_to_c_3

plt.close('all')

c_p_aif= np.array([-2.19962928e-03,  2.15229337e-03,  5.00234997e-05,  5.15988989e-03,
        2.70068054e-01,  9.26311987e-01,  2.64564450e+00,  3.00162108e+00,
        2.84837762e+00,  2.51718248e+00,  2.10865919e+00,  1.79821298e+00,
        1.60299221e+00,  1.49910023e+00,  1.41013864e+00,  1.32236252e+00,
        1.26104970e+00,  1.20402978e+00,  1.16223861e+00,  1.12146056e+00,
        1.07728176e+00,  1.04559582e+00,  1.01135273e+00,  9.70996005e-01,
        9.61448602e-01,  9.26552385e-01,  9.11220125e-01,  9.17914712e-01,
        8.98841757e-01,  8.79464485e-01,  8.66798051e-01,  8.52346341e-01])

n = c_p_aif.size
t_res = 39.62
t_res_interp = t_res/3
t = np.linspace(t_res/2., t_res*(n-0.5), n)
t_interp=interpolate_time_series(t_res_interp,t)
print(t)
print(t_interp)


#interpolate AIF
f=interp1d(t,c_p_aif,kind='quadratic', bounds_error=False, fill_value=c_p_aif[0])
c_p_aif_interp=f(t_interp)
#plt.plot(t,c_p_aif,'.')
#plt.plot(t_interp,c_p_aif_interp,'-')

hct=0.42
irf_model='patlak'
pk_pars = {'vp': 0.01, 'ps': 1e-3, 've': 0.2}

st=time.time()
for m in range(1000):
    c, ct = pkp_to_c(t, c_p_aif, pk_pars, hct, irf_model, options=None)
print(time.time()-st)
st=time.time()
for m in range(1000):
    c2, ct2 = pkp_to_c_2(t, t_interp, c_p_aif_interp, pk_pars, hct, irf_model, options=None)
print(time.time()-st)
st=time.time()
for m in range(1000):
    c3, ct3 = pkp_to_c_3(t, t_interp, c_p_aif_interp, pk_pars, hct, irf_model, options=None)
print(time.time()-st)

print(np.sqrt(np.sum((ct-ct2)**2)/ct.size))
print(np.sqrt(np.sum((ct3-ct2)**2)/ct.size))

plt.plot(t,ct,'-')
plt.plot(t,ct2,'o')
plt.plot(t,ct2,'x')
# ax.plot(t_interp,'-')

# #create interpolation function and get interpolated signal 
# f=interp1d(t1,s1,kind='quadratic', bounds_error=False, fill_value=s1[0])
# s2=f(t2)
# plt.plot(t2,s2)
# print(s2)


# #c, ct = pkp_to_c(t, c_p_aif, pk_pars, hct, irf_model, options=None):




# #generate example signal
# s1=[100, 105, 95, 115, 145, 180, 200, 190, 187, 180]
# ax=plt.plot(t1,s1,'o')

# #create interpolation function and get interpolated signal 
# f=interp1d(t1,s1,kind='quadratic', bounds_error=False, fill_value=s1[0])
# s2=f(t2)
# plt.plot(t2,s2)
# print(s2)
# return t_interp