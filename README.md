# South Edinburgh Perfusion Analysis Library (SEPAL)

**Please note: This library is also hosted in the [OSIPI DCE-DSC-MRI_CodeCollection repository](https://github.com/OSIPI/DCE-DSC-MRI_CodeCollection), where unit tests and perfusion code by other authors can also be found.**

Python library for simulating and fitting DCE-MRI data. It permits arbitrary combinations of pulse sequence, pharmacokinetic model, water exchange model, etc. The code is a work-in-progress, has not been extensively tested and is not recommended or approved for use.

Created 28 September 2020
@authors: Michael Thrippleton  
@email: m.j.thrippleton@ed.ac.uk  
@institution: University of Edinburgh, UK

### Functionality:
- Enhancement-to-concentration conversion (assuming fast water exchange)
- Fit tissue concentration using pharmacokinetic model
- Fit signal enhancement using pharmacokinetic model
- Pharmacokinetic models: steady-state, Patlak, extended Tofts, Tofts, 2CXM, 2CUM
- AIFs: patient-specific (measured), Parker, bi-exponential Parker
- Fitting free AIF time delay parameter
- Relaxivity models: linear
- Signal models: spoiled gradient echo
- Water exchange models: FXL, NXL, NXL_be
- T1 fitting using variable flip angle method, IR-SPGR and DESPOT1-HIFI

### Not yet implemented/limitations:
- Generally untested. Not optimised for speed or robustness.
- Additional pharmacokinetic models (add by inheriting from PkModel class)
- Additional relaxivity models (add by inheriting from CRModel class)
- Additional water exchange models, e.g. 3S2X, 2S1X (add by inheriting from WaterExModel class)
- Additional signal models (add by inheriting from SignalModel class)
- R2/R2* effects not included in fitting of enhancement curves (but is included for enhancement-to-concentration conversion)
- Compartment-specific relaxivity parameters/models
- Fitting free water exchange parameters
- Special model implementations, e.g. linear and graphical versions of Patlak model

### TODO:
- linear Patlak model
- fast C calculation for SPGR with r2=0
- inversion recovery T1 measurement
