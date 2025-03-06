# South Edinburgh Perfusion+ Analysis Library (SEPAL)

**Please note: This library is also hosted in the [OSIPI DCE-DSC-MRI_CodeCollection repository](https://github.com/OSIPI/DCE-DSC-MRI_CodeCollection), where unit tests and perfusion code by other authors can also be found.**

Python library for simulating and fitting DCE- and other quantitative MRI data. It permits arbitrary combinations of pulse sequence, pharmacokinetic model, water exchange model, etc. The code is a work-in-progress, has not been extensively tested and is not recommended or approved for clinical use.

Created 28 September 2020  
@authors: Michael Thrippleton  
@email: m.j.thrippleton@ed.ac.uk  
@institution: University of Edinburgh, UK
---

### Installation:
pip install sepal

### Use:
Most functionality is demonstrated in Jupyter notebook format in ./demo 

### Functionality:
- Enhancement-to-concentration conversion (assuming fast water exchange)
- Fit tissue concentration using pharmacokinetic model
- Fit signal enhancement using pharmacokinetic model
- Pharmacokinetic models: steady-state, Patlak, extended Tofts, Tofts, 2CXM, 2CUM
- Patlak fitting with multiple linear regression
- AIFs: including patient-specific (measured), Parker, bi-exponential Parker, Georgiou
- Fitting free AIF time delay parameter
- Relaxivity models: linear
- Signal models: spoiled gradient echo, inversion-recovery spin-echo
- Water exchange models: FXL, NXL, NXL_be
- T1 fitting using variable flip angle method, IR-SPGR, DESPOT1-HIFI and inversion recovery 
- T2(*) fitting for multi-TE acquisitions
- MTR and MTSat calculation

### Not yet implemented/limitations:
- Additional pharmacokinetic models (add by inheriting from PkModel class)
- Additional relaxivity models (add by inheriting from CRModel class)
- Additional water exchange models, e.g. 3S2X, 2S1X (add by inheriting from WaterExModel class)
- Additional signal models (add by inheriting from SignalModel class)
- R2/R2* effects not included in fitting of enhancement curves (but is included for enhancement-to-concentration conversion)
- Compartment-specific relaxivity parameters/models
- Fitting water exchange parameters
---

### Updates
Release 1.3.2 - further comments added to dce_fit demo notebook  
Release 1.3.1 - Checking ve+vp<=1 can now be disabled  
Release 1.2.1 - "Georgiou" AIF added to *aifs* module  
Release 1.1.1 - *roi_measure* exclude NaNs when calculating percentiles  
Release 1.1.0 - *roi_measure* modified to generate more ROI statistics  
Release 1.0.3 - Add exception handling for some zero/negative inputs.  
Release 1.0.2 - Changed AIF interpolation method to linear to further reduce oscillations.  
Release 1.0.1 - Changed AIF interpolation method to avoid oscillations. Added demo notebook on interpolation.

