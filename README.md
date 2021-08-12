# Python code for fitting DCE-MRI data
This repository contains code for simulating and fitting DCE-MRI data. It is primarily an exercise in writing code that is reasonably general, e.g. permits arbitrary combinations of sequence, pharmacokinetic model, water exchange model, etc.  
This repository is a work-in-progress, has not been extensively tested and is not recommended or approved for use.

### Functionality:
- Enhancement-to-concentration conversion (assuming fast water exchange)
- Fit tissue concentration using pharmacokinetic model
- Fit signal enhancement using pharmacokinetic model
- Pharmacokinetic models: steady-state, Patlak, extended Tofts, Tofts, 2CXM, 2CUM
- AIFs: patient-specific (measured), Parker, bi-exponential Parker
- Relaxivity models: linear
- Signal models: spoiled gradient echo
- Water exchange models: FXL, NXL, NXL_be
- T1 fitting using variable flip angle method

### Not yet implemented/limitations:
- Generally untested. Not optimised for speed or robustness.
- Additional pharmacokinetic models (add by inheriting from pk_model class)
- Non-linear relaxivity models (add by inheriting from c_to_r_model class)
- Additional AIF functions (add by inheriting from aif class)
- Additional water exchange models, e.g. 3S2X, 2S1X (add by inheriting from water_ex_model class)
- Additional signal models (add by inheriting from signal_model class)
- R2/R2* effects not included in fitting of enhancement curves (but is included for enhancement -> concentration)
- Compartment-specific relaxivity parameters/models
- Fitting a time delay
- Fitting water exchange parameters
- Special model implementations, e.g. linear(ised) versions of pharmacokinetic models.
- T1 fitting using other techniques

### To do list (short-term):
- restructure and add to osipi