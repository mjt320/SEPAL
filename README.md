# Python code for fitting DCE-MRI data
Note: this repository is a work-in-progress, un-tested and not recommended for use.

### Functionality:
- Enhancement-to-concentration conversion
- Fit concentration using pharmacokinetic model
- Fit enhancement using pharmacokinetic model
- Pharmacokinetic models: steady-state, Patlak, extended Tofts, Tofts, 2CXM, 2CUM
- Individual and Parker AIFs
- Linear relaxivity model
- SPGR signal model
- Water exchange models: FXL, NXL, NXL_be

### Not yet implemented/limitations:
- Generally untested!
- Other pharmacokinetic models (add by inheriting pk_model class)
- Non-linear relaxivity models add by inheriting c_to_r_model class)
- Other population AIFs (add by inheriting aif class)
- Other water exchange models, e.g. 3S2X, 2S1X (add by inheriting water_ex_model class)
- Signal models for other sequences (add by inheriting signal_model class)
- Model constraints
- R2* effects not included
- Compartment-specific relaxivity parameters/models
- Fitting a time delay
- Fitting water exchange parameters
- Raising exceptions e.g. for non-physical values
- Finding global minimum
- Special model implementations, e.g. linear(ised) versions
- Generalise T1 fitting code to use other methods than VFA, e.g. IR, DESPOT1-HIFI
- Find global minimum

### To do list (short-term):
- docstrings and comments
- finish demo notebooks (water_ex, signal_models, relax)
- optimise tolerances, iterations etc.


