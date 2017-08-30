# EPIFIL
Matlab code used to fit the lymphatic filariasis model EPIFIL to baseline infection data and simulate mass drug administration and vector control interventions

BaselineFitting: contains the code to run the baseline fitting routine given mf age profile and ABR data

Interventions: contains the code to simulate interventions using the fitted parameters

CommonFunctions: contains several small functions which are used repeatedly in both the fitting and intervention code

The entry point to each routine is a file in each folder that starts with the word "Main" (Main_ToRunBaselineFitting.m and Main_ToRunIntv.m). One output .mat file is saved by running the baseline fitting procedure which contains the fitted parameters and model-generated age prevalence curves, and another .mat file is saved by running the intervention procedure which contains the simulated mf prevalence over time. 
