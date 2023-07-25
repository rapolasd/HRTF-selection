# HRTF selection
Scripts to select a good and a bad non-individual HRTFs based on an auditory model. 

 - Requires Auditory Modelling Toolbox (AMT): http://amtoolbox.org
 - Python script requires the following libraries:
    - numpy
    - pandas
    - seaborn
    - matplolib
    - scipy
    - pathlib
    - os
    - warnings

## Example use:

   1. In MATLAB, load AMT: ```amt_start```
   2. Run: ```daugintis2023_hrtfSelection({'P0001_Windowed_48kHz.sofa','P0019_Windowed_48kHz.sofa'},'./test_hrtfs');```
      - To run the prediction for all the HRTFs in the folder as subjects, an empty cell array can be supplied, e.g.: ```daugintis2023_hrtfSelection({},'./test_hrtfs');```
   3. After MATLAB function is executed, run Python script: ```python daugintis2023_hrtfSelection.py``` 
      - ```no_plots``` flag can be added to prevent the script from plotting the modelled distributions and the selection matrix.
   4. HRTF selection is stored in a table **selected_HRTFs.csv**

## Notes

The MATLAB function creates a new folder, named ```model_predictions```, where it stores the tables of modelled errors for each subject, separately. The python script goes through every file in this folder and selects the good/bad HRTFs for each subject. The selection table is saved in the parent folder. If plotting is allowed (default behaviour), the script creates  ```selection_figures``` folder and stores the plots there. Both MATLAB and Python functions overwrite the old files in the corresponding folders with the same names, so make sure to save them elsewhere if you want to keep them. Furthermore, for the plotting to work, the files in the ```model_prediction``` folder must be from the same model run.
