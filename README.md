# HRTF selection
Scripts to select a good and a bad non-individual HRTFs based on an auditory model. 

 - Requires Auditory Modelling Toolbox (AMT): http://amtoolbox.org
 - Python script requires the following libraries:
    - numpy
    - pandas
    - seaborn
    - matplolib
    - scipy

## Example use:

   1. In MATLAB, load AMT: ```amt_start```
   2. Run: ```daugintis2023_hrtfSelection({'P0001_Windowed_48kHz.sofa','P0019_Windowed_48kHz.sofa'},'./test_hrtfs');```
   3. After MATLAB function is executed, run Python script: ```python daugintis2023_hrtfSelection.py``` 
      - ```no_plots``` flag can be added to prevent the script from plotting the modelled distributions and the selection matrix.
   4. HRTF selection is stored in a table **selected_HRTFs.csv** 
