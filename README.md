# HRTF selection
Scripts to select good and bad non-individual Head-Related Transfer Functions (HRTFs) based on the [barumerli2023](http://amtoolbox.org/amt-1.5.0/doc/models/barumerli2023.php) auditory model [[1](#barumerli2023)]. The procedure is detailed in [[2](#daugintis2023)].

 - Requires Auditory Modelling Toolbox (AMT): http://amtoolbox.org.
 - Python script requires the following libraries:
    - numpy
    - pandas
    - seaborn
    - matplotlib
    - scipy
    - pathlib

 - [test_hrtfs](https://github.com/rapolasd/HRTF-selection/tree/main/test_hrtfs) folder contains a few HRTFs from the [SONICOM](https://www.axdesign.co.uk/tools-and-devices/sonicom-hrtf-dataset) database [[3](#engel2023)] that can be used to test the selection method.

## Example use:

   1. In MATLAB, load AMT: ```amt_start```
   2. Run: ```daugintis2023_hrtfSelection({'P0001_Windowed_48kHz.sofa','P0019_Windowed_48kHz.sofa'},'./test_hrtfs');```
      - To run the prediction for all the HRTFs in the folder as subjects, an empty cell array can be supplied, e.g.: ```daugintis2023_hrtfSelection({},'./test_hrtfs');```
   3. After the MATLAB function is executed, run Python script: ```python daugintis2023_hrtfSelection.py```
      - ```no_plots``` flag can be added to prevent the script from plotting the modelled distributions and the selection matrix.
   4. HRTF selection is stored in a table **selected_HRTFs.csv**

## Notes

The MATLAB function creates a new folder, named ```model_predictions```, where it stores the tables of modelled errors for each subject, separately. The Python script goes through every file in this folder and selects the good/bad HRTFs for each subject. The selection table is saved in the parent folder. If plotting is allowed (default behaviour), the script creates  the ```selection_figures``` folder and stores the plots there. Both MATLAB and Python functions overwrite the old files in the corresponding folders with the same names, so make sure to save them elsewhere if you want to keep them. Furthermore, for the plotting to work, the files in the ```model_prediction``` folder must be from the same model run.

### Optional usage

Optionally, MATLAB function can take column arrays of corresponding target azimuths and elevations to be used in the ```barumerli2023``` model, e.g.: ```daugintis2023_hrtfSelection({},'./test_hrtfs', 'target_az', [-5; 0; 5], 'target_el', [10; 10; 10]);```. This overwrites the default selection of directions within ±30° azimuth and ±11.5° elevation. **NB**: the Python script is not adapted to select HRTFs based on directions other than the default ones and might show inconsistent behaviour.

## References

<a id="barumerli2023">[1]</a> R. Barumerli, P. Majdak, M. Geronazzo, D. Meijer, F. Avanzini, and R. Baumgartner, "A Bayesian model for human directional localization of broadband static sound sources," *Acta Acust.*, vol. 7, p. 12, 2023, doi: [10.1051/aacus/2023006](https://doi.org/10.1051/aacus/2023006).

<a id="daugintis2023">[2]</a> R. Daugintis, R. Barumerli, L. Picinali, and M. Geronazzo, “Classifying Non-Individual Head- Related Transfer Functions with A Computational Au- ditory Model: Calibration And Metrics,” in *Proc. IEEE Int. Conf. Acoust. Speech Signal Process. (ICASSP)*, 2023, doi: [10.1109/ICASSP49357.2023.10095152](https://doi.org/10.1109/ICASSP49357.2023.10095152).

<a id="engel2023">[3]</a> I. Engel, R. Daugintis, T. Vicente, A. O. T. Hogg, J. Pauwels, A. J. Tournier, and L. Picinali, “The SONICOM HRTF dataset,” *J. Audio Eng. Soc. (AES)*, vol. 71, no. 5, pp. 241–253, 2023, doi: [10.17743/jaes.2022.0066](https://doi.org/10.17743/jaes.2022.0066).
