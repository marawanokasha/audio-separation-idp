# Audio-separation-idp

Use the main_* files in the following directories to run the code, a directory 
data should be found outside of the current working directory containing the audio
data to run on:

* grid_code
* sisec_code
* mass_code

The scatt_functions folder is provided by Joan Bruna and is available from 
here: https://github.com/joanbruna/scattnmf. It contains functions for applying
scattering transform and phase recovery for audio reconstruction

The spams library found in the lib folder has to be recompiled on every machine.
We were only able to get it running on Ubuntu using gcc
