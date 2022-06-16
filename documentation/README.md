# GROMORA Documentation

This folder contains a specific documentation for the GROMORA project. 

This documentation tries to make a comprehensive description of the GROMORA
routine which is constantly evolving. Therefore, it might sometimes contain
error or lack updates. In this case, please refer to the author.

It completes the 2 user guides located at
[UserGuideGROMORA](https://git.iap.unibe.ch/IAP_MCH/UserGuideGROMORA.git) also on
the IAP git server.

For the calibration and integration routine, the following documentation files are the most important:

* [main](main.md): describes the main script for launching GROMORA calibration and integration.
* [run_calibration](run_calibration.md): describes the calibration sub-routine.
* [run_integration](run_integration.md): describes the integration sub-routine.
* [calibrationTool](calibrationTool.md): describes the *calibrationTool* structure.
* [level1](level1.md): describes the structure and content of the calibrated (level 1a) and integrated (level 1b) output files.
* [quality_control_calibration](quality_control_calibration.md): describes the quality control done during the calibration. 

In [example_level1](example_level1), you can find some examples of level 1a and
1b netCDF files as well as the standard plots for calibrated and integrated
data. 

For the retrieval routine, the following documentation files are available:
* [retrieval_param](retrieval_param.md): describes the *retrieval_param* dictionary.

In additon, an example of level 2 daily file and diagnostics plots can be seen in [example_level2](example_level2). 


