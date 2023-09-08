# ReadME

This repository contains the code necessary to recreate the figures in the paper "Diagnostic Testing and the evolution of detection avoidance by pathogens"

Here we investigate the evolution of detection avoidance in response to different diagnostic testing regimes. 

Before running any of the python files please ensure the src files are correctly compiled, and the necessary data folders exist. The data folders can be created by running setup.sh from the main folder of the repository.

Figure 1 from the manuscript is created using the file "heatmaps_perfect_both.py", which uses the derived fitness gradients and multiple wrapper files which solve the ordinary differential equations of the various models.

Figure 2 in the manuscript is created running the file "Both_models_evolutionary_analysis_perfect.py", which uses the same wrapper files and fitness functions.

Figure 3 is created using the file "Both_models_evolutionary_analysis_imperfect.py".

Figure 4 is created using the file "Both_gillespie_imperfect_large_pop.py".

Figure S1 from the supplementary material is created using the file "Both_gillespie_imperfect.py"

Figure S2 from the supplementary material is created using the file "Checking_Model_B_over_larger_range_of_zeta_values.py". 
This may produce a different result from the manuscript as the number of repeats is determined by the available nubmer of cpu cores accessible to python. 
Warning, this file produces many large .csv files, each roughly 1GB in size.