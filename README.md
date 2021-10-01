**multiple_animals_movement_model**


This repository contains the code to replicate the results of the paper **Multiple animals movement model**.

The code is written in *R* and *Julia*. To estimate the models you have to run the *Julia* code  *ModelEstimate.jl*. Inside the script, the parameter **DIR_CODE** must be equal directory path where this repository is downloaded. A file *Rdata. is saved in the directory *output*. Notice that the code need the library *HierarchicalMultivariateAnimalMovement* that is part of this repository.
The parameter **imod** is used to obtain the results of the three models shown in the paper:

* **imod = 1** model M3
* **imod = 2** model M1
* **imod = 3** model M2


Once the model is estimated, the plots and tables present in the paper, with the exception of the one with the models comparisons, can be reproduced using the *R* code *PlotsAndTables.R*. Figures are saved in the directory *output/Figures/*, while the tables are printed in the *R* console. As for *ModelEstimate.jl*., **DIR_CODE** must be set as the directory path where this repository is downloaded.
