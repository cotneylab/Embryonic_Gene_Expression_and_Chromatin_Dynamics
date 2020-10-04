Files used to construct WGCNA from RNA-seq of craniofacial samples and NCCs cell line. 

Used the Step-by-Step network from Horvath Tutorial.

WGCNA_step.sh: slurm script that runs wgcna_step.R and produces a TOM and gene tree data.

MergeModules.R: contains functions to be used once WGCNA_step.sh has run. These functions merge the eigenmodules and continue to follow the step-by-step network construction from the Horvath Tutorial.

Build_Intermod_Network.R: R script that was used to create the inter-module plot from the outputs of MergeModules.R
