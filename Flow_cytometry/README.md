# FLOW CYTOMETRY ANALYSIS

Don't forget that these analyses come from one replicate of a 3-replicate experiment, and as such, they are meant to be representative of the results you see consistently across replicates. These type of representative figures are supposed to be accompained by a summary plot describing all measures obtained from the 3 replicates (such as a bar plot with percentages, with a statical analysis stating whether the results are significant).

All scripts require as input a Metadata table, indicating condition/sample and its corresponding file name, like this:

| condition | file_name             |
|-----------|------------------------|
| control   | example_control.csv   |
| wildtype  | example_wildtype.csv  |
| knockout  | example_knockout.csv  |


**Instructions to run:**
For all scripts, you can print a help message by typing:
```bash
Rscript <script name> --help
```
This will tell you the required arguments and their types.



## Script 1 - 2-AXIS PLOT - - - - 

The script 1 takes as input a metadata table pointing to csv tables of fluorescence, assuming compensation and gating of singlets (single cells) has been done inside the software of the cytometer (such as Accuri C6). You can follow this [video tutorial](https://www.youtube.com/watch?v=TdzVxcGLVzU&t=37s) to see how the export of the data was done (minutes from 0 to 4:16).


Example command:
```bash
Rscript 1_flow_cytometry_from_csv.R example_metadata.txt FL1.A FL4.A NKX2.5 TNNT2
```

Output: 
the script generate a pdf file named cytometry_hexplot.pdf (you can find the plot generated using the example data provided whithin this folder).

<img width="800" alt="hexplot" src="https://github.com/user-attachments/assets/eaac6fa0-beea-4d4a-8365-7c35f6e0bdb0" />


*The example data for this analysis corresponds to cardiomyocytes derived from human pluripotent stem cells, at day 15 (previous to metabolic selection).*

## Script 2 - 2-AXIS PLOT - - - - 

The starting point of this script is a metadata table pointing to .fcs files that are exported from the cytometer software (such as Accuri C6).   
NOTES:  
First, you will most likely have to change the parameters of the gates created, as gates change for different cell types. You can always add more gates if needed following the code logic on the first 3 gates.  
Second, I recommend you run this script from RStudio. For some reason, the first time running from command line there are some adjustments to make. 


Example command:
```bash
Rscript 2_flow_cytometry_from_fcs.R example_metadata_fcs.txt yes no FL1-A FL2-A 
```

Output: 
the script generates 3 pdf files:  
-the first two to visualize where the filtering gates are being placed (gate1_cells and gate2_singlets)  
-the last one to point and quantify the percentage of a specific population (Population_of_study)  
Also a quality report is generated when quality-filter argument is set to "yes".

<img width="850" alt="population_of_study" src="https://github.com/user-attachments/assets/105d2a69-040e-4451-bdd4-416189111b56" />  

*The example data for this analysis corresponds to day 3.5 of cardiac differentiation, when a Mesoderm Progenitor population forms (compared to pluripotent stem cells).*  
