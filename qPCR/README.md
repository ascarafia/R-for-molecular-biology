# qPCR and statistical analysis

This script was developed for the analysis of Step One Plus results, and following the concept of Ariel Waisman [video tutorials](https://www.youtube.com/watch?v=hPO7ptWOT7M). It is a 4 part series on RT-qPCR data analysis I highly recommend.



## General instructions



1- Do not use underscores for your gene names.
2- Keep names of genes and samples consistent between different runs.
3- Export "Results" and "Amplification Data" slots as .xls files (video 1, minute 1:50).
4- Process the table using free software LinRegPCR (following the instructions in the videos)

5- If there is a mistake in a well and you want to eliminate it from the analysis, erase it from the slot "Amplification Data\_compact", column "N0"



\## For one replicate - - - - - - - 



INSTRUCTIONS:

"Sample Name" in "Results" slot must reflect cell line and treatment separated by a space. Leave it blank if it is a control (no sample)



example:

WT day0

WT day2

KO day0

KO day2



IPUT:
The script takes as input the amplification data (.xls file) of only one biological replicate (n=1) and the name of 1 or 2 housekeeping genes.

The mean of the technical replicates is calculated for all genes.

If 2 housekeeping are used, the geometrical mean of both is calculated.
Data is normalized with the housekeeping genes provided



OUTPUTS:
The script provides one table of normalized gene expresión data and a plot of normalized expresión for visualization. 



IMPORTANT NOTE: This is done for 1 replicate of a 3-replicate experiment. The plot is meant to ve a visual help for you to get a broad idea how the experiment went. Assuming you calculate the values for all 3 biological replicates and include the statistical analysis, plotting normalized expresión of genes is perfectly reasonable and correct! However, many scientist struggle to quickly interpret them. For that reason, I suggest once you have all 3 replicates you make a final plot with the second script that calculates and plots the mean gene expression relative to a control, including the statistics.

