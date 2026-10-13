# qPCR and statistical analysis

This script was developed for the analysis of Step One Plus results, and following the concept of Ariel Waisman [video tutorials](https://www.youtube.com/watch?v=hPO7ptWOT7M). It is a 4 part series on RT-qPCR data analysis I highly recommend.





## General instructions



1- Do not use underscores for your gene names.


2- Keep names of genes and samples consistent between different runs.


3- The name "Samples" reffer to a cell line, such as the Wild Type line, or the control line.


4- Export "Results" and "Amplification Data" slots as .xls files (video 1, minute 1:50).


5- Process the table using free software LinRegPCR (following the instructions in the videos)


6- If there is a mistake in a well and you want to eliminate it from the analysis, erase it from the slot "Amplification Data\_compact", column "N0"


7- "Sample Name" in "Results" slot must reflect cell line and treatment separated by a space. Leave it blank if it is a qPCR control (no sample).
example:

|         |
|---------|
| WT day0 |
| WT day2 |
| KO day0 |
| KO day2 |





## For one replicate - - - - - - -



INPUT:
The script takes as input the amplification data (.xls file) of one or two files, only for one biological replicate (n=1) and the name of 1 or 2 housekeeping genes.

Remember, it is not good practice for a given sample, to split the measurement of a gene from its corresponding housekeeping in two different runs (different multiwell plates). But it can happen when measuring a large number of samples+genes that you need such a design. If you must do so, remember to run a control: the exact same sample and set of primers in both plates, so you can compare the measurements and make sure they are consistent.


Inside the script:

First, the mean of the technical replicates is calculated for all genes.
If 2 housekeeping genes are used, the geometrical mean of both is calculated.
Finally, data is normalized with the housekeeping gene/s provided.

Arguments needed:
```bash
Rscript <script name> --file <name_file1.xls> --file2 <name_file2.xls> --housekeeping <GENE1> --housekeeping2 <GENE2> --timeseries [yes/no]
```

Example run:
```
Rscript qPCR_one_replicate.R --file 2107_N1_real1.xls --file2 2107_N1_real2.xls --housekeeping RPL7 --housekeeping2 HPRT1 --timeseries yes
```

OUTPUTS:
The script provides one table of normalized gene expresión data and a plot of normalized expresión for visualization, in log scale.

If you select timeseries to be "yes" a dashed line Will be added to your plot, connecting the dots to make visualization easier.

IMPORTANT NOTE: This is done for 1 replicate of a 3-replicate experiment. The plot is meant to be a visual help for you to get a general idea how the experiment went. Assuming you calculate the values for all 3 biological replicates and include the statistical analysis, plotting normalized expresión of genes is perfectly reasonable and correct! However, many scientist struggle to quickly interpret them. For that reason, I suggest once you have all 3 replicates you make a final plot. You can do so with the second script that calculates and plots the mean gene expression relative to a control, including the statistical values.




## For three replicates - - - - - - -


INSTRUCTIONS:
The script takes as inputs two files and two arguments as follows:

Arguments needed:
```bash
Rscript <script name> --files <file_metadata.xls> --experiment <experiment_metadata.xls> --relative [control/cell_line] --timeseries [yes/no] --include_all [yes/no]
```

**file 1:** qPCR_file_data.txt

This file must be tab delimited and contain the path to the files and the replcate it belongs to:

|  File_name | Replicate |
|------------|-----------|
| file1.xls  |     1     |
| file2.xls  |     1     |
| file3.xls  |     2     |
| file4.xls  |     2     |


**file 2:** qPCR_experimental_data.txt

This file must be tab delimited and indicate the folowwing fields: cell lines (samples), conditions (time, concentration, treatment, etc.), genes and housekeepings used in the experiment. Example:

|    TYPE       |    VALUE  |
|---------------|-----------|
| CellLine      |     WT    |
| CellLine      |     KO    |
| Conditions    |     T0    | 
| Conditions    |     T14   |
| Genes         |    GENE1  |
| Genes         |    GENE2  |
| Genes         |    GENE3  |
| Genes         |    GENE4  |
| Housekeeping  |    HPRT1  |
| Housekeeping  |    RPL7   |


The script will only keep the information you provide in this table. Any extra genes or cell lines (simples) present in the qPCR files will be omitted. Importantly, only samples that have 3 replicates or more will be calculated, the rest will be dismissed. You can change this inside the function "calcular_resumen", by changing the parameter filter(N >= 3) to 2, if you have 2 replicates for one gene you want to include (although this is not recommended for statistical purposes).
It is important the order you choose for cell lines and conditions, since this Will be the order for alculations and plotting. 


**argument 1** --relative [control/cell_line]:

Expression will be calculated relative to the first condition annotated in the experiment metadata file, assuming this is an untreated sample or a time 0 sample.
You must choose between calculating the expresión of all lines relative either to the control cell line (the first in the experimental metadata file) or keep expresión relative whithin each cell line to the first condition. The graph included bellow as example, is relativized to the day 0 (T0) of the control cell line (WT) so that it is clear how the mean expresión of all 3 replicates relate to the behaviour of the wild type cell line.


**argument 2** --timeseries [yes/no]:

When plotting a time series, it is Good practice to include a line connecting the measurements to make the Reading of values easier.

**argument 3** --include_all [yes/no]

Sometimes it helps to plot each individual replicate along the mean of the measurements to prove consistency in your experiments. You can choose to include or not all measurements.


OUTPUT:
The script generates two tables and one plot (for now).

1 "table_relative_expression.txt": this will be used to perform statistical analysis.

2 "table_relative_summary.txt": with the calculated mean of each gene, cell line and condition, plus the standard error of the mean, standard deviation and confidence interval.

3 "plot_relative_expression.pdf": plot of the mean mentioned above plus an errorbar of the standard error of the mean.



