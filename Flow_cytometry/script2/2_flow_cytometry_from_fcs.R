##############################################
#---- Adding help message ---- 
# Ignore this first section to run from R Studio

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

suppressPackageStartupMessages({
  library("optparse")
})


option_list <- list(make_option("--meta_file", type = "character", metavar = "FILE", 
                                help = "Metadata table: sample and psth to the .csv file."),
                    make_option("--compensate", type = "character", default = "no", metavar = "yes/no", 
                                help = "Perform compensation (default: %default)"),
                    make_option("--quality", type = "character", default = "no", metavar = "yes/no", 
                                help = "Perform quality filtering (default: %default)"),
                    make_option("--channel1", type = "character", metavar = "STR", 
                                help = "Name of channel 1, with . instead of - (e.g. FL1.A)"),
                    make_option("--channel2", type = "character", metavar = "STR", 
                                help = "Name of channel 2, with . instead of - (e.g. FL2.A)")
)

parser <- OptionParser(
  usage = "Rscript 1_flow_cytometry_from_csv.R --meta_file FILE --channel1 CHANNEL --channel2 CHANNEL --gene1 NAME --gene2 NAME \n",
  description = "Description: 
  Makes a plot based on the .csv files provided (1 per sample).
  The .csv files are exported from cytometry software after manually gating, compensating and exporting.",
  option_list = option_list
)
opt <- parse_args(parser)

#---- Setting libraries ----
# List of packages that will be used
cran_packages <- c("dplyr", "tidyr", "ggplot2", "optparse")
bioc_packages <- c("flowCore","openCyto","ggcyto","flowStats","flowAI")

# Make sure BiocManager is installed, install bioc packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}

# Install missing Bioconductor packages
for (package in bioc_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package, update = FALSE, ask = FALSE, quiet = TRUE)
  }
}

for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# Load packages
suppressPackageStartupMessages({
  library("ggplot2")
  library("dplyr")
  library("tidyr")
  library("optparse")
  library("flowCore")
  library("openCyto")
  library("ggcyto")
  library("flowStats")
  library("flowAI")
})

options(scipen = 999)
options(stringsAsFactors = FALSE)


#---- Setting functions ----

import_fcsfiles <- function(myfile){
  metadata <- read.delim(myfile)
  samples <- metadata[,1]
  cytofiles <- metadata[,2]
  fcsset <- flowCore::read.flowSet(cytofiles, path=".")
  pData(fcsset)$name <- c(samples)
  return(fcsset)
}


compensation <- function(cytoset, afirmation){
  if (afirmation == "no"){
    fcs_comp = cytoset
  } else if (afirmation == "yes"){
    spill_list = spillover(cytoset[[1]]) 
    comp_matrix = spill_list[lengths(spill_list) != 0][[1]]
    names_list <- data.frame(names = colnames(cytoset), 
                             num = seq(1:length(colnames(cytoset))))
    current_indices <- as.numeric(colnames(comp_matrix))
    new_names <- names_list$names[match(current_indices, names_list$num)]
    colnames(comp_matrix) <- new_names
    rownames(comp_matrix) <- new_names
    fcs_comp <- compensate(cytoset, comp_matrix)
  } else {
    print("compensation parameter must be lowercase yes/no")
  }
  return(fcs_comp)
}


quality_filter <- function(cytoset, afirmation){
  if (afirmation == "no"){
    fcsset_qc = cytoset
  } else if (afirmation == "yes"){
    fcsset_qc = flow_auto_qc(cytoset)
  } else {
    print("quality-control parameter must be lowercase yes/no")
  }
  return(fcsset_qc)
}



#---- Uncomment this part to run from R studio ----

#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)

#meta_file <- "example_metadata_fcs.txt"
#compensate <- "yes"
#quality <- "yes"
#channel1 <- "FL1.A"
#channel2 <- "FL2.A"

#---- IMPORT SAMPLES TO CYTOSET AND PROCESS ----

# Import fcs files as cytoset
fcsset <- import_fcsfiles(meta_file = opt$meta_file)

# Compensate if needed
fcsset_comp <- compensation(fcsset, compensate = opt$compensate)

# Check quality control and filter
fcsset_qc <- quality_filter(fcsset_comp, quality = opt$quality)


#---- GATES AND VISUALIZATION PLOTS ----

# create gating object
gatingset <- GatingSet(fcsset_qc)

# Number of subplots to define width of images
length(pData(fcsset)$name)
number <- sum(4+3*(length(pData(fcsset)$name)-1))

# You might need to adjust the gatings to meet your own experiment
# The workflow is: 1 define a gate, 2 visualize in plot, 3 add to the gatingset object

# GATE 1
# gate 1: filter the cells and remove debris
cells <- polygonGate(filterId = "cells",
                     "FSC-A"=c(27e5,5e6,1e7,125e5,13e6,1e7,5e6,4e6),
                     "SSC-A"=c(1e5,0,3e5,1e6,2e6,25e5,15e5,9e5))

# Visualize gate 1 in the graph to see if it matches the desired population
pdf("gate1_cells.pdf", height =4, width = number)
ggcyto(gatingset, aes(x="FSC-A",y="SSC-A"), subset="root")+
  geom_hex(bins = 250)+
  geom_gate(cells)+
  ggcyto_par_set(limits = list(y = c(0, 6000000))) +
  labs(caption = "Make sure the gate matches the cell population. \n If not, adjust cell polygonal gate (gate 1)")+
  theme_bw()
dev.off()

# Incorporate gate 1 (cells) into the set
gs_pop_add(gatingset, cells, parent = "root", name = "cells")
recompute(gatingset) 
#gs_get_pop_paths(gatingset) # uncomment to see nested population paths


# GATE 2
# gate 2: filter singlets from doublets
singlets <- polygonGate(filterId = "singlets",
                        "FSC-A"=c(4000000,9700000,8000000,3000000),
                        "FSC-H"=c(2200000,5500000,6500000,3500000)) 

# Visualize gate 2 in the graph to see if it matches the desired population
pdf("gate2_singlets.pdf", height =4, width = number)
ggcyto(gatingset, aes(x="FSC-A",y="FSC-H"), subset="cells")+
  geom_hex(bins = 250)+
  geom_gate(singlets)+
  ggcyto_par_set(limits = list(y = c(0, 6000000))) +
  labs(caption = "Make sure the gate matches the cell population. \n If not, adjust singlet polygonal gate (gate 2)")+
  theme_bw()
dev.off()

# Incorporate gate 2 (singlets) into the set
gs_pop_add(gatingset, singlets, parent = "root", name = "singlets")
recompute(gatingset) 
#gs_get_pop_paths(gatingset) # uncomment to see nested population paths


# GATE 3

# Transform the channels that will be plotted to log scale 
# m parameter is data specific. Adjust it to your own.
transformation <- estimateLogicle(gatingset[[1]], channels = c(opt$channel1, opt$channel2), m = 7)
transformed_gateset <- transform(gatingset, transformation)

# gate 3: population of interest
# Adjust values to match population in your plot
# also adjust position in geom_stats location to move the label
mp <- polygonGate(filterId = "Mesoderm Progenitor",
                  "FL1-A"=c(3.8, 4.2, 4.7, 5.1, 5.1, 4.7, 4.3),
                  "FL2-A"=c(3.6, 3.0, 3.0, 3.5, 4.3, 4.7, 4.7)) 


# Visualize gate 3 in the data
pdf("Population_of_study.pdf", height =4, width = number)
ggcyto(transformed_gateset, aes(.data[[ch1]], .data[[ch2]]), subset="singlets")+
  geom_hex(bins = 80) +
  geom_gate(mp)+
  geom_stats(location = "fixed", 
             adjust = c(5,2.5), 
             type = "percent",
             digits = 1)+
  ggcyto_par_set(limits = list(x = c(0,6), y=c(0,6))) +
  theme_bw()
dev.off()



