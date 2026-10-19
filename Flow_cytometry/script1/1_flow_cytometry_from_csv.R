##############################################
#---- Adding help message ---- 
# Ignore this first section to run from R Studio

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

suppressPackageStartupMessages({
  library("optparse")
})


option_list <- list(make_option("--meta_file", type = "character", metavar = "FILE", help = "Metadata table: sample and psth to the .csv file."),
                    make_option("--channel1", type = "character", metavar = "STR", help = "Name of channel 1, with . instead of - (e.g. FL1.A)"),
                    make_option("--channel2", type = "character", metavar = "STR", help = "Name of channel 2, with . instead of - (e.g. FL2.A)"),
                    make_option("--gene1", type = "character", metavar = "STR", help = "Name of gene or protein measured in channel 1 (e.g. GFP)"),
                    make_option("--gene2", type = "character", metavar = "STR", help = "Name of gene or protein measured in channel 2 (e.g. GAPDH)")
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

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(rlang)
})

options(scipen = 999)

#---- Setting functions ----

tidying_tables <- function(myfile, ch1, ch2, name1, name2){
  metadata <- read.delim(myfile)
  conditions <- metadata[,1]
  files <- metadata[,2]
  tidy_data <- data.frame()
  for(i in seq_along(files)){
    tidy_table <- read.delim(files[i], sep=",") %>%
      select(all_of(c(ch1, ch2))) %>% 
      rename(!!name1 := !!ch1, !!name2 := !!ch2) %>%
      mutate(condition = conditions[i])
    tidy_data <- rbind(tidy_data, tidy_table)
  }
  return(tidy_data)
}

classify_populations <- function(df, thresh1 = 100000, thresh2 = 100000, name1, name2) {
  df %>% mutate(population = case_when(
    df[,1] >= thresh1 & df[,2] >= thresh2 ~ paste0(!!name1,"+\n",!!name2,"+"),
    df[,1] >= thresh1 & df[,2] <  thresh2 ~ paste0(!!name1,"+"),
    df[,1] <  thresh1 & df[,2] >= thresh2 ~ paste0(!!name2,"+"),
    TRUE ~ ""
  )
  )
}

build_labels <- function(df, name1, name2, x_range = c(100, 10000000), y_range = c(100, 10000000)) {
  all_populations <- c(
    paste0(name1, "+\n", name2, "+"),
    paste0(name1, "+"),
    paste0(name2, "+"),
    ""
  )
  
  summary_df <- df %>% group_by(condition) %>%
    count(population) %>% ungroup() %>%
    complete(condition, population = all_populations, fill = list(n = 0)) %>%
    group_by(condition) %>%
    mutate(percent = paste0(round(100 * n / sum(n), 1),"%")) %>% ungroup()
  
  positions <- tibble(
    population = all_populations,
    !!name1 := c(y_range[2], y_range[2], y_range[1], y_range[1]),
    !!name2 := c(x_range[2], x_range[1], x_range[2], x_range[1]),
    hjust    = c(1, 0, 1, 0),   
    vjust    = c(1, 1, 0, 0))
  
  label_summary <- summary_df %>%
    full_join(positions, by = "population")
  
  return(label_summary)
}


#---- Code to run from within R ----
# Uncomment the code to run from R Studio

# Seting the working directory 
#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)

#meta_file <- "name_of_metadata_file.txt"
#channel1 = "name of first channel"
#channel2 = "name of second channel"
#gene1 = "name of the first gene"
#gene2 = "name of the second gene"


#---- IMPORT AND PROCESS DATA ---- 

data <- tidying_tables(myfile = opt$meta_file, 
                       ch1 = opt$channel1, 
                       ch2 = opt$channel2, 
                       name1 = opt$gene1, 
                       name2 = opt$gene2)

data <- classify_populations(data, thresh1 = 100000, thresh2 = 100000, name1 = opt$gene1, name2 = opt$gene2) 

label <- build_labels(data, name1 = opt$gene1, name2 = opt$gene2)

#---- THEME SETTING ----
theme_set(theme_bw()+
            theme(axis.title.y=element_text(size = 12),
                  axis.text.y=element_text(size = 12),
                  axis.title.x = element_text(size = 12),
                  axis.text.x = element_text(size = 12),
                  legend.text = element_text(size = 10),
                  legend.title = element_text(size = 10)))

my_pal <- colorRampPalette(c('#5e4fa2ff','#4d71b2ff','#358abcff',
                             '#5db2acff','#7ecaa5ff','#ddf19aff','#fff9b6ff','#fed27fff','#fca65dff',
                             '#f7824cff','#ed6346ff','#ca354dff','#9e0142ff'))(250)

#---- MAKING THE PLOT ----

pdf("cytometry_hexplot.pdf", height =4, width =10)

ggplot(data, aes(.data[[opt$gene1]], .data[[opt$gene2]]))+
  facet_wrap(~condition)+
  geom_hex(binwidth = 0.05, alpha = 0.9)+
  guides(fill = guide_colourbar(barwidth = 0.7, barheight = 15))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_x_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_color_gradientn(colors = my_pal) +
  scale_fill_gradientn(colors = my_pal) +
  geom_text(data = label, aes(x = .data[[opt$gene2]], y = .data[[opt$gene1]], hjust = hjust, vjust = vjust,
                              label = paste0(population, "\n", percent)), size = 3)+
  coord_cartesian(xlim=c(100,10000000),ylim=c(100,10000000))+
  geom_vline(xintercept=100000,linetype="dashed",linewidth = 0.3, color = "gray35")+
  geom_hline(yintercept=100000,linetype="dashed",linewidth = 0.3, color = "gray35")

dev.off()
