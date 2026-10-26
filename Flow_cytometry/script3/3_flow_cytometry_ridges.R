#---- SETTING HELP MESSAGE AND OPTIONS -----

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

suppressPackageStartupMessages({
  library("optparse")
})

option_list <- list(make_option("--file", type = "character", metavar = "FILE",
                                help = "Metadata with files to be imported"),
                    make_option("--channel", type = "character", metavar = "STR",
                                help = "Name of the channel to be plotted"),
                    make_option("--threshold", type = "numeric", default = 10000, metavar = "NUM",
                                help = "Threshold to separate positive from negative cells (default: %default)"),
                    make_option("--palette", type = "character", default = "viridis", metavar = "warm/viridis",
                                help = "Color palette to use for the plot (default: %default)")
)


parser <- OptionParser(
  usage = "\n \n Rscript 3_flow_cytometry_ridges.R --file FILE --channel CHANNEL --threshold NUMBER --palette [warm/viridis] \n",
  description = "Description: 
  The script takes as input cytometry data tables (.csv) that have been previously filtered and compensated,
  and generates a ridgeline plot. Percentages are calculated based on the cutoff value (threshold) selected with the 
  \"threshold\" option provided, for the selected channel.",
  option_list = option_list
)

opt <- parse_args(parser)

if (is.null(opt$file) || is.null(opt$channel)) {
  print_help(parser)
  quit(status = 1)
}


#---- SETTING LIBARIES AND OPTIONS ----

cran_packages <- c("dplyr","dplyr","ggridges","scales","viridis","forcats")

for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

suppressPackageStartupMessages({
  library("dplyr")
  library("ggplot2")
  library("ggridges")
  library("scales")
  library("viridis")
  library("forcats")
  
})

options(scipen = 999)
options(stringsAsFactors = FALSE)


#---- FUNCTIONS ----

import_tables <- function(metadata_file, channel){
  metadata  <- read.delim(metadata_file, header = TRUE)
  file_list <- metadata[, 1]
  condition <- metadata[, 2]
  tidy_data <- data.frame()
  for (i in seq_along(file_list)) {
    data <- read.csv(file_list[i]) %>%
      select(all_of(channel)) %>%
      mutate(treatment = condition[i])
    tidy_data <- rbind(tidy_data, data)
  }
  tidy_data <- tidy_data %>%
    filter(.data[[channel]] != 0, .data[[channel]] <= 100000)
  tidy_data$treatment <- factor(tidy_data$treatment, levels = unique(tidy_data$treatment))
  return(tidy_data)
}

calculate_percent <- function(data, channel, thresh1 = 10000) {
  data %>%
    group_by(treatment) %>%
    dplyr::summarise(percentage = round(mean(.data[[channel]] > thresh1) * 100, 2)) %>%
    mutate(negative = 100 - percentage)
}


#---- CODE TO RUN FROM R STUDIO ----
# Uncomment to run from R Studio

#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)

#metafile <- path to metadata file
#channel = quoted name of channel
#thresh1 = number to use as threshold

#data <- import_tables(metafile, channel = channel) 
#data_percentage <- calculate_percent(data, channel = channel, thresh1 = thresh1)


#---- PAREPARE THE DATA ----

data <- import_tables(opt$file, channel = opt$channel) 
data_percentage <- calculate_percent(data, channel = opt$channel, thresh1 = opt$threshold)
write.table(data_percentage, "percentages.txt", quote = F, row.names = F)

#---- PLOT SETTING ----

theme_set(theme_bw()+
            theme(plot.title = element_text(size = 18),
                  plot.subtitle = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 12, vjust = 1, hjust=1),
                  panel.border = element_rect(colour = "black", fill = NA),
                  strip.background=element_blank(),
                  strip.text = element_text(size=11, face="bold")))


my_pal <- colorRampPalette(c('gray35','#5e4fa2ff','#4d71b2ff','#358abcff',
                             '#5db2acff','#7ecaa5ff','#ddf19aff','#fff9b6ff',
                             '#fca65dff','#ed6346ff','#9e0142ff'))(12)

my_viridis <- viridis(length(unique(data$treatment)))
my_viridis[1] <- "gray35"

treat <- length(unique(data$Treatment))

if (treat <=3 ){
  he = 4
} else if (treat <=6){
  he = 7
} else {
  he = 10
}


#---- MAKING THE PLOT ----

plot <- ggplot(data, aes(x = FL2.A, y = fct_rev(treatment), fill = treatment, group = treatment))+
  geom_density_ridges(alpha = 0.5)+
  geom_vline(xintercept = 11000,linetype = "dashed", color = "gray35")+
  scale_x_continuous(limits = c(800, 100000), trans =log10_trans(),
                     labels = trans_format("log10", math_format(10^.x))) +
  geom_text(data = data_percentage, aes(x = 80000, hjust = 0.5, vjust = -1,
                label = paste0(percentage, "%")), size = 4)+
  geom_text(data = data_percentage, aes(x = 800, hjust = 0, vjust = -1,
                label = paste0(negative, "%")), size = 4)+
  ylab("") +
  xlab("Fluorescence (GFP)")


if (opt$palette == "warm"){
  plot <- plot + scale_fill_manual(values = c(my_pal[1], my_pal[5:14]))
} else if (opt$palette == "viridis"){
  plot <- plot + scale_fill_manual(values = my_viridis)
}

pdf("ridgeline_plot.pdf", height = he, width = 6.5)
print(plot)
dev.off()

