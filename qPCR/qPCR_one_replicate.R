



#---- SETTING HELP MESSAGE AND OPTIONS -----

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

suppressPackageStartupMessages({
  library("optparse")
})

option_list <- list(make_option("--file", type = "character", metavar = "FILE",
                                help = "Path to qPCR .xls file"),
                    make_option("--file2", type = "character", default = NA, metavar = "FILE",
                                help = "Path to second qPCR .xls file, if needed"),
                    make_option("--housekeeping", type = "character", metavar = "STR",
                                help = "Housekeeping gene to use for normalization (e.g. RPL7)"),
                    make_option("--housekeeping2", type = "character", default = NA, metavar = "STR",
                                help = "Second housekeeping gene to use for normalization (e.g. GAPDH)"),
                    make_option("--timeseries", type = "character", default = "no", metavar = "yes/no",
                                help = "Connect time points with a line (default: %default)")
)


parser <- OptionParser(
  usage = "Rscript qPCR_one_replicate.R --file FILE --housekeeping STR [options] \n",
  description = "Description: \n Calculates normalized expression (log10 scale) from LinRegPCR output. 
  Calculates averages between experimental replicates. 
  Normalices and plots expression normalized to 1 or 2 housekeeping genes.",
  option_list = option_list
)
opt <- parse_args(parser)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$file) || is.null(opt$housekeeping)) {
  print_help(parser)
  quit(status = 1)
}


#---- SETTING LIBRARIES AND OPTIONS -----

cran_packages <- c("dplyr", "tidyr", "ggplot2", "readxl")

for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("readxl")
})

options(scipen = 999)
options(digits=5)


#---- FUNCTIONS ----
`%notin%` = Negate(`%in%`)

importar_tabla <- function(path){
  archivo <- read_xls(path, sheet = "Amplification Data_compact", skip = 3) %>%
    select(name, Amplicon, Cq, N0) 
  archivo[,-c(1,2,3)][archivo[, -c(1,2,3)] < 0] <- 0 # change negative amplification values for 0
  
  muestras <- read_xls(path, sheet = "Results", skip = 7) %>%
    select(Well, `Target Name`, `Sample Name`) %>%
    rename(Sample = `Sample Name`) %>% 
    rename(Target = `Target Name`) %>% 
    unite(name, Well:Target, sep = "_", remove = FALSE)
  muestras$Well <- factor(muestras$Well, levels=c(muestras$Well))
  muestras <-  filter(muestras, muestras$name %in% archivo$name)
  
  misdatos <- merge(archivo, muestras, by = 'name')
  misdatos$Sample[is.na(misdatos$Sample)] <- "B"
  return(misdatos)
}

promedios <- function(dframe){
  dframe <- dframe %>% group_by(Amplicon, Sample) %>%
    filter(Sample != "B")  %>%
    mutate(PROMEDIO = mean(N0, na.rm = TRUE))
  return(dframe)
}

media_geometrica <- function(dframe, housekeepings){
  geomean <- filter(dframe, Amplicon %in% c(housekeepings)) %>% 
    filter(Sample != "B") 
  if(length(housekeepings) == 2){
    mediaGeom <- geomean %>% group_by(Sample) %>%
      summarise(GEO_MEAN = exp(mean(log(PROMEDIO))))  
  }else if (length(housekeepings) == 1){
    mediaGeom <- geomean %>% ungroup() %>% mutate(GEO_MEAN = PROMEDIO) %>%
      select(Sample, GEO_MEAN)  
  }
  return(mediaGeom)
}

expresion_normalizada <- function(dframe){
  expresion <- dframe %>% mutate(EXP_NORM = PROMEDIO / GEO_MEAN) %>%
    filter(row_number() %% 2 == 1) %>%
    return(expresion)
}

analiza_mi_real <- function(camino1, camino2 = NA, hk1, hk2 = NA){
  if(is.na(hk2)){
    housekeepings <- c(hk1)
  } else {
    housekeepings <- c(hk1, hk2)
  }
  if(is.na(camino2)){
    primero <- importar_tabla(camino1)
    primero <- with(primero, primero[order(Well),])
  } else {
    tabla1 <- importar_tabla(camino1)
    tabla2 <- importar_tabla(camino2)
    primero <- rbind(with(tabla1, tabla1[order(Well),]), with(tabla2, tabla2[order(Well),]))
  }
  segundo <- promedios(primero)
  tercero <- media_geometrica(segundo, housekeepings)
  cuarto <- full_join(segundo, tercero)
  quinto <- expresion_normalizada(cuarto)
  quinto <- quinto %>% filter(Amplicon %notin% housekeepings )
  return(quinto)
}

#---- SETTING THEME ----

theme_set(theme_bw()+
            theme(plot.title = element_text(size = 18),
                  plot.subtitle = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 12, angle = 40, vjust = 1, hjust=1),
                  panel.border = element_rect(colour = "black", fill = NA),
                  strip.background=element_blank(),
                  strip.text = element_text(size=11, face="bold")))


#---- UNCOMMENT TO RUN LOCALLY ON RStudio ----

# Uncomment to run from RStudio
#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)

#path <- "path to the file"
#path2 <- "if your samples are across two runs, path to second file"
#hk1 <- "name of housekeeping gene" #do not use spaces for gene names
#hk2 <- "name of second housekeeping gene" #do not use spaces for gene names
#data <- analiza_mi_real(camino1 = path, camino2 = path2, hk1= hk1, hk2= hk2)

#---- CALCULATE ----

data <- analiza_mi_real(camino1 = opt$file, camino2 = opt$file2, hk1= opt$housekeeping, hk2= opt$housekeeping2)

data <- data %>% separate(Sample, into = c("Cell line", "Treatment"), sep = " ", remove = F) %>%
  select(-name)
write.table(data, "table_normalized_expression.txt",sep = "\t", quote = F, row.names = F)
data$Treatment <- factor(data$Treatment,
  levels = unique(data$Treatment)[order(readr::parse_number(as.character(unique(data$Treatment))))]
)


# Number of subplots to define width of images

genenum <- length(unique(data$Amplicon))
condit <- length(unique(data$Treatment))
linenum <- length(unique(data$`Cell line`))

if (genenum <=3 ){
  he = 3
} else if (genenum <=6){
  he = 6
} else if (genenum <=9){
  he = 9
}

if (genenum >= 3){
  wi = (9 + condit/5)
} else if (genenum == 2){
  wi = (6 + condit/5)
} else {
  wi = (3 + condit/5)
}

#Define colors for custom palette

my_pal <- colorRampPalette(c("gray35","lightslateblue","darkorchid1"))(linenum+1)


#---- mAKING THE PLOT ----


plot <- ggplot(data, aes(x=Treatment,y=EXP_NORM, fill=`Cell line`, color=`Cell line`, group=`Cell line`))+
  facet_wrap(~Amplicon, scales = "free_y")+
  geom_point(shape=21, size = 4, alpha=0.7,position = position_dodge2(width = 0.4))+
  scale_color_manual(values = my_pal[1:linenum])+
  scale_fill_manual(values = my_pal[1:linenum])+
  scale_y_log10()+
  ylab(paste0("Normalized expression \n(log scale)"))

if (opt$timeseries == "yes") {
  plot <- plot +
    geom_line(linetype = "dashed", linewidth = 0.7, alpha = 0.6,
              position = position_dodge2(width = 0.4))
}

pdf("plot_normalized_expression.pdf", height = he, width = wi)
print(plot)
dev.off()
