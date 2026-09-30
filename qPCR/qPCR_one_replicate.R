
#----- Setting libraries and options -----

cran_packages <- c("optparse", "dplyr", "tidyr", "ggplot2", "readxl")

for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

suppressPackageStartupMessages({
  library("optparse")
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  library("readxl")
})

options(scipen = 999)
options(digits=5)


# - - - - 
option_list <- list(make_option(c("--file"), action = "store", type = "character"),
                    make_option(c("--file2"), action = "store", type = "character", default = NA),
                    make_option(c("--housekeeping"), action = "store", type = "character"),
                    make_option(c("--housekeeping2"), action = "store", type = "character", default = NA)
)

opt <- parse_args(OptionParser(option_list = option_list))


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


analiza_mis_reals <- function(camino1, camino2 = NA, hk1, hk2 = NA){
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

#------- Theme set ----------

theme_set(theme_bw()+
            theme(plot.title = element_text(size = 18),
                  plot.subtitle = element_text(size = 14),
                  axis.title.y = element_text(size = 14),
                  axis.text.y = element_text(size = 12),
                  axis.title.x = element_text(size = 15),
                  axis.text.x = element_text(size = 14),
                  panel.border = element_rect(colour = "black", fill = NA),
                  strip.background=element_blank(),
                  strip.text = element_text(size=11, face="bold")))

my_pal <- colorRampPalette(c('#5e4fa2ff','#358abcff','#5db2acff',
                             '#7ecaa5ff','#ddf19aff','#fff9b6ff',
                             '#fca65dff','#ed6346ff','#9e0142ff'))(9)

# - - - - Code to run from RStudio  - - - - -----
# Uncomment to run from RStudio
#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)

#path <- "path to the file"
#path2 <- "if your samples are across two runs, path to second file"
#hk1 <- "name of housekeeping gene" #do not use spaces for gene names
#hk2 <- "name of second housekeeping gene" #do not use spaces for gene names
#data <- analiza_mis_reals(camino1 = path, camino2 = path2, hk1= hk1, hk2= hk2)


data <- analiza_mis_reals(camino1 = opt$file, camino2 = opt$file2, hk1= opt$housekeeping, hk2= opt$housekeeping2)

data <- data %>% separate(Sample, into = c("Cell_Line", "Treatment"), sep = " ", remove = F)
write.csv(data, "Normalized_expression_table.cvs", quote = F, row.names = F)
data$Treatment <- factor(data$Treatment,
  levels = unique(data$Treatment)[order(readr::parse_number(as.character(unique(data$Treatment))))]
)

# Number of subplots to define width of images
genenum <- length(unique(data$Amplicon))
linenum <- length(unique(data$Cell_Line))

if (genenum <=3 ){
  he = 4
} else if (genenum <=6){
  he = 8
} else if (genenum <=9){
  he = 12
}

if (genenum >= 3){
  wi = 12
} else if (genenum == 2){
  wi = 8
} else {
  wi = 4
}


pdf("Normalized_expression_plot.pdf", height = he, width = wi)
ggplot(data, aes(x=Treatment,y=EXP_NORM, fill=Cell_Line, color=Cell_Line))+
  facet_wrap(~Amplicon, scales = "free_y")+
  geom_point(shape=21, size = 4, alpha=0.7)+
  scale_color_manual(values = my_pal)+
  scale_fill_manual(values = my_pal)+
  scale_y_log10()
dev.off()

