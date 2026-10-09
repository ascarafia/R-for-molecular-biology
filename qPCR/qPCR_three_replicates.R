
#----- INSTRUCCIONES: -----

work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(work.space)

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
options(scipen = 999)
options(digits=5)


exp_data <- "qPCR_experimental_data.txt"
file_data <- "qPCR_file_data.txt"

#---- FUNCIONES ----
`%notin%` = Negate(`%in%`)

importar_archivo <- function(path){
  archivo <- read_xls(path, sheet = "Amplification Data_compact", skip = 3) %>%
    select(name, Amplicon, Cq, N0) 
  archivo[,-c(1,2,3)][archivo[, -c(1,2,3)] < 0] <- 0 # change negative amplification values for 0
  
  muestras <- read_xls(path, sheet = "Results", skip = 7) %>%
    select(Well, `Target Name`, `Sample Name`) %>%
    rename(sample = `Sample Name`) %>% 
    rename(target = `Target Name`) %>% 
    unite(name, Well:target, sep = "_", remove = FALSE)
  muestras$Well <- factor(muestras$Well, levels=c(muestras$Well))
  muestras <-  filter(muestras, muestras$name %in% archivo$name)
  
  misdatos <- merge(archivo, muestras, by = 'name')
  misdatos$sample[is.na(misdatos$sample)] <- "B"
  
  return(misdatos)
}

importar_metadata <- function(myfile){
  metadata <- read.delim(myfile)
  files <- metadata[,1]
  replicate <- metadata[,2]
  tabla <- data.frame()
  for(i in seq_along(files)){
    archivo <- importar_archivo(files[i]) %>%
      mutate(rep = replicate[i])
    tabla <- rbind(tabla, archivo)
  }
  return(tabla)
}

promedios <- function(dframe){
  dframe <- dframe %>% group_by(Amplicon, sample, rep) %>%
    filter(sample != "B")  %>% group_by(Amplicon, sample, rep) %>%
    summarise(PROMEDIO = mean(N0, na.rm = TRUE), .groups = "drop_last") %>%
    separate(sample, into = c("sample", "condition"), sep = " ")
  return(dframe)
}

media_geometrica <- function(dframe, housekeepings){
  geomean <- filter(dframe, Amplicon %in% c(housekeepings)) 
  if(length(housekeepings) == 2){
    mediaGeom <- geomean %>% group_by(sample, condition, rep) %>%
      summarise(GEO_MEAN = exp(mean(log(PROMEDIO))), .groups = "drop_last")  
  }else if (length(housekeepings) == 1){
    mediaGeom <- geomean %>% group_by(sample, condition, rep) %>%
      summarise(GEO_MEAN = mean(PROMEDIO), .groups = "drop_last")
  }else{
    print("Select only 1 or 2 Housekeeping genes!")
  }
  return(mediaGeom)
}

expresion_normalizada <- function(dframe, housekeepings){
  expresion <- dframe %>% #filter(Amplicon %notin% c(housekeepings)) %>%
    mutate(EXP_NORM = PROMEDIO / GEO_MEAN) 
    return(expresion)
}

analiza_mis_reals <- function(metafile, expfile){
  experimento <- read.delim(expfile)
  experimento <- split(experimento$VALUE, 
                       factor(experimento$TYPE, levels = unique(experimento$TYPE)))
  cell_lines <- experimento$CellLine
  conditions <- experimento$Conditions
  genes <- experimento$Genes
  housekeepings <- experimento$Housekeeping
  
  primero <- importar_metadata(metafile)
  segundo <- promedios(primero) 
  tercero <- segundo %>%
    filter(sample %in% cell_lines) %>% 
    filter(condition %in% conditions)
  cuarto <- media_geometrica(tercero, housekeepings)
  quinto <- merge(tercero, cuarto) %>%
    filter(Amplicon %in% genes)
  sexto <- expresion_normalizada(quinto, housekeepings)
  
  sexto$sample <- factor(sexto$sample, levels = c(cell_lines))
  sexto$condition <- factor(sexto$condition, levels = c(conditions))
  sexto$Amplicon <- factor(sexto$Amplicon, levels = c(genes))
  
  return(sexto)
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

relati_log <- function(dframe){
  dframe$EXP_NORM[dframe$EXP_NORM == 0] <- 0.00001 #for logs sake
  dframe <- dframe %>% mutate(LOG_NORM = log2(EXP_NORM))
  logarithm <- dframe %>% group_by(Amplicon, sample, condition) %>% mutate(AVG_LOG = mean(LOG_NORM))
  logarithm <- logarithm %>% group_by(Amplicon) %>% mutate(REL_EXP = LOG_NORM - first(AVG_LOG))
  logarithm <- logarithm %>% group_by(Amplicon, sample) %>% mutate(FOLD_CH = LOG_NORM - first(AVG_LOG))
  return(logarithm)
}

calculo_resumen <- function(dframe, relative){
  logdata <- relati_log(dframe)
  if (relative == "control"){
    final <- summarySE(logdata, measurevar = "REL_EXP", 
             groupvars = c("Amplicon", "sample", "condition"),na.rm=TRUE) %>% 
             filter(N >= 3)
  } else if (relative == "cell_line"){
    final <- summarySE(logdata, measurevar = "FOLD_CH", 
             groupvars = c("Amplicon", "sample", "condition"),na.rm=TRUE) %>% 
             filter(N >= 3) %>% dplyr::rename(REL_EXP = FOLD_CH)
  } else {
    print("Choose relativization between: cell_line or control")
  }
  detach("package:plyr", unload = TRUE)
  return(final)
}
  
#-------

data <- analiza_mis_reals(file_data, exp_data)
resumen <- calculo_resumen(data, "control")


ggplot(resumen, aes(x=condition, y=REL_EXP, fill=sample, color=sample)) +
  geom_point(size=4,position = position_dodge2(width = 0.4), shape = 21) +
  facet_wrap(~ Amplicon, scales = "free_y")+
  scale_fill_manual(values = c("gray35",  "#66BD63", "#1A9850"))+
  scale_color_manual(values = c("black",  "black", "black"))





