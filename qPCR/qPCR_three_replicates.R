
#----- INSTRUCCIONES: -----

work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(work.space)

library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
options(scipen = 999)
options(digits=5)

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
  expresion <- dframe %>% filter(Amplicon %notin% c(housekeepings)) %>%
    mutate(EXP_NORM = PROMEDIO / GEO_MEAN) 
    return(expresion)
}

analiza_mis_reals <- function(metafile, hk1, hk2 = NA){
  if(is.na(hk2)){
    housekeepings <- c(hk1)
  } else {
    housekeepings <- c(hk1, hk2)
  }
  
  primero <- importar_metadata(metafile)
  segundo <- promedios(primero)
  tercero <- media_geometrica(segundo, housekeepings)
  cuarto <- merge(segundo, tercero)
  quinto <- expresion_normalizada(cuarto, housekeepings)
  return(quinto)
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

#-------


meta <- "qPCR_file_metadata.txt"
test1 <- importar_metadata(meta)
test2 <- promedios(test1)
test3 <- media_geometrica(test2, c("HPRT1", "RPL7"))
test4 <- merge(test2, test3)
test5 <- expresion_normalizada(test4, c("HPRT1", "RPL7"))
test6 <- analiza_mis_reals(meta, "HPRT1", "RPL7")


test6 <- test6 %>% filter(Amplicon != "LINC881e12")
test7 <- summarySE(test6, "EXP_NORM", 
                   groupvars = c("sample", "condition", "Amplicon"), 
                   na.rm=TRUE) %>% filter(N >= 3)
test7$condition <- factor(test7$condition, levels= c(unique(test7$condition)))

ggplot(test7, aes(x=condition, y=EXP_NORM, fill=sample))+
  geom_point(shape=21, position = position_jitter(width=0.2))+
  facet_wrap(~Amplicon, scale = "free_y")+
  scale_y_log10()
