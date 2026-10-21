#---- SETTING HELP MESSAGE AND OPTIONS -----

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

suppressPackageStartupMessages({
  library("optparse")
})


option_list <- list(make_option("--file", type = "character", metavar = "FILE",
                                help = "Table with"),
                    make_option("--paired", type = "character", default = "no", metavar = "yes/no",
                                help = "Type of data between paired or unpaired (default: %default)"),
                    make_option("--compare_between", type = "character", metavar = "Cell_line/Treatment",
                                help = "The comparison between cell lines within a treatment or between treatments for each cell line")
)

parser <- OptionParser(
  usage = "Rscript qPCR_statistical_analysis.R --file FILE --paired [yes/no] \n",
  description = "Description: 
  The script takes gene expression values of at least 3 replicates and performs statistical analysis. 
  Takes as input the table saved from script qPCR_three_replicates.R and uses the log-transformed expression data
  to perform ANOVA or T-test based on the number of samples and type of data (paired or not).
  If there are no differences found, it prints out a message.
  
  If statistical differences are found, before giving the results, the corresponding suppositions are evaluated. 
  However, if suppositions are not met, non-parametric tests are performed instead:
 
  Kruskal-Wallis instead of one-way ANOVA (independent samples).
  Friedman instead of ANOVA (paired samples).
  Mann-Whitney U-test instead of T-test (independent samples).
  Wilcoxon instead of T-test (for paired samples).",
  option_list = option_list
)

opt <- parse_args(parser)

if (is.null(opt$file) || is.null(opt$compare_by)) {
  print_help(parser)
  quit(status = 1)
}


#---- SETTING LIBARIES AND OPTIONS -----

cran_packages <- c("stats","car","sjstats","dplyr","tidyr","ggplot2")

for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

suppressPackageStartupMessages({
  library("stats")
  library("car")
  library("sjstats")
  library("dplyr")
  library("tidyr")
  library("ggplot2")
  
})

options(scipen = 999)
options(digits=5)

#---- FUNCIONES ----
#---- UNCOMMENT TO RUN LOCALLY ----

#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)


#---- STATISTICAL ANALYSIS ------

read_in_data <- function(path){
  data <- read.table(path, header=T) %>% select(Cell_line, Treatment, Amplicon, rep, LOG_NORM) %>%
    group_by(Cell_line, Treatment, Amplicon) %>% 
    mutate(total = n()) %>%  
    dplyr::filter(total >= 3) %>%
    ungroup()
  return(data)
}

data <- read_in_data("table_relative_expression.txt")
data_nkx <- filter(data, Amplicon == "NKX2-5")

if(opt$compare_between == "Cell_line"){
  compare_between = "Cell_line"
  inside_of = "Treatment"
} else if (opt$compare_between == "Treatment") {
  compare_between = "Treatment"
  inside_of = "Cell_line"
}

# data %>% group_by(Amplicon) %>% group_modify(~ _ _ _ _ _(.x))

statistical_analysis <- function(data_by_gene){
  formula1 <- reformulate(c("rep",compare_between), response = "LOG_NORM")
  
  final_table = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(final_table) <- c("test", "p_value", "significance")
  # ANOVA
  data_by_gene$rep <- factor(data_by_gene$rep)
  model <- aov(formula1, data = data_by_gene)
  aov_result <- broom::tidy(model)
  p1 <- aov_result$p.value[aov_result$term == compare_between]
  aov_asterisk <- as.character(cut(p1, breaks = c(0, 0.001, 0.01, 0.05, 1),
                                   labels = c("***", "**", "*", "")))
  final_table <- rbind(final_table, data.frame(test = "ANOVA", p_value = p1, significance = aov_asterisk))
  # SHAPIRO (normality test)
  residuos_aov <- model$residuals
  shapi <- shapiro.test(residuos_aov)
  shapiro_result <- broom::tidy(shapi)
  p2 <- shapiro_result$p.value
  if (p2 < 0.1){
    print("WARNING: Normality of residuals was not met! Proceeding with the analysis.")
  }
  # LEVENE (homoscedasticity)
  formula2 <- reformulate(compare_between, response = "LOG_NORM")
  leve <- leveneTest(formula2, data = data_by_gene)
  leve_result <- broom::tidy(leve)
  p3 <- leve_result$p.value
  if (p3 < 0.1){
    print("WARNING: Independence of variances was not met! This data is not suitable for analysis with parametric tests.
          Proceeding with non-parametric statistical analysis.")
  }
  
  return(final_table)
}

test <- data %>% group_by(Amplicon, Treatment) %>% group_modify(~ statistical_analysis(.))


# Now testing tukey test loop for multiple comparisons!
result <- data.frame()
for (i in seq_along(unique(data_nkx$Treatment))){
  data_nkx_byday <- data_nkx[data_nkx$Treatment == unique(data_nkx$Treatment)[i],]
  tuky <- TukeyHSD(aov(LOG_NORM ~ Cell_line, data = data_nkx_byday))
  tuky_result <- broom::tidy(tuky)
  tuky_result$testing <- unique(data_nkx$Treatment)[i]
  result <- rbind(result, tuky_result)
}


