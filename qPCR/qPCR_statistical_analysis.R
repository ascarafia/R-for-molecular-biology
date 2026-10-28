#---- SETTING HELP MESSAGE AND OPTIONS ----

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
  usage = "Rscript qPCR_statistical_analysis.R --file FILE --paired [yes/no] --compare_between [Cell_line/Treatment] \n",
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

if (is.null(opt$file) || is.null(opt$compare_between)) {
  print_help(parser)
  quit(status = 1)
}


#---- SETTING LIBARIES AND OPTIONS ----

cran_packages <- c("broom","stats","car","sjstats","rstatix","dplyr","tidyr","ggplot2")

for (package in cran_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

suppressPackageStartupMessages({
  library("broom")
  library("stats")
  library("car")
  library("sjstats")
  library("rstatix")
  library("dplyr")
  library("tidyr")
  library("ggplot2")
})

options(scipen = 999)
options(digits=5)

#---- FUNCIONES ----

# read in data
read_in_data <- function(path){
  data <- read.table(path, header=T) %>% select(Cell_line, Treatment, Amplicon, rep, LOG_NORM) %>%
    group_by(Cell_line, Treatment, Amplicon) %>% 
    mutate(total = n()) %>%  
    dplyr::filter(total >= 3) %>%
    ungroup()
  data$Treatment <- factor(data$Treatment, levels = unique(data$Treatment))
  data$rep <- factor(data$rep)
  return(data)
}

if (opt$compare_between == "Cell_line") {
  compare_between <- "Cell_line"
  inside_of <- "Treatment"
} else if (opt$compare_between == "Treatment") {
  compare_between <- "Treatment"
  inside_of <- "Cell_line"
}

# Set transform p-values into asterisks
astrx <- function(p) {
  as.character(cut(p, breaks = c(0, 0.001, 0.01, 0.05, 1),
                   labels = c("***", "**", "*", ""), include.lowest = TRUE))
}

# Select formula for statistical analyses
build_formula <- function() {
  model_terms <- if (opt$paired == "yes") c("rep", compare_between) else compare_between
  reformulate(model_terms, response = "LOG_NORM")
}

# ANOVA test (when 2 samples, equal results to T-test)
anova_test <- function(model) {
  aov_result <- broom::tidy(model)
  aov_p <- aov_result$p.value[aov_result$term == compare_between]
  data.frame(test = "ANOVA", contrast = compare_between,
             p_value = aov_p, significance = astrx(aov_p))
}

# Tukey test loop for multiple comparisons!
multiple_testing_tukey <- function(model) {
  tukey_result <- broom::tidy(TukeyHSD(model, which = compare_between))
  data.frame(test = "TUKEY", contrast = tukey_result$contrast,
             p_value = tukey_result$adj.p.value,
             significance = astrx(tukey_result$adj.p.value))
}

# Non-parametric, independent samples: Mann-Whitney (2 groups) or Kruskal-Wallis + Dunn (3 or more)
kruskal_dunn_test <- function(data_by_gene) {
  formula_groups <- reformulate(compare_between, response = "LOG_NORM")
  
  if (n_distinct(data_by_gene[[compare_between]]) == 2) {
    mw_p <- wilcox.test(formula_groups, data = data_by_gene)$p.value
    return(data.frame(test = "MANN-WHITNEY", contrast = compare_between,
                      p_value = mw_p, significance = astrx(mw_p)))
  }
  
  dunn <- rstatix::dunn_test(data_by_gene, formula_groups, p.adjust.method = "holm")
  data.frame(test = "DUNN",
             contrast = paste(dunn$group2, dunn$group1, sep = "-"),  # same direction as Tukey
             p_value = dunn$p.adj, significance = astrx(dunn$p.adj))
}

# Statistical analysis  
statistical_analysis <- function(data_by_gene, key = NULL) {
  label <- paste(vapply(key, as.character, character(1)), collapse = " | ")
  
  empty_table <- data.frame(test = character(), contrast = character(),
                            p_value = numeric(), significance = character(),
                            shapiro_p = numeric(), levene_p = numeric())
  
  n_groups <- n_distinct(data_by_gene[[compare_between]])
  if (n_groups < 2) return(empty_table)   # nothing to compare in this subset
  
  formula_groups <- reformulate(compare_between, response = "LOG_NORM")
  
  # LEVENE (homogeneity of variances)
  levene_p <- suppressWarnings(leveneTest(formula_groups, data = data_by_gene))[1, "Pr(>F)"]
  shapiro_p <- NA_real_
  
  if (is.na(levene_p) || levene_p >= 0.1) {
    # ---- parametric branch ----
    model <- aov(build_formula(), data = data_by_gene)
    # SHAPIRO (normality of the residuals)
    shapiro_p <- tryCatch(shapiro.test(residuals(model))$p.value,
                          error = function(e) NA_real_)
    if (!is.na(shapiro_p) && shapiro_p < 0.1) {
      message("WARNING [", label, "]: normality of residuals not met (Shapiro p = ",
              signif(shapiro_p, 3), "). Proceeding with the analysis.")
    }
    
    if (n_groups == 2) {
      final_table <- anova_test(model)
    } else {
      final_table <- multiple_testing_tukey(model)
    }
    
  } else {
    # ---- non-parametric branch ----
    message("WARNING [", label, "]: variances not homogeneous (Levene p = ",
            signif(levene_p, 3), "). Proceeding with non-parametric tests.")
    
    if (opt$paired == "yes") {
      # Friedman / Wilcoxon signed-rank still to be written
      message("  Paired non-parametric tests are not implemented yet.")
      final_table <- data.frame(test = "paired non-parametric (not implemented)",
                                contrast = compare_between,
                                p_value = NA_real_, significance = NA_character_)
    } else {
      final_table <- kruskal_dunn_test(data_by_gene)
    }
  }
  
  final_table$shapiro_p <- shapiro_p
  final_table$levene_p  <- levene_p
  final_table
}


#---- UNCOMMENT TO RUN LOCALLY ----

#work.space <- dirname(rstudioapi::getSourceEditorContext()$path)
#setwd(work.space)

# still is process . . . 

#---- STATISTICAL ANALYSIS ----

data <- read_in_data("table_relative_expression.txt")

results <- data %>%
   group_by(across(all_of(c("Amplicon", inside_of)))) %>%
   group_modify(~ statistical_analysis(.x, .y))

write.table(results, "table_statistical_analysis.txt", quote=F)

#---- TRYING TO PLOT ----
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

num_lines <- length(unique(data$Cell_line))
my_pal <- colorRampPalette(c("gray35","lightslateblue","blueviolet","darkorchid1"))(num_lines+2)

genenum <- length(unique(data$Amplicon))
condit <- length(unique(data$Treatment))

# height of plot is calculated by number of genes to be plotted
if (genenum <=3 ){
  he = 3
} else if (genenum <=6){
  he = 6
} else if (genenum <=9){
  he = 9
}

# width is calculated by number of genes and amount of conditions
if (genenum >= 3){
  wi = 9 
} else if (genenum == 2){
  wi = 6
} else {
  wi = 3 
}

sum_data <- read.table("table_relative_summary.txt", header = T)
sum_data$Treatment <- factor(sum_data$Treatment, levels = unique(sum_data$Treatment))
sum_data$Cell_line <- factor(sum_data$Cell_line, levels = unique(sum_data$Cell_line))


plot <- ggplot(sum_data, aes(x=Treatment, y=REL_EXP, fill=Cell_line, color=Cell_line, group=Cell_line)) +
  facet_wrap(~ Amplicon, scales = "free_y")+
  geom_point(size = 4, alpha = 0.8, shape = 21,
             position = position_dodge2(width = 0.4)) +
  geom_errorbar(aes(ymin=REL_EXP-se, ymax=REL_EXP+se),
                color="black", linewidth=0.2,
                width=0.4, position = position_dodge2(width = 1))+
  geom_line(linetype = "dashed", linewidth = 0.7, alpha = 0.6,
            position = position_dodge2(width = 0.4))+
  geom_text(data=results, inherit.aes = F, aes(x = Treatment, y = 10, label = significance),
            position = position_stack(), color = 'black', size = 8)+
  ylab(paste0("Fold Change \n(relative to control)"))+
  scale_fill_manual(values = my_pal[1:num_lines])+
  scale_color_manual(values = my_pal[1:num_lines])


pdf("plot_statistics.pdf", height = he, width = wi)
print(plot)
dev.off()
