library(CrossTalkeR)
library(igraph)
#install.packages("remotes")
library(stringr)
library(tibble)
#library(pandoc)
library(EnhancedVolcano)
#install.packages("devtools")
#devtools::install_github("CostaLab/CrossTalkeR", ref="dev")
#MF 
#if (!require("BiocManager", quietly = TRUE)) {
#install.packages("BiocManager")
#}
#BiocManager::install("clusterProfiler")


#output <- ("/home/larissa/Documents/Masterarbeit/MF_results/crosstalker/no_cycling_pval005/pval0_05")
#output <- ("D:\\studium\\Masterarbeit\\MF_results\\crosstalker\\all_celltypes\\pval0_05")
output <- ("C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\MF_results\\crosstalker\\final_filtering")

#tmp <- read.csv("/home/larissa/Documents/Masterarbeit/MF_results/liana/u10_filtered/WT_lr_ready.csv")
#tmp2 <- read.csv("/home/larissa/Documents/Masterarbeit/MF_results/liana/u10_filtered/MF_lr_ready.csv")
#tmp$gene_A <- gsub('_', '+', tmp$gene_A)
#tmp$gene_B <- gsub('_', '+', tmp$gene_B)
##tmp$source <- gsub('_', '+', tmp$source)
##tmp$target <- gsub('_', '+', tmp$target)

#tmp2$gene_A <- gsub('_', '+', tmp2$gene_A)
#tmp2$gene_B <- gsub('_', '+', tmp2$gene_B)
#tmp2$source <- gsub('_', '+', tmp2$source)
#tmp2$target <- gsub('_', '+', tmp2$target)


#write.csv(tmp, "/home/larissa/Documents/Masterarbeit/MF_results/liana/u10_filtered/WT_lr_ready_plus.csv")
#write.csv(tmp2, "/home/larissa/Documents/Masterarbeit/MF_results/liana/u10_filtered/MF_lr_ready_plus.csv")

#paths <- c(
#  "WT" = "/home/larissa/Documents/Masterarbeit/MF_results/liana/final_filtering/WT_lr_ready.csv",
#  "MF" = "/home/larissa/Documents/Masterarbeit/MF_results/liana/final_filtering/MF_lr_ready.csv"
#)

#paths <- c(
#  "WT" = "D:\\studium\\Masterarbeit\\MF_results\\liana\\final_filtering\\WT_lr_ready.csv",
#  "MF" = "D:\\studium\\Masterarbeit\\MF_results\\liana\\final_filtering\\MF_lr_ready.csv"
#)

paths <- c(
  "WT" ="C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\MF_results\\liana\\final_filtering\\WT_lr_ready.csv",
  "MF" = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\MF_results\\liana\\final_filtering\\MF_lr_ready.csv"
)

data <- generate_report(paths,
                        out_path = paste0(output, "/"),
                        out_file = "MF_results.html",
                        output_fmt = "html_document",
                        report = TRUE,
                        filtered_net = TRUE,
                        comparison = list(c("MF", "WT"))
                        )

#RCC

#output <- ("/home/larissa/Documents/Masterarbeit/RCC_results/crosstalker/no_subset_glom_removed")
#output <- ("D:\\studium\\Masterarbeit\\RCC_results\\crosstalker\\no_subset_glom_removed")
output <- ("C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\crosstalker\\final_filtering")

#tmp <- read.csv("RCC_results/liana_log_norm/no_subsets/N_lr_ready.csv")
#tmp2 <- read.csv("RCC_results/liana_log_norm/no_subsets/T_lr_ready.csv")
#tmp$gene_A <- gsub('_', '+', tmp$gene_A)
#tmp$gene_B <- gsub('_', '+', tmp$gene_B)
#tmp$source <- gsub('_', '+', tmp$source)
#tmp$target <- gsub('_', '+', tmp$target)

#tmp2$gene_A <- gsub('_', '+', tmp2$gene_A)
#tmp2$gene_B <- gsub('_', '+', tmp2$gene_B)
#tmp2$source <- gsub('_', '+', tmp2$source)
#tmp2$target <- gsub('_', '+', tmp2$target)


#write.csv(tmp, "RCC_results/liana_log_norm/no_cycling_pval005/N_lr_ready_plus.csv")
#write.csv(tmp2, "RCC_results/liana_log_norm/no_cycling_pval005/T_lr_ready_plus.csv")

paths <- c(
  "N" = "RCC_results/liana_log_norm/final_filtering/N_lr_ready.csv",
  "T" = "RCC_results/liana_log_norm/final_filtering/T_lr_ready.csv"
)

data <- generate_report(paths,
  out_path = paste0(output, "/"),
  out_file = "RCC_results.html",
  output_fmt = "html_document",
  report = TRUE,
  filtered_net = TRUE,
  comparison = list(c("T", "N"))
)


##################################
#plots


#RCC <- readRDS("/home/larissa/Documents/Masterarbeit/RCC_results/crosstalker/no_subsets/LR_data_final.Rds")
RCC <- readRDS('C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\crosstalker\\final_filtering\\LR_data_final.Rds')

pdf(file= "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subsets/T_x_N_LR_scores_hist.pdf", width = 10, height = 12)
hist(RCC@tables[["T_x_N"]][["LRScore"]], breaks= "Scott")
dev.off()


pdf(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\plots\\crosstalker\\final_filtering\\CCI_T_x_N_RCC.pdf", width = 10, height = 12)
#pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/CCI_T_x_N_RCC.pdf", width = 10, height = 12)
plot_cci(RCC@graphs$T_x_N,
             colors = RCC@colors,
             coords = RCC@coords[igraph::V(RCC@graphs$T_x_N)$name, ],
             plt_name = "T_x_N",
             pg_node_size_high = 30,
             node_label_size = 0.9,
             node_label_position = 1.5,
             pg = RCC@rankings[["T_x_N"]]$Pagerank)

dev.off()



svg(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\plots\\crosstalker\\final_filtering\\CCI_T_x_N_RCC_smaller_preedit2.svg", width = 10, height = 10, , onefile = TRUE, family = "sans")
par(mar = c(5, 5, 5, 15), xpd = TRUE) 
plot_cci(RCC@graphs$T_x_N,
  plt_name = "",
  emax = NULL,
  leg = TRUE,
  low = 0,
  high = 0 / 100,
  ignore_alpha = FALSE,
  log = TRUE,
  efactor = 3,
  vfactor = 12,
  vnames = TRUE,
  pg = RCC@rankings[["T_x_N"]]$Pagerank,
  vnamescol = NULL,
  colors = RCC@colors,
  coords =  RCC@coords[igraph::V(RCC@graphs$T_x_N)$name, ],
  col_pallet = c("#89bbdb", "white", "#f74222"),
  standard_node_size = 20,
  pg_node_size_low = 10,
  pg_node_size_high = 60,
  arrow_size = 0.4,
  arrow_width = 0.7,
  node_label_position = 1.25,
  node_label_size = 0.5
)
dev.off()




pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/CCI_T_RCC.pdf", width = 10, height = 12)
new_plot_cci(RCC@graphs$T,
             colors = RCC@colors,
             coords = RCC@coords[igraph::V(RCC@graphs$T)$name, ],
             plt_name = "T",
             pg_node_size_high = 30,
             node_label_size = 0.7,
             node_label_position = 1.25,
             pg = RCC@rankings[["T"]]$Pagerank)

dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/CCI_N_RCC.pdf", width = 10, height = 12)
new_plot_cci(RCC@graphs$N,
             colors = RCC@colors,
             coords = RCC@coords[igraph::V(RCC@graphs$N)$name, ],
             plt_name = "N",
             pg_node_size_high = 30,
             node_label_size = 0.7,
             node_label_position = 1.25,
             pg = RCC@rankings[["N"]]$Pagerank)

dev.off()

pdf(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\plots\\crosstalker\\no_subset_glom_removed\\CCI_T_x_N_filtered_RCC.pdf", width = 10, height = 12)
#pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/CCI_T_x_N_filtered_RCC.pdf", width = 10, height = 12)
plot_cci(RCC@graphs$T_x_N_filtered,
             colors = RCC@colors,
             coords = RCC@coords[igraph::V(RCC@graphs$T_x_N_filtered)$name, ],
             plt_name = "T_x_N_filtered",
             pg_node_size_high = 30,
             node_label_size = 0.6,
             node_label_position = 1.4,
             pg = RCC@rankings[["T_x_N_filtered"]]$Pagerank)

dev.off()


svg(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\plots\\crosstalker\\no_subset_glom_removed\\CCI_T_x_N_RCC_filtered_smaller_preedit.svg", width = 4.5, height = 4.5, , onefile = TRUE, family = "sans")
plot_cci(RCC@graphs$T_x_N_filtered,
  plt_name = "",
  emax = NULL,
  leg = FALSE,
  low = 0,
  high = 0 / 100,
  ignore_alpha = FALSE,
  log = TRUE,
  efactor = 3,
  vfactor = 12,
  vnames = TRUE,
  pg = RCC@rankings[["T_x_N_filtered"]]$Pagerank,
  vnamescol = NULL,
  colors = RCC@colors,
  coords =  RCC@coords[igraph::V(RCC@graphs$T_x_N_filtered)$name, ],
  col_pallet = c("#306180", "white", "#f74222"),
  standard_node_size = 20,
  pg_node_size_low = 10,
  pg_node_size_high = 60,
  arrow_size = 0.4,
  arrow_width = 0.7,
  node_label_position = 1.25,
  node_label_size = 0.5
)
dev.off()


#MF <- readRDS("/home/larissa/Documents/Masterarbeit/MF_results/crosstalker/final_filtering/LR_data_final.Rds")
MF <- readRDS("C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\MF_results\\crosstalker\\final_filtering\\LR_data_final.Rds")


pdf(file= "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/all_celltypes/pval0_05/MF_x_WT_LR_scores_hist.pdf", width = 10, height = 12)
hist(MF@tables[["MF_x_WT"]][["LRScore"]],
       breaks= "Scott"
)
dev.off()

#pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/final_filtering/CCI_MF_x_WT_MF.pdf", width = 10, height = 12)
pdf(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\MF_results\\plots\\crosstalker\\final_filtering\\CCI_MF_x_WT_MF_edit.pdf", width = 10, height = 12)
plot_cci(MF@graphs$MF_x_WT,
             colors = MF@colors,
             leg = TRUE,
             coords = MF@coords[igraph::V(MF@graphs$MF_x_WT)$name, ],
             plt_name = "MF_x_WT",
             node_label_size = 0.7,
             pg_node_size_high = 30,
             pg = MF@rankings[["MF_x_WT"]]$Pagerank)

dev.off()

svg(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\MF_results\\plots\\crosstalker\\final_filtering\\CCI_MF_x_WT_MF_smaller_preedit2.svg", width = 10, height = 10, , onefile = TRUE, family = "sans")
par(mar = c(5, 5, 5, 15), xpd = TRUE) 
plot_cci(MF@graphs$MF_x_WT,
  plt_name = "",
  emax = NULL,
  leg = TRUE,
  low = 0,
  high = 0 / 100,
  ignore_alpha = FALSE,
  log = TRUE,
  efactor = 3,
  vfactor = 12,
  vnames = TRUE,
  pg = MF@rankings[["MF_x_WT"]]$Pagerank,
  vnamescol = NULL,
  colors = MF@colors,
  coords =   MF@coords[igraph::V(MF@graphs$MF_x_WT)$name, ],
  col_pallet = c("#306180", "white", "#f74222"),
  standard_node_size = 20,
  pg_node_size_low = 10,
  pg_node_size_high = 60,
  arrow_size = 0.4,
  arrow_width = 0.7,
  node_label_position = 1.25,
  node_label_size = 0.5
)
dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_0001/CCI_MF_MF.pdf", width = 10, height = 12)
new_plot_cci(MF@graphs$MF,
             colors = MF@colors,
             coords = MF@coords[igraph::V(MF@graphs$MF)$name, ],
             plt_name = "MF",
             node_label_size = 0.7,
             pg_node_size_high = 30, 
             pg = MF@rankings[["MF"]]$Pagerank)

dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_0001/CCI_WT_MF.pdf", width = 10, height = 12)
new_plot_cci(MF@graphs$WT,
             colors = MF@colors,
             coords = MF@coords[igraph::V(MF@graphs$WT)$name, ],
             plt_name = "WT",
             node_label_size = 0.7,
             pg_node_size_high = 30, 
             pg = MF@rankings[["WT"]]$Pagerank)

dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_0001/CCI_MF_x_WT_filtered_MF.pdf", width = 10, height = 12)
new_plot_cci(MF@graphs$MF_x_WT_filtered,
             colors = MF@colors,
             coords = MF@coords[igraph::V(MF@graphs$MF_x_WT_filtered)$name, ],
             plt_name = "MF_x_WT_filtered",
             node_label_size = 0.7,
             pg_node_size_high = 30,
             pg = MF@rankings[["MF_x_WT_filtered"]]$Pagerank)

dev.off()

names(MF@graphs)

############################
#volcano 

pdf(file = "C:\\Users\\laris\\Documents\\Studium\\Masterarbeit\\RCC_results\\plots\\crosstalker\\no_subset_glom_removed\\Volcano_T_x_N.pdf", width = 10, height = 12)
#pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/Volcano_T_x_N.pdf", width = 14, height = 12)
EnhancedVolcano(RCC@stats$T_x_N,
                lab = RCC@stats$T_x_N$columns_name,
                x = "lodds",
                y = "p",
                pCutoff = 0.05,
                labSize = 4.0,
                colAlpha = 0.3,
                ylim = c(0, 10))

dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_05/Volcano_MF_x_WT.pdf", width = 14, height = 12)
EnhancedVolcano(MF@stats$MF_x_WT,
                lab = MF@stats$MF_x_WT$columns_name,
                x = "lodds",
                y = "p",
                pCutoff = 0.05,
                labSize = 4.0,
                colAlpha = 0.4
                #ylim = c(0, 10)
                )

dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/Heatmap_T_x_N.pdf", width = 14, height = 10)
plot_Heatmap(graph = RCC@graphs$T_x_N, weight = "LRScore")
dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/Heatmap_T.pdf", width = 14, height = 10)
plot_Heatmap(graph = RCC@graphs$T, weight = "LRScore")
dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/RCC_results/plots/crosstalker/no_subset_with_cycling/Heatmap_N.pdf", width = 14, height = 10)
plot_Heatmap(graph = RCC@graphs$N, weight = "LRScore")
dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_05/Heatmap_MF_x_WT.pdf", width = 14, height = 10)
plot_Heatmap(graph = MF@graphs$MF_x_WT, weight = "LRScore")
dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_0001/Heatmap_MF.pdf", width = 14, height = 10)
plot_Heatmap(graph = MF@graphs$MF, weight = "LRScore")
dev.off()

pdf(file = "/home/larissa/Documents/Masterarbeit/MF_results/plots/crosstalker/u10_filtered/pval0_0001/Heatmap_WT.pdf", width = 14, height = 10)
plot_Heatmap(graph = MF@graphs$WT, weight = "LRScore")
dev.off()

####################

