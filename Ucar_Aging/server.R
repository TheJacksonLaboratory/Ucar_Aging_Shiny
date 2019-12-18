# 
# R Shiny application server logic for [provide paper citation here!]
# 
# Author: Dave Mellert
# Institution: The Jackson Laboratory
# Date: July 11, 2019


library(shiny)
library(ggplot2)
library(dplyr)
library(pals)

#Data import
input_data <- readRDS("./fig4_barplot_data.RDS") #precalculated ATAC-seq and RNA-seq data
gene_list <- readRDS("./gene_list.RDS") #list of genes, derived from input data
input_data$sex <- factor(input_data$sex)
input_data$age.group <- factor(input_data$age.group)
flow_data <- readRDS('./flow_data.RDS')
flow_data$sex <- factor(flow_data$sex)
ts_ATAC_lookup <- readRDS('./ts_ATAC_lookup_pruned.RDS') #is peak closest to gene part of a temporal peak set?
ts_RNA_lookup <- readRDS('./ts_RNA_lookup_wGene.RDS') #is transcript part of a temporally regulated set?

#Menu options
data_source_options <- c("ATAC-seq", "RNA-seq")
plot_options <- c("boxplot by age group", "scatterplot")
cell_type_options <- c("CD4T", "CD4TNAIVE",
                       "CD4TCENTRALMEMORY", "CD4TEFFECTORMEMORY",
                       "CD4TEFFECTOR", "CD8T",
                       "CD8TNAIVE", "CD8TCENTRALMEMORY",
                       "CD8TEFFECTORMEMORY", "CD8TEFFECTOR",
                       "CD14", "CD19", "OTHER")

#Figure color options
colors_sex <- c("Females"="#ED6A6F","F"="#ED6A6F","Males"="#829CBD","M"="#829CBD")
colors_sex_darken1 <- sapply(colors_sex,function(h) colorRampPalette(c(h,"black"))(10)[3])
colors_sex_lighten1 <- sapply(colors_sex,function(h) colorRampPalette(c(h,"white"))(10)[3])
darker <- function(clr,step=3) {
  colorRampPalette(c(clr,"black"))(10)[step]
}


############################
#Function definitions
############################

#stats
boxplot_stat <- function(data){
  data_Y <- data[data$age.group == "HY",]
  data_M <- data[data$age.group == "HM",]
  data_O <- data[data$age.group == "HO",]
  
  FYO_filter <- (data$age.group == "HY" & data$sex == "Females") | 
                (data$age.group == "HO" & data$sex == "Females")
  data_FYO <- data[FYO_filter,]
  MYO_filter <- (data$age.group == "HY" & data$sex == "Males") |
                (data$age.group == "HO" & data$sex == "Males")
  data_MYO <- data[MYO_filter,]
  
  w_test_Y <- wilcox.test(adj.score ~ sex, data = data_Y)
  w_test_M <- wilcox.test(adj.score ~ sex, data = data_M)
  w_test_O <- wilcox.test(adj.score ~ sex, data = data_O)
  w_test_FYO <- wilcox.test(adj.score ~ age.group, data = data_FYO)
  w_test_MYO <- wilcox.test(adj.score ~ age.group, data = data_MYO)  
  
  MvF_HY <- w_test_Y$p.value
  MvF_HM <- w_test_M$p.value
  MvF_HO <- w_test_O$p.value
  F_YvO <- w_test_FYO$p.value
  M_YvO <- w_test_MYO$p.value
  
  wilcox_out <- data.frame(MvF_HY, MvF_HM, MvF_HO, F_YvO, M_YvO)
  row.names(wilcox_out) <- c("Wilcoxon test p-value")
  colnames(wilcox_out) <- c("M vs. F (HY)",
                            "M vs. F (HM)",
                            "M vs. F (HO)",
                            "HY vs. HO (F)",
                            "HY vs. HO (M)")
  return(wilcox_out)
}

scatter_stat <- function(data, column_name="adj.score"){
  lm_out <- lm(data[,column_name] ~ age:sex + sex, data = data)
  lm_sum <- summary(lm_out)
  row.names(lm_sum$coefficients) <- c("delete_me",
                                      "Effect of sex (M)",
                                      "Effect of age in females",
                                      "Effect of age in males")
  return(lm_sum$coefficients[2:4,])
}

#reactive text
create_ts_message <- function(gene, data_type){
  if (data_type == 'ATAC-seq'){
    ATAC_lookup_rows <- ts_ATAC_lookup[ts_ATAC_lookup$GeneName == gene,]
    if (nrow(ATAC_lookup_rows) == 0){
      output_text <- ""
    } else if (nrow(ATAC_lookup_rows) == 1){
      output_text <- paste0('* This is a temporal gene in the ATAC-seq data in ', ATAC_lookup_rows$sex)
    } else if (nrow(ATAC_lookup_rows) == 2){
      output_text <- "* This is a temporal gene in the ATAC-seq data in both Males and Females"
    } else {
      output_text <- "AN ERROR HAS OCCURRED"
    }
  } else if (data_type == 'RNA-seq'){
    RNA_lookup_rows <- ts_RNA_lookup[ts_RNA_lookup$GeneName == gene,]
    if (nrow(RNA_lookup_rows) == 0){
      output_text <- ""
    } else if (nrow(RNA_lookup_rows) == 1){
      output_text <- paste0('* This is a temporal gene in the RNA-seq data in ', RNA_lookup_rows$sex)
    } else if (nrow(RNA_lookup_rows) == 2){
      output_text <- "* This is a temporal gene in the RNA-seq data in both Males and Females"
    } else {
      output_text <- "AN ERROR HAS OCCURRED"
    }
  }
  return(output_text)
}

#plots
plot_box <- function(graph_data, ylab, data_var){
  ggplot(graph_data, aes(age.group, adj.score, fill = sex, color = sex)) +
    geom_point(size = 1.8,
               position = position_jitterdodge(dodge.width = 0.75,
                                               jitter.width = 0.1)) +
    geom_boxplot(size = 1.0, alpha = 0.75, outlier.size = 0) +
    labs(x = "Age group",
         y = ylab,
         title = data_var) +
    facet_wrap(~GeneName) +
    scale_fill_manual(values = colors_sex, guide = guide_legend(title = NULL)) +
    scale_color_manual(values = colors_sex_darken1, guide = F) +
    theme_bw() +
    theme(aspect.ratio = 1,
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "italic"),
          axis.line = element_line(size = 1.2),
          axis.ticks = element_line(size = 1.5),
          text = element_text(size = 20))
}

plot_scatter <- function(graph_data, ylab, data_var){
  ggplot(graph_data, aes(age,adj.score,fill=sex,color=sex)) +
    geom_point(size = 1.8,shape = 21,alpha = 0.75,stroke = 0.5) +
    geom_smooth(size=1,alpha = 0.75,
                method = "lm",
                se = T,
                show.legend = F,
                level = .68, #shaded area in regressions are 1 standard error
                fullrange = T) +
    labs(x = "Age, yo",
         y = ylab,
         title = data_var) +
    facet_wrap(~GeneName) +
    scale_fill_manual(values = colors_sex_lighten1,
                      guide = guide_legend(title = NULL, override.aes = list(size=2))) +
    scale_color_manual(values = sapply(colors_sex, darker, 2), guide=F) +
    theme_bw() +
    theme(aspect.ratio = 1,
          strip.background = element_blank(),
          strip.text = element_text(hjust = 0, face = "italic"),
          axis.line = element_line(size = 1.2),
          axis.ticks = element_line(size = 1.5),
          text = element_text(size = 20))
}

plot_flow <- function(flow_data, y){ggplot(flow_data, aes_string(x="age",y=y,fill="sex",color="sex")) +
  geom_point(size = 1.8,shape = 21,alpha = 0.75,stroke = 0.5) +
  geom_smooth(size=1,alpha = 0.75,
              method = "lm",
              se = T,
              show.legend = F,
              level = .68, #shaded area in regressions are 1 standard error
              fullrange = T) +
  labs(x = "Age, yo",
       y = "% of PBMCs",
       title = y) +
  scale_fill_manual(values = colors_sex_lighten1,
                    guide = guide_legend(title = NULL, override.aes = list(size=2))) +
  scale_color_manual(values = sapply(colors_sex, darker, 2), guide=F) +
  theme_bw() +
  theme(aspect.ratio = 1,
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0, face = "italic"),
        axis.line = element_line(size = 1.2),
        axis.ticks = element_line(size = 1.5),
        text = element_text(size = 20))
}


###########################
#Define server logic by tab
###########################

shinyServer(function(input, output) {

  #TAB: GENE-----------------
  
  #Input variables
  output$gene_choice = renderUI({
    selectizeInput("gene_var", "Gene", gene_list, options = list(maxOptions = 100000))
  })
  output$data_choice = renderUI({
    selectInput("data_var", "Data Type", data_source_options)
  })
  output$plot_choice = renderUI({
    selectInput("plot_var", "Plot Type", plot_options)
  })
  
  #Main Plot
  output$gene_main_plot <- renderPlot({
    graph_data <- input_data %>% filter(DataSource == input$data_var, GeneName == input$gene_var)
    if (input$data_var == "ATAC-seq") {
      ylab = "Normalized accessibility"
    } else if (input$data_var == "RNA-seq") {
      ylab = "Normalized expression"
    }
    if (input$plot_var == "boxplot by age group"){
      plot_box(graph_data, ylab, input$data_var)
    } else if (input$plot_var == "scatterplot"){
      plot_scatter(graph_data, ylab, input$data_var)
    }
  })

  #stats table
  output$gene_stats_table <- renderTable({
    stat_data <- input_data %>% filter(DataSource == input$data_var, GeneName == input$gene_var)
    if (input$plot_var == "boxplot by age group"){
      boxplot_stat(stat_data) 
    } else if (input$plot_var == "scatterplot"){
      scatter_stat(stat_data)
    }
  },
  rownames=TRUE,
  digits=5)
  
  #reactive text
  output$ts_report_out <- renderText({
    create_ts_message(input$gene_var, input$data_var)
  })
  
  #TAB: FLOW DATA--------------
 
  #Input variables
  output$cell_choice = renderUI({
    selectInput("cell_var", "Cell Type", cell_type_options)
  })
  
  #Main plot
  output$flow_main_plot <- renderPlot({
    plot_flow(flow_data, y = input$cell_var)
  })
  
  #stats table
  output$flow_stats_table <- renderTable({
    scatter_stat(flow_data, input$cell_var)
    },
  rownames=TRUE,
  digits=5)
})
