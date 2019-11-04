# 
# R Shiny application UI definition for [provide paper citation here!]
# 
# Author: Dave Mellert
# Institution: The Jackson Laboratory
# Date: July 11, 2019

library(shiny)

shinyUI(
  navbarPage("Ucar Human Immune Aging Data",
    tabPanel("Explore Data by Gene",
      fluidPage(
        titlePanel("Gene Selection"),
        
        # Sidebar with dropdowns
        sidebarLayout(
          #gene selection
          sidebarPanel(
             uiOutput("gene_choice"),
             uiOutput("data_choice"),
             uiOutput("plot_choice")
             ),
          
          # Plots in main panels
          mainPanel(
             plotOutput("gene_main_plot"),
             tableOutput('gene_stats_table'),
             #textOutput('ts_report_out')
             span(textOutput('ts_report_out'), style="color:red")
             )
          )
        )
      ),
    tabPanel("Explore Flow Cytometry Data",
      fluidPage(
        titlePanel("Flow Cytometry Data"),
        sidebarLayout(
          sidebarPanel(
            uiOutput("cell_choice")
          ),
          mainPanel(
            plotOutput("flow_main_plot"),
            tableOutput('flow_stats_table')
          )
        )
      )
    )
  )
)
