library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
RNAconf = readRDS("RNAconf.rds")
RNAdef  = readRDS("RNAdef.rds")



ATACconf = readRDS("ATACconf.rds")
ATACdef  = readRDS("ATACdef.rds")



Peaksconf = readRDS("Peaksconf.rds")
Peaksdef  = readRDS("Peaksdef.rds")



### Start server code 
shinyUI(fluidPage( 
### HTML formatting of error messages 
 
tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))), 
 
   
### Page title 
titlePanel("Lmxb1 multiome"),  
navbarPage( 
  NULL,  
 navbarMenu("RNA-SEQ",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("RNAa1drX", "X-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                           selected = RNAdef$dimred[1]), 
            selectInput("RNAa1drY", "Y-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                        selected = RNAdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("RNAa1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.RNAa1togL % 2 == 1", 
          selectInput("RNAa1sub1", "Cell information to subset:", 
                      choices = RNAconf[grp == TRUE]$UI, 
                      selected = RNAdef$grp1), 
          uiOutput("RNAa1sub1.ui"), 
          actionButton("RNAa1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("RNAa1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("RNAa1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.RNAa1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("RNAa1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("RNAa1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("RNAa1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("RNAa1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("RNAa1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("RNAa1inp1", "Cell information:", 
                           choices = RNAconf$UI, 
                           selected = RNAdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("RNAa1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.RNAa1tog1 % 2 == 1", 
              radioButtons("RNAa1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("RNAa1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("RNAa1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("RNAa1oup1.ui"))), 
        downloadButton("RNAa1oup1.pdf", "Download PDF"), 
        downloadButton("RNAa1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("RNAa1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("RNAa1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("RNAa1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.RNAa1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("RNAa1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("RNAa1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("RNAa1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("RNAa1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.RNAa1tog2 % 2 == 1", 
              radioButtons("RNAa1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("RNAa1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("RNAa1oup2.ui"))), 
        downloadButton("RNAa1oup2.pdf", "Download PDF"), 
        downloadButton("RNAa1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("RNAa1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("RNAa1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("RNAa2drX", "X-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                           selected = RNAdef$dimred[1]), 
            selectInput("RNAa2drY", "Y-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                        selected = RNAdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("RNAa2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.RNAa2togL % 2 == 1", 
          selectInput("RNAa2sub1", "Cell information to subset:", 
                      choices = RNAconf[grp == TRUE]$UI, 
                      selected = RNAdef$grp1), 
          uiOutput("RNAa2sub1.ui"), 
          actionButton("RNAa2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("RNAa2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("RNAa2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.RNAa2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("RNAa2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("RNAa2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("RNAa2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("RNAa2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("RNAa2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("RNAa2inp1", "Cell information:", 
                           choices = RNAconf$UI, 
                           selected = RNAdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("RNAa2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.RNAa2tog1 % 2 == 1", 
              radioButtons("RNAa2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("RNAa2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("RNAa2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("RNAa2oup1.ui"))), 
        downloadButton("RNAa2oup1.pdf", "Download PDF"), 
        downloadButton("RNAa2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("RNAa2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("RNAa2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("RNAa2inp2", "Cell information:", 
                           choices = RNAconf$UI, 
                           selected = RNAdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("RNAa2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.RNAa2tog2 % 2 == 1", 
              radioButtons("RNAa2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("RNAa2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("RNAa2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("RNAa2oup2.ui"))), 
        downloadButton("RNAa2oup2.pdf", "Download PDF"), 
        downloadButton("RNAa2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("RNAa2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("RNAa2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("RNAa3drX", "X-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                           selected = RNAdef$dimred[1]), 
            selectInput("RNAa3drY", "Y-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                        selected = RNAdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("RNAa3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.RNAa3togL % 2 == 1", 
          selectInput("RNAa3sub1", "Cell information to subset:", 
                      choices = RNAconf[grp == TRUE]$UI, 
                      selected = RNAdef$grp1), 
          uiOutput("RNAa3sub1.ui"), 
          actionButton("RNAa3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("RNAa3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("RNAa3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.RNAa3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("RNAa3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("RNAa3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("RNAa3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("RNAa3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("RNAa3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("RNAa3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("RNAa3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.RNAa3tog1 % 2 == 1", 
              radioButtons("RNAa3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("RNAa3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("RNAa3oup1.ui"))), 
        downloadButton("RNAa3oup1.pdf", "Download PDF"), 
        downloadButton("RNAa3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("RNAa3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("RNAa3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("RNAa3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("RNAa3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.RNAa3tog2 % 2 == 1", 
              radioButtons("RNAa3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("RNAa3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("RNAa3oup2.ui"))), 
        downloadButton("RNAa3oup2.pdf", "Download PDF"), 
        downloadButton("RNAa3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("RNAa3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("RNAa3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("RNAb2drX", "X-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                           selected = RNAdef$dimred[1]), 
           selectInput("RNAb2drY", "Y-axis:", choices = RNAconf[dimred == TRUE]$UI, 
                       selected = RNAdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("RNAb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.RNAb2togL % 2 == 1", 
         selectInput("RNAb2sub1", "Cell information to subset:", 
                     choices = RNAconf[grp == TRUE]$UI, 
                    selected = RNAdef$grp1), 
         uiOutput("RNAb2sub1.ui"), 
         actionButton("RNAb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("RNAb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("RNAb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.RNAb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("RNAb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("RNAb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("RNAb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("RNAb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("RNAb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("RNAb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("RNAb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("RNAb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.RNAb2tog1 % 2 == 1", 
         radioButtons("RNAb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("RNAb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("RNAb2oup1.ui"), 
       downloadButton("RNAb2oup1.pdf", "Download PDF"), 
       downloadButton("RNAb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("RNAb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("RNAb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("RNAb2oup2.ui"), 
       downloadButton("RNAb2oup2.pdf", "Download PDF"), 
       downloadButton("RNAb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("RNAb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("RNAc1inp1", "Cell information (X-axis):", 
                   choices = RNAconf[grp == TRUE]$UI, 
                   selected = RNAdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("RNAc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("RNAc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("RNAc1pts", "Show data points", value = FALSE), 
       actionButton("RNAc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.RNAc1togL % 2 == 1", 
         selectInput("RNAc1sub1", "Cell information to subset:", 
                     choices = RNAconf[grp == TRUE]$UI, 
                     selected = RNAdef$grp1), 
         uiOutput("RNAc1sub1.ui"), 
         actionButton("RNAc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("RNAc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("RNAc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.RNAc1tog % 2 == 1", 
         sliderInput("RNAc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("RNAc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("RNAc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("RNAc1oup.ui"),  
            downloadButton("RNAc1oup.pdf", "Download PDF"),  
            downloadButton("RNAc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("RNAc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("RNAc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("RNAc2inp1", "Cell information to plot (X-axis):", 
                  choices = RNAconf[grp == TRUE]$UI, 
                  selected = RNAdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("RNAc2inp2", "Cell information to group / colour by:", 
                  choices = RNAconf[grp == TRUE]$UI, 
                  selected = RNAdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("RNAc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("RNAc2flp", "Flip X/Y", value = FALSE), 
      actionButton("RNAc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.RNAc2togL % 2 == 1", 
        selectInput("RNAc2sub1", "Cell information to subset:", 
                    choices = RNAconf[grp == TRUE]$UI, 
                    selected = RNAdef$grp1), 
        uiOutput("RNAc2sub1.ui"), 
        actionButton("RNAc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("RNAc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("RNAc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.RNAc2tog % 2 == 1", 
        radioButtons("RNAc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("RNAc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("RNAc2oup.ui"),  
           downloadButton("RNAc2oup.pdf", "Download PDF"),  
           downloadButton("RNAc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("RNAc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("RNAc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("RNAd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(RNAdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("RNAd1grp", "Group by:", 
                    choices = RNAconf[grp == TRUE]$UI, 
                    selected = RNAconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("RNAd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("RNAd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("RNAd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("RNAd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("RNAd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.RNAd1togL % 2 == 1", 
          selectInput("RNAd1sub1", "Cell information to subset:", 
                      choices = RNAconf[grp == TRUE]$UI, 
                      selected = RNAdef$grp1), 
          uiOutput("RNAd1sub1.ui"), 
          actionButton("RNAd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("RNAd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("RNAd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.RNAd1tog % 2 == 1", 
          radioButtons("RNAd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("RNAd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("RNAd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("RNAd1oupTxt")), 
             uiOutput("RNAd1oup.ui"), 
             downloadButton("RNAd1oup.pdf", "Download PDF"), 
             downloadButton("RNAd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("RNAd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("RNAd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("ATAC-SEQ",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("ATACa1drX", "X-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                           selected = ATACdef$dimred[1]), 
            selectInput("ATACa1drY", "Y-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                        selected = ATACdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("ATACa1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.ATACa1togL % 2 == 1", 
          selectInput("ATACa1sub1", "Cell information to subset:", 
                      choices = ATACconf[grp == TRUE]$UI, 
                      selected = ATACdef$grp1), 
          uiOutput("ATACa1sub1.ui"), 
          actionButton("ATACa1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("ATACa1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("ATACa1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.ATACa1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("ATACa1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("ATACa1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("ATACa1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("ATACa1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("ATACa1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("ATACa1inp1", "Cell information:", 
                           choices = ATACconf$UI, 
                           selected = ATACdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("ATACa1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.ATACa1tog1 % 2 == 1", 
              radioButtons("ATACa1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("ATACa1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("ATACa1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("ATACa1oup1.ui"))), 
        downloadButton("ATACa1oup1.pdf", "Download PDF"), 
        downloadButton("ATACa1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("ATACa1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("ATACa1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("ATACa1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.ATACa1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("ATACa1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("ATACa1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("ATACa1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("ATACa1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.ATACa1tog2 % 2 == 1", 
              radioButtons("ATACa1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("ATACa1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("ATACa1oup2.ui"))), 
        downloadButton("ATACa1oup2.pdf", "Download PDF"), 
        downloadButton("ATACa1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("ATACa1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("ATACa1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("ATACa2drX", "X-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                           selected = ATACdef$dimred[1]), 
            selectInput("ATACa2drY", "Y-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                        selected = ATACdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("ATACa2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.ATACa2togL % 2 == 1", 
          selectInput("ATACa2sub1", "Cell information to subset:", 
                      choices = ATACconf[grp == TRUE]$UI, 
                      selected = ATACdef$grp1), 
          uiOutput("ATACa2sub1.ui"), 
          actionButton("ATACa2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("ATACa2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("ATACa2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.ATACa2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("ATACa2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("ATACa2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("ATACa2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("ATACa2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("ATACa2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("ATACa2inp1", "Cell information:", 
                           choices = ATACconf$UI, 
                           selected = ATACdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("ATACa2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.ATACa2tog1 % 2 == 1", 
              radioButtons("ATACa2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("ATACa2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("ATACa2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("ATACa2oup1.ui"))), 
        downloadButton("ATACa2oup1.pdf", "Download PDF"), 
        downloadButton("ATACa2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("ATACa2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("ATACa2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("ATACa2inp2", "Cell information:", 
                           choices = ATACconf$UI, 
                           selected = ATACdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("ATACa2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.ATACa2tog2 % 2 == 1", 
              radioButtons("ATACa2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("ATACa2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("ATACa2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("ATACa2oup2.ui"))), 
        downloadButton("ATACa2oup2.pdf", "Download PDF"), 
        downloadButton("ATACa2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("ATACa2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("ATACa2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("ATACa3drX", "X-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                           selected = ATACdef$dimred[1]), 
            selectInput("ATACa3drY", "Y-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                        selected = ATACdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("ATACa3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.ATACa3togL % 2 == 1", 
          selectInput("ATACa3sub1", "Cell information to subset:", 
                      choices = ATACconf[grp == TRUE]$UI, 
                      selected = ATACdef$grp1), 
          uiOutput("ATACa3sub1.ui"), 
          actionButton("ATACa3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("ATACa3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("ATACa3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.ATACa3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("ATACa3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("ATACa3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("ATACa3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("ATACa3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("ATACa3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("ATACa3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("ATACa3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.ATACa3tog1 % 2 == 1", 
              radioButtons("ATACa3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("ATACa3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("ATACa3oup1.ui"))), 
        downloadButton("ATACa3oup1.pdf", "Download PDF"), 
        downloadButton("ATACa3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("ATACa3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("ATACa3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("ATACa3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("ATACa3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.ATACa3tog2 % 2 == 1", 
              radioButtons("ATACa3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("ATACa3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("ATACa3oup2.ui"))), 
        downloadButton("ATACa3oup2.pdf", "Download PDF"), 
        downloadButton("ATACa3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("ATACa3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("ATACa3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("ATACb2drX", "X-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                           selected = ATACdef$dimred[1]), 
           selectInput("ATACb2drY", "Y-axis:", choices = ATACconf[dimred == TRUE]$UI, 
                       selected = ATACdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("ATACb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.ATACb2togL % 2 == 1", 
         selectInput("ATACb2sub1", "Cell information to subset:", 
                     choices = ATACconf[grp == TRUE]$UI, 
                    selected = ATACdef$grp1), 
         uiOutput("ATACb2sub1.ui"), 
         actionButton("ATACb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("ATACb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("ATACb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.ATACb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("ATACb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("ATACb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("ATACb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("ATACb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("ATACb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("ATACb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("ATACb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("ATACb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.ATACb2tog1 % 2 == 1", 
         radioButtons("ATACb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("ATACb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("ATACb2oup1.ui"), 
       downloadButton("ATACb2oup1.pdf", "Download PDF"), 
       downloadButton("ATACb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("ATACb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("ATACb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("ATACb2oup2.ui"), 
       downloadButton("ATACb2oup2.pdf", "Download PDF"), 
       downloadButton("ATACb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("ATACb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("ATACc1inp1", "Cell information (X-axis):", 
                   choices = ATACconf[grp == TRUE]$UI, 
                   selected = ATACdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("ATACc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("ATACc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("ATACc1pts", "Show data points", value = FALSE), 
       actionButton("ATACc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.ATACc1togL % 2 == 1", 
         selectInput("ATACc1sub1", "Cell information to subset:", 
                     choices = ATACconf[grp == TRUE]$UI, 
                     selected = ATACdef$grp1), 
         uiOutput("ATACc1sub1.ui"), 
         actionButton("ATACc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("ATACc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("ATACc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.ATACc1tog % 2 == 1", 
         sliderInput("ATACc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("ATACc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("ATACc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("ATACc1oup.ui"),  
            downloadButton("ATACc1oup.pdf", "Download PDF"),  
            downloadButton("ATACc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("ATACc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("ATACc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("ATACc2inp1", "Cell information to plot (X-axis):", 
                  choices = ATACconf[grp == TRUE]$UI, 
                  selected = ATACdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("ATACc2inp2", "Cell information to group / colour by:", 
                  choices = ATACconf[grp == TRUE]$UI, 
                  selected = ATACdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("ATACc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("ATACc2flp", "Flip X/Y", value = FALSE), 
      actionButton("ATACc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.ATACc2togL % 2 == 1", 
        selectInput("ATACc2sub1", "Cell information to subset:", 
                    choices = ATACconf[grp == TRUE]$UI, 
                    selected = ATACdef$grp1), 
        uiOutput("ATACc2sub1.ui"), 
        actionButton("ATACc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("ATACc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("ATACc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.ATACc2tog % 2 == 1", 
        radioButtons("ATACc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("ATACc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("ATACc2oup.ui"),  
           downloadButton("ATACc2oup.pdf", "Download PDF"),  
           downloadButton("ATACc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("ATACc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("ATACc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("ATACd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(ATACdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("ATACd1grp", "Group by:", 
                    choices = ATACconf[grp == TRUE]$UI, 
                    selected = ATACconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("ATACd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("ATACd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("ATACd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("ATACd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("ATACd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.ATACd1togL % 2 == 1", 
          selectInput("ATACd1sub1", "Cell information to subset:", 
                      choices = ATACconf[grp == TRUE]$UI, 
                      selected = ATACdef$grp1), 
          uiOutput("ATACd1sub1.ui"), 
          actionButton("ATACd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("ATACd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("ATACd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.ATACd1tog % 2 == 1", 
          radioButtons("ATACd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("ATACd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("ATACd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("ATACd1oupTxt")), 
             uiOutput("ATACd1oup.ui"), 
             downloadButton("ATACd1oup.pdf", "Download PDF"), 
             downloadButton("ATACd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("ATACd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("ATACd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

navbarMenu("Peaks",### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("CellInfo vs GeneExpr"), 
    h4("Cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("Peaksa1drX", "X-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                           selected = Peaksdef$dimred[1]), 
            selectInput("Peaksa1drY", "Y-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                        selected = Peaksdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("Peaksa1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.Peaksa1togL % 2 == 1", 
          selectInput("Peaksa1sub1", "Cell information to subset:", 
                      choices = Peaksconf[grp == TRUE]$UI, 
                      selected = Peaksdef$grp1), 
          uiOutput("Peaksa1sub1.ui"), 
          actionButton("Peaksa1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("Peaksa1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("Peaksa1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.Peaksa1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("Peaksa1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("Peaksa1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("Peaksa1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("Peaksa1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("Peaksa1txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("Peaksa1inp1", "Cell information:", 
                           choices = Peaksconf$UI, 
                           selected = Peaksdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("Peaksa1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.Peaksa1tog1 % 2 == 1", 
              radioButtons("Peaksa1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("Peaksa1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("Peaksa1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("Peaksa1oup1.ui"))), 
        downloadButton("Peaksa1oup1.pdf", "Download PDF"), 
        downloadButton("Peaksa1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("Peaksa1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("Peaksa1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)), br(), 
        actionButton("Peaksa1tog9", "Toggle to show cell numbers / statistics"), 
        conditionalPanel( 
          condition = "input.Peaksa1tog9 % 2 == 1", 
          h4("Cell numbers / statistics"), 
          radioButtons("Peaksa1splt", "Split continuous cell info into:", 
                       choices = c("Quartile", "Decile"), 
                       selected = "Decile", inline = TRUE), 
          dataTableOutput("Peaksa1.dt") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("Peaksa1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("Peaksa1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.Peaksa1tog2 % 2 == 1", 
              radioButtons("Peaksa1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("Peaksa1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("Peaksa1oup2.ui"))), 
        downloadButton("Peaksa1oup2.pdf", "Download PDF"), 
        downloadButton("Peaksa1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("Peaksa1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("Peaksa1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
  ### Tab1.a2: cellInfo vs cellInfo on dimRed 
  tabPanel( 
    HTML("CellInfo vs CellInfo"), 
    h4("Cell information vs cell information on dimension reduction"), 
    "In this tab, users can visualise two cell informations side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("Peaksa2drX", "X-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                           selected = Peaksdef$dimred[1]), 
            selectInput("Peaksa2drY", "Y-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                        selected = Peaksdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("Peaksa2togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.Peaksa2togL % 2 == 1", 
          selectInput("Peaksa2sub1", "Cell information to subset:", 
                      choices = Peaksconf[grp == TRUE]$UI, 
                      selected = Peaksdef$grp1), 
          uiOutput("Peaksa2sub1.ui"), 
          actionButton("Peaksa2sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("Peaksa2sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("Peaksa2tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.Peaksa2tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("Peaksa2siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("Peaksa2psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("Peaksa2fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("Peaksa2asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("Peaksa2txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information 1"), 
        fluidRow( 
          column( 
            6, selectInput("Peaksa2inp1", "Cell information:", 
                           choices = Peaksconf$UI, 
                           selected = Peaksdef$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("Peaksa2tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.Peaksa2tog1 % 2 == 1", 
              radioButtons("Peaksa2col1", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("Peaksa2ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("Peaksa2lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("Peaksa2oup1.ui"))), 
        downloadButton("Peaksa2oup1.pdf", "Download PDF"), 
        downloadButton("Peaksa2oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("Peaksa2oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("Peaksa2oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Cell information 2"), 
        fluidRow( 
          column( 
            6, selectInput("Peaksa2inp2", "Cell information:", 
                           choices = Peaksconf$UI, 
                           selected = Peaksdef$meta2) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("Peaksa2tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.Peaksa2tog2 % 2 == 1", 
              radioButtons("Peaksa2col2", "Colour (Continuous data):", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("Peaksa2ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("Peaksa2lab2", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("Peaksa2oup2.ui"))), 
        downloadButton("Peaksa2oup2.pdf", "Download PDF"), 
        downloadButton("Peaksa2oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("Peaksa2oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("Peaksa2oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
   
  ### Tab1.a3: geneExpr vs geneExpr on dimRed 
  tabPanel( 
    HTML("GeneExpr vs GeneExpr"), 
    h4("Gene expression vs gene expression on dimension reduction"), 
    "In this tab, users can visualise two gene expressions side-by-side ", 
    "on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("Peaksa3drX", "X-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                           selected = Peaksdef$dimred[1]), 
            selectInput("Peaksa3drY", "Y-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                        selected = Peaksdef$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, actionButton("Peaksa3togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.Peaksa3togL % 2 == 1", 
          selectInput("Peaksa3sub1", "Cell information to subset:", 
                      choices = Peaksconf[grp == TRUE]$UI, 
                      selected = Peaksdef$grp1), 
          uiOutput("Peaksa3sub1.ui"), 
          actionButton("Peaksa3sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("Peaksa3sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, actionButton("Peaksa3tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.Peaksa3tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("Peaksa3siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("Peaksa3psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("Peaksa3fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("Peaksa3asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("Peaksa3txt", "Show axis text", value = FALSE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Gene expression 1"), 
        fluidRow( 
          column( 
            6, selectInput("Peaksa3inp1", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("Peaksa3tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.Peaksa3tog1 % 2 == 1", 
              radioButtons("Peaksa3col1", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("Peaksa3ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("Peaksa3oup1.ui"))), 
        downloadButton("Peaksa3oup1.pdf", "Download PDF"), 
        downloadButton("Peaksa3oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("Peaksa3oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("Peaksa3oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression 2"), 
        fluidRow( 
          column( 
            6, selectInput("Peaksa3inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("Peaksa3tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.Peaksa3tog2 % 2 == 1", 
              radioButtons("Peaksa3col2", "Colour:", 
                           choices = c("White-Red", "Blue-Yellow-Red", 
                                       "Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("Peaksa3ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("Peaksa3oup2.ui"))), 
        downloadButton("Peaksa3oup2.pdf", "Download PDF"), 
        downloadButton("Peaksa3oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("Peaksa3oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("Peaksa3oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 
 ### Tab1.b2: Gene coexpression plot 
 tabPanel( 
   HTML("Gene coexpression"), 
   h4("Coexpression of two genes on reduced dimensions"), 
   "In this tab, users can visualise the coexpression of two genes ", 
   "on low-dimensional representions.", 
   br(),br(), 
   fluidRow( 
     column( 
       3, h4("Dimension Reduction"), 
       fluidRow( 
         column( 
           12, selectInput("Peaksb2drX", "X-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                           selected = Peaksdef$dimred[1]), 
           selectInput("Peaksb2drY", "Y-axis:", choices = Peaksconf[dimred == TRUE]$UI, 
                       selected = Peaksdef$dimred[2])) 
       ) 
     ), # End of column (6 space) 
     column( 
       3, actionButton("Peaksb2togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.Peaksb2togL % 2 == 1", 
         selectInput("Peaksb2sub1", "Cell information to subset:", 
                     choices = Peaksconf[grp == TRUE]$UI, 
                    selected = Peaksdef$grp1), 
         uiOutput("Peaksb2sub1.ui"), 
         actionButton("Peaksb2sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("Peaksb2sub1non", "Deselect all groups", class = "btn btn-primary") 
       ) 
     ), # End of column (6 space) 
     column( 
       6, actionButton("Peaksb2tog0", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.Peaksb2tog0 % 2 == 1", 
         fluidRow( 
           column( 
             6, sliderInput("Peaksb2siz", "Point size:", 
                            min = 0, max = 4, value = 1.25, step = 0.25), 
             radioButtons("Peaksb2psz", "Plot size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE), 
             radioButtons("Peaksb2fsz", "Font size:", 
                          choices = c("Small", "Medium", "Large"), 
                          selected = "Medium", inline = TRUE) 
           ), 
           column( 
             6, radioButtons("Peaksb2asp", "Aspect ratio:", 
                             choices = c("Square", "Fixed", "Free"), 
                             selected = "Square", inline = TRUE), 
             checkboxInput("Peaksb2txt", "Show axis text", value = FALSE) 
           ) 
         ) 
       ) 
     )  # End of column (6 space) 
   ),   # End of fluidRow (4 space) 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", h4("Gene Expression"), 
       selectInput("Peaksb2inp1", "Gene 1:", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
               title = "Gene expression to colour cells by", 
               content = c("Select gene to colour cells by gene expression", 
                          paste0("- Gene expression are coloured in a ", 
                                 "White-Red colour scheme which can be ", 
                                 "changed in the plot controls"))), 
       selectInput("Peaksb2inp2", "Gene 2:", choices=NULL) %>% 
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene expression to colour cells by", 
                content = c("Select gene to colour cells by gene expression", 
                            paste0("- Gene expression are coloured in a ", 
                                   "White-Blue colour scheme which can be ", 
                                   "changed in the plot controls"))), 
       actionButton("Peaksb2tog1", "Toggle plot controls"), 
       conditionalPanel( 
         condition = "input.Peaksb2tog1 % 2 == 1", 
         radioButtons("Peaksb2col1", "Colour:", 
                      choices = c("Red (Gene1); Blue (Gene2)", 
                                  "Orange (Gene1); Blue (Gene2)", 
                                  "Red (Gene1); Green (Gene2)", 
                                  "Green (Gene1); Blue (Gene2)"), 
                      selected = "Red (Gene1); Blue (Gene2)"), 
         radioButtons("Peaksb2ord1", "Plot order:", 
                      choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                      selected = "Max-1st", inline = TRUE) 
       ) 
     ), # End of column (6 space) 
     column( 
       6, style="border-right: 2px solid black", 
       uiOutput("Peaksb2oup1.ui"), 
       downloadButton("Peaksb2oup1.pdf", "Download PDF"), 
       downloadButton("Peaksb2oup1.png", "Download PNG"), br(), 
       div(style="display:inline-block", 
           numericInput("Peaksb2oup1.h", "PDF / PNG height:", width = "138px", 
                        min = 4, max = 20, value = 8, step = 0.5)), 
       div(style="display:inline-block", 
           numericInput("Peaksb2oup1.w", "PDF / PNG width:", width = "138px", 
                        min = 4, max = 20, value = 10, step = 0.5)) 
     ), # End of column (6 space) 
     column( 
       3, uiOutput("Peaksb2oup2.ui"), 
       downloadButton("Peaksb2oup2.pdf", "Download PDF"), 
       downloadButton("Peaksb2oup2.png", "Download PNG"), 
       br(), h4("Cell numbers"), 
       dataTableOutput("Peaksb2.dt") 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("Violinplot / Boxplot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("Peaksc1inp1", "Cell information (X-axis):", 
                   choices = Peaksconf[grp == TRUE]$UI, 
                   selected = Peaksdef$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("Peaksc1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell Info / Gene to plot", 
                content = c("Select cell info / gene to plot on Y-axis", 
                            "- Can be continuous cell information (e.g. nUMIs / scores)", 
                            "- Can also be gene expression")), 
       radioButtons("Peaksc1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("Peaksc1pts", "Show data points", value = FALSE), 
       actionButton("Peaksc1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = "input.Peaksc1togL % 2 == 1", 
         selectInput("Peaksc1sub1", "Cell information to subset:", 
                     choices = Peaksconf[grp == TRUE]$UI, 
                     selected = Peaksdef$grp1), 
         uiOutput("Peaksc1sub1.ui"), 
         actionButton("Peaksc1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("Peaksc1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(), br(), 
       actionButton("Peaksc1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.Peaksc1tog % 2 == 1", 
         sliderInput("Peaksc1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("Peaksc1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("Peaksc1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("Peaksc1oup.ui"),  
            downloadButton("Peaksc1oup.pdf", "Download PDF"),  
            downloadButton("Peaksc1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("Peaksc1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("Peaksc1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 
### Tab1.c2: Proportion plot 
tabPanel( 
  HTML("Proportion plot"), 
  h4("Proportion / cell numbers across different cell information"), 
  "In this tab, users can visualise the composition of single cells based on one discrete ", 
  "cell information across another discrete cell information. ",  
  "Usage examples include the library or cellcycle composition across clusters.", 
  br(),br(), 
  fluidRow( 
    column( 
      3, style="border-right: 2px solid black", 
      selectInput("Peaksc2inp1", "Cell information to plot (X-axis):", 
                  choices = Peaksconf[grp == TRUE]$UI, 
                  selected = Peaksdef$grp2) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to plot cells by",  
               content = c("Select categorical cell information to plot cells by", 
                           "- Plotted as the X-axis of the proportion plot")), 
      selectInput("Peaksc2inp2", "Cell information to group / colour by:", 
                  choices = Peaksconf[grp == TRUE]$UI, 
                  selected = Peaksdef$grp1) %>%  
        helper(type = "inline", size = "m", fade = TRUE, 
               title = "Cell information to group / colour cells by", 
               content = c("Select categorical cell information to group / colour cells by", 
                           "- Proportion / cell numbers are shown in different colours")), 
      radioButtons("Peaksc2typ", "Plot value:", 
                   choices = c("Proportion", "CellNumbers"), 
                   selected = "Proportion", inline = TRUE), 
      checkboxInput("Peaksc2flp", "Flip X/Y", value = FALSE), 
      actionButton("Peaksc2togL", "Toggle to subset cells"), 
      conditionalPanel( 
        condition = "input.Peaksc2togL % 2 == 1", 
        selectInput("Peaksc2sub1", "Cell information to subset:", 
                    choices = Peaksconf[grp == TRUE]$UI, 
                    selected = Peaksdef$grp1), 
        uiOutput("Peaksc2sub1.ui"), 
        actionButton("Peaksc2sub1all", "Select all groups", class = "btn btn-primary"), 
        actionButton("Peaksc2sub1non", "Deselect all groups", class = "btn btn-primary") 
      ), br(), br(), 
      actionButton("Peaksc2tog", "Toggle graphics controls"), 
      conditionalPanel( 
        condition = "input.Peaksc2tog % 2 == 1", 
        radioButtons("Peaksc2psz", "Plot size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE), 
        radioButtons("Peaksc2fsz", "Font size:", 
                     choices = c("Small", "Medium", "Large"), 
                     selected = "Medium", inline = TRUE)) 
    ), # End of column (6 space) 
    column(9, uiOutput("Peaksc2oup.ui"),  
           downloadButton("Peaksc2oup.pdf", "Download PDF"),  
           downloadButton("Peaksc2oup.png", "Download PNG"), br(), 
           div(style="display:inline-block", 
               numericInput("Peaksc2oup.h", "PDF / PNG height:", width = "138px", 
                            min = 4, max = 20, value = 8, step = 0.5)), 
           div(style="display:inline-block", 
               numericInput("Peaksc2oup.w", "PDF / PNG width:", width = "138px", 
                            min = 4, max = 20, value = 10, step = 0.5)) 
    )  # End of column (6 space) 
  )    # End of fluidRow (4 space) 
),     # End of tab (2 space) 
 
  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("Bubbleplot / Heatmap"), 
    h4("Gene expression bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. library / cluster).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("Peaksd1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(Peaksdef$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("Peaksd1grp", "Group by:", 
                    choices = Peaksconf[grp == TRUE]$UI, 
                    selected = Peaksconf[grp == TRUE]$UI[1]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        radioButtons("Peaksd1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("Peaksd1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("Peaksd1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("Peaksd1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("Peaksd1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = "input.Peaksd1togL % 2 == 1", 
          selectInput("Peaksd1sub1", "Cell information to subset:", 
                      choices = Peaksconf[grp == TRUE]$UI, 
                      selected = Peaksdef$grp1), 
          uiOutput("Peaksd1sub1.ui"), 
          actionButton("Peaksd1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("Peaksd1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), br(), 
        actionButton("Peaksd1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.Peaksd1tog % 2 == 1", 
          radioButtons("Peaksd1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("Peaksd1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("Peaksd1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("Peaksd1oupTxt")), 
             uiOutput("Peaksd1oup.ui"), 
             downloadButton("Peaksd1oup.pdf", "Download PDF"), 
             downloadButton("Peaksd1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("Peaksd1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("Peaksd1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ), 

   
br(),br(),br(),br(),br() 
))) 
 
 
 
 