library(shiny)
library(sleuth)
library(DT)
library(plotly)
library(data.table)

so <- sleuth_load(file.path("data/sleuth_object.so"))
wald_test <- colnames(design_matrix(so))[2]
condition <- colnames(so$sample_to_covariates)[2]
sleuth_table <- sleuth_results(so, wald_test, 'wt', show_all = FALSE)
default_columns <- c("target_id", "qval", "b")
if("ens_gene" %in% colnames(sleuth_table)){
    default_columns <- c(default_columns, "ens_gene")
} 

# Define UI -----
ui <- fluidPage(
    navbarPage("Differential analysis",
        tabPanel("Results",
            sidebarLayout(
                sidebarPanel(
                    helpText(h3("Differential analysis results"),
			p("target_id: transcript name"),
			p("pval: p-value of the chosen model"),
			p("qval: false discovery rate adjusted p-vaue"),
			p("b: beta value (effect size). Technically a biased estimator of the fold change"),
			p("se_b: standard error of the beta"),
			p("mean_obs: mean of natural log counts of observations"),
			p("var_obs: variance of observation"),
			p("tech_var: technical variance of observation from the bootstraps"),
			p("sigma_sq: raw estimator of the variance once the technical variance has been reemoved"),
			p("smooth_sigma_sq: smooth regression fit for the shrinkage estimation"),
			p("final_sigma_sq: max(sigma_sq, smooth_sigma_sq); used for covariance estimation of beta"),
                        p("If gene names were added, then"),
                        p("ens_gene: gene name in ensembl"),
                        p("ext_gene: external gene name")
		    ),
		    checkboxGroupInput("show_vars", "Table columns to display", colnames(sleuth_table), default_columns) 
                ),
 
                mainPanel(
                    DT::dataTableOutput("sleuth_results"),
                    downloadButton("downloadData", "Download", "show_vars")
                )
            )
        ),

        tabPanel("Bootstrap",
            sidebarLayout(
                sidebarPanel(
                    selectInput("transcript",
                        label = "Transcript",
                        choices = sleuth_table[1],
                        selected = ""),
		    if("ext_gene" %in% colnames(sleuth_table)){
                        checkboxInput("show_genes", "Show genes", value = FALSE)	
                    }
                ),
        
                mainPanel(plotOutput("bootstrap", width = "800px", height = "600px"))
            )
        ),

        tabPanel("PCA",
            sidebarLayout(
                sidebarPanel(
                    helpText("Hover over each point to know its corresponding sample.")
                ),

                mainPanel(plotlyOutput("pca", width = "800px", height = "600px"))
            )
        ),

        tabPanel("Volcano",
            sidebarLayout(
                sidebarPanel(
                    textInput("target_id", "Search for transcript")
                ),

                mainPanel(plotlyOutput("volcano", width = "800px", height = "600px"))
            )
        ),

        tabPanel("Loadings",
            sidebarLayout(
                sidebarPanel(
                    helpText("Loadings plots show the transcripts that exhibit the most variability for each principal component"),
                    selectInput("pc_choice",
                        label = "Principal Component",
                        choices = c(1, 2, 3, 4, 5),
                        selected = "")	
                ),

                mainPanel(plotOutput("pc_loadings", width = "800px", height = "600px"))
            )
        ),

        tabPanel("Sample Heatmap",
            sidebarLayout(
                sidebarPanel(
                    helpText("This is a sample to sample heatmap or distance matrix.")
                ),

                mainPanel(plotOutput("heatmap", width = "800px", height = "600px"))
            )
        )
    )
)


# Define server logic -----
server <- function(input, output) {
    output$sleuth_results <- DT::renderDataTable({
        DT::datatable(sleuth_table[, input$show_vars, drop = FALSE], filter = 'top')
    })

    output$downloadData <- downloadHandler(
        filename = "sleuth_results.csv",
        content = function(file) {
           write.csv(sleuth_table[input[["sleuth_results_rows_all"]], input$show_vars, drop = FALSE], file)
        }  
    )    

    output$bootstrap <- renderPlot({
        p <-plot_bootstrap(so, input$transcript)
        if(input$show_genes) {
            p$labels$subtitle <- sleuth_table[sleuth_table$target_id == input$transcript, ]$ext_gene
        } else {
            p$labels$subtitle <- NULL
        }
        p
    })

    output$pca <- renderPlotly({
        p <- plot_pca(so, color_by = condition)
        p <- p + aes(text=paste0("Sample: ", sample))
        p <- ggplotly(p, tooltip="text")
    })

    volcano_p <- plot_volcano(so, wald_test)
    volcano_p$data$significant <- as.character(volcano_p$data$significant)
    output$volcano <- renderPlotly({
	if(!is.null(input$target_id) && input$target_id != "") {
            volcano_p$data <- volcano_p$data %>% filter(target_id %like% input$target_id)
        }
	plot_ly(volcano_p$data, x=~b, y=~-log10(qval), type="scattergl", color=~significant, colors=c("black", volcano_p$plot_env$sig_color), alpha=volcano_p$plot_env$point_alpha, text=~target_id, hoverinfo = 'text')
    })

    output$pc_loadings <- renderPlot({
        plot_loadings(so, pc_input = as.numeric(input$pc_choice))
    })

    output$heatmap <- renderPlot({
        plot_sample_heatmap(so)
    })
}


# Run the app -----
shinyApp(ui = ui, server = server)
