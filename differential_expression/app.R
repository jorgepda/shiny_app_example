library(shiny)
library(sleuth)
library(DT)
library(plotly)
library(data.table)
library(dplyr)

so <- sleuth_load(file.path("data/sleuth_object.so"))
wald_test <- colnames(design_matrix(so))[2]
condition <- colnames(so$sample_to_covariates)[2]
sleuth_table <- sleuth_results(so, wald_test, 'wt', show_all = FALSE)
#sleuth_table %>% 
# mutate_if(is.numeric, round, digits = 5)
sleuth_table <- data.frame(lapply(sleuth_table, function(y) if(is.numeric(y)) signif(y, 5) else y))
table_columns <- c("target_id: transcript name" = "target_id",
    "pval: p-value of the chosen model" = "pval",
    "qval: false discovery rate adjusted p-vaue" = "pval",
    "b: beta value (effect size). Technically a biased estimator of the fold change" = "b",
    "se_b: standard error of the beta" = "se_b",
    "mean_obs: mean of natural log counts of observations" = "mean_obs",
    "var_obs: variance of observation" = "var_obs",
    "tech_var: technical variance of observation from the bootstraps" = "tech_var",
    "sigma_sq: raw estimator of the variance once the technical variance has been reemoved" = "sigma_sq",
    "smooth_sigma_sq: smooth regression fit for the shrinkage estimation" = "smooth_sigma_sq",
    "final_sigma_sq: max(sigma_sq, smooth_sigma_sq); used for covariance estimation of beta" = "final_sigma_sq")
default_columns <- c("target_id", "qval", "b")
if("ens_gene" %in% colnames(sleuth_table)){
    table_columns <- c(table_columns,c("ens_gene: gene name in ensembl" = "ens_gene",
    	"ext_gene: external gene name" = "ext_gene"))
    default_columns <- c(default_columns, "ext_gene")
} 

# Define UI -----
ui <- fluidPage(
    navbarPage("Differential analysis",
        tabPanel("Results",
            sidebarLayout(
                sidebarPanel(
                    helpText(h3("Differential analysis resuts")),
		    checkboxGroupInput("show_vars", "Table columns to display", table_columns, default_columns) 
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
                    },
                    helpText(p("Boxplots showing transcript abundances and technical variation for each sample"))
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
        if("ext gene" %in% colnames(sleuth_table) && input$show_genes) {
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
