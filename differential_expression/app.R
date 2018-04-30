library(shiny)
library(sleuth)
library(DT)

so <- sleuth_load(file.path("data/sleuth_object.so"))
wald_test <- colnames(design_matrix(so))[2]
sleuth_table <- sleuth_results(so, wald_test, 'wt', show_all = FALSE)


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
			p("final_sigma_sq: max(sigma_sq, smooth_sigma_sq); used for covariance estimation of beta")
			),
			checkboxGroupInput("show_vars", "Table columns to display", list("target_id", "pval", "qval", "b", "se_b", "mean_obs", "var_obs", "tech_var", "sigma_sq", "smooth_sigma_sq", "final_sigma_sq"),list("target_id", "qval", "b")) 
                ),
 
                mainPanel(DT::dataTableOutput("sleuth_results"), downloadButton("downloadData", "Download", "show_vars"))
            )
        ),

        tabPanel("Bootstrap",
            sidebarLayout(
                sidebarPanel(
                    selectInput("transcript",
                        label = "Transcript",
                        choices = sleuth_table[1],
                        selected = "")	
                ),
        
                mainPanel(plotOutput("bootstrap", width = "800px", height = "600px"))
            )
        ),

        tabPanel("PCA",
            sidebarLayout(
                sidebarPanel(
                    helpText("This is a PCA plot")
                ),

                mainPanel(plotOutput("pca", width = "800px", height = "600px"))
            )
        ),

        tabPanel("Volcano",
            sidebarLayout(
                sidebarPanel(
                    helpText("This is a volcano plot")
                ),

                mainPanel(plotOutput("volcano", width = "800px", height = "600px"))
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
        )
    )
)


# Define server logic -----
server <- function(input, output) {
    output$sleuth_results <- DT::renderDataTable({
        DT::datatable(sleuth_table[, input$show_vars, drop = FALSE])
    })

    output$downloadData <- downloadHandler(
        filename = "sleuth_results.csv",
        content = function(file) {
           write.csv(sleuth_table[, input$show_vars, drop = FALSE], file)
        }  
    )    

    output$bootstrap <- renderPlot({
        plot_bootstrap(so, input$transcript)
    })

    output$pca <- renderPlot({
        plot_pca(so, color_by = "condition_1")
    })

    output$volcano <- renderPlot({
        plot_volcano(so, wald_test)
    })

    output$pc_loadings <- renderPlot({
        plot_loadings(so, pc_input = as.numeric(input$pc_choice))
        #plot_loadings(so)
    })
}


# Run the app -----
shinyApp(ui = ui, server = server)
