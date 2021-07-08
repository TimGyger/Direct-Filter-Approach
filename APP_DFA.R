###############################################################################################################
### Direct Filter Approach
###############################################################################################################

# clear objects from workspace
rm(list = ls())

###########
### Load Packages
###########

listOfPackages <- c('tidyr', 'scales', 'shiny', 'tidyverse', 'plotly','PerformanceAnalytics',
                    'data.table', 'zoo', 'dplyr', 'shinyjs', 'shinycssloaders','tidyquant')

for (i in listOfPackages){
  if(! i %in% installed.packages()){
    install.packages(i, dependencies = TRUE)
  }
  require(i, character.only = TRUE)
}
select <- dplyr::select


# Functions
source("~/Desktop/Quantitative Investment/Functions.R",encoding="latin1")

# Data Load
source("~/Desktop/Quantitative Investment/Data_Load_DFA.R",encoding="latin1")

# Parameter
source("~/Desktop/Quantitative Investment/Parameter.R",encoding="latin1")

# Watch-List
Current.Stocks <- c("AAPL")


# Define UI for app that draws a histogram ----
ui <- fluidPage(
  useShinyjs(),
  # App title ----
  titlePanel("Direct Filter Approach"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 3,h3("Choose Data"),
      
      # Input: Slider for the number of bins ----
      selectizeInput("Stock", label = 'Choose Stock from Watchlist:', choices = Current.Stocks, selected = 1),
      checkboxInput("Stock.Choice","Choose arbitrary Stock",value = FALSE),
      textInput("Stock1", label = 'Choose Stock:', value = "", width = NULL, placeholder = NULL),
      selectizeInput("Freq", label = 'Choose Data-Frequency:', choices = c("4-Daily","3-Daily","2-Daily","Daily"), selected = "Daily"),
      dateRangeInput('dateRange',
                     label = 'Period:',
                     start = "2010-01-01", end = Sys.Date()), 
      actionButton("Reset.Date","Reset Enddate to today"), hr(),
      actionButton("Reload.Button","Reload App")
      ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(
        tabPanel("Stock Data", br(),
                  h2("Chart"),
                  plotlyOutput("Stock_plot") %>% withSpinner(color="#0dc5c1"),
                  h2("Log-Return"),
                  plotlyOutput("Log_plot") %>% withSpinner(color="#0dc5c1"), hr(),
                  h2("Log-Return Periodogram"),
                  plotlyOutput("LogP_plot") %>% withSpinner(color="#0dc5c1")),
        tabPanel("Analysis", 
                 tabsetPanel(
                   tabPanel("Inputs", br(),
                            p("Link to Theory & R-Package:"),
                            uiOutput("tab"), br(),
                            splitLayout(numericInput("L","Filter-Length:",value = 100, min = 3, max = 200),
                                        numericInput("lambda","Time-Parameter:",value = 0, min = 0, max = 40, step = 0.2),
                                        numericInput("eta","Smoothness-Parameter:",value = 0, min = 0, max = 40, step = 0.2),
                                        numericInput("k","Trend-Length in Months:",value = 12, min = 1, step = 1),
                                        checkboxGroupInput("Boolean", label = "Optional:", selected = NULL,
                                                           inline = FALSE, width = NULL, choiceNames = c("Amplitude Constraint","Shift Constraint"), choiceValues = c("i1","i2")),
                                        cellWidths = c("15%","15%","25%","25%","20%")), hr(),
                            checkboxInput("DF","Double Filter",value = TRUE),
                            splitLayout(numericInput("L1","Filter-Length:",value = 3, min = 3, max = 200),
                                        numericInput("lambda1","Time-Parameter:",value = 0, min = 0, max = 40, step = 0.2),
                                        numericInput("eta1","Smoothness-Parameter:",value = 0, min = 0, max = 40, step = 0.2),
                                        numericInput("k1","Trend-Length in Months:",value = 12, min = 1, step = 1),
                                        checkboxGroupInput("Boolean1", label = "Optional:", selected = NULL,
                                                           inline = FALSE, width = NULL, choiceNames = c("Amplitude Constraint","Shift Constraint"), choiceValues = c("i1","i2")),
                                        cellWidths = c("15%","15%","25%","25%","20%")), hr(),
                            h2("Log-Return Filter"),
                            plotlyOutput("LogDFA_plot") %>% withSpinner(color="#0dc5c1"), br(),
                            h2("Trend Signals"),
                            plotlyOutput("Trends") %>% withSpinner(color="#0dc5c1"),
                            h2("Error"),
                            tableOutput("Error_Table") %>% withSpinner(color="#0dc5c1")),
                   tabPanel("Analysis", 
                            h2("Amplitude and Shift"),
                            splitLayout(plotlyOutput("LogA_plot") %>% withSpinner(color="#0dc5c1"),
                            plotlyOutput("LogS_plot") %>% withSpinner(color="#0dc5c1"),cellWidths = c("50%","50%")),
                            h2("Periodogram"),
                            plotlyOutput("LogPDFA_plot") %>% withSpinner(color="#0dc5c1"),
                            h2("Filter-Coefficients"),
                            plotlyOutput("Coefficient_DFA") %>% withSpinner(color="#0dc5c1")),
                   tabPanel("Out of Sample Performance", br(),
                            dateInput(
                              "date.OOS",
                              "Start of Out of Sample Period:",
                              value = NULL,
                              min = NULL,
                              max = NULL),
                            h2("Log-Return and Filter"),
                            plotlyOutput("Out_of_Sample_p") %>% withSpinner(color="#0dc5c1"),
                            h2("Yield and Trends"),
                            plotlyOutput("TrendsM_Out") %>% withSpinner(color="#0dc5c1"), hr(),
                            h1("Performance"),
                            h2("Cross-Correlation"),
                            checkboxInput("OS.Cross","Out of Sample",value = T),
                            plotlyOutput("Cross_Corr") %>% withSpinner(color="#0dc5c1"),
                            h3("Interpretation"),
                            textOutput("Cross_Corr.Text") %>% withSpinner(color="#0dc5c1"),
                            h2("Performance"),
                            checkboxInput("OS.Perf","Out of Sample",value = T),
                            checkboxInput("Fee","Fee per each Transaction",value = F),
                            numericInput("Fee.num","Fee in Percent:",value = 1, min = 0, max = 100),
                            plotlyOutput("Performance_Chart") %>% withSpinner(color="#0dc5c1"),
                            plotlyOutput("OutPerformance_Chart",height = "200px") %>% withSpinner(color="#0dc5c1")))),
        tabPanel("Aktuelle Resultate", 
                 h2("Filter over all Data-Frequencies"),
                            checkboxInput("Zoom","Zoom"),  
                            plotlyOutput("Filter_Combo") %>% withSpinner(color="#0dc5c1"),
                            h2("Signalised Trends of selected Frequenzy"),
                            selectizeInput("Freq1", label = 'Choose Data-Frequency:', choices = c("4-Daily","3-Daily","2-Daily","Daily","Average"), selected = "Average"),
                            plotlyOutput("Trend_Durchschnitt") %>% withSpinner(color="#0dc5c1"),
                            h2("Results"),
                            tableOutput("Resultate") %>% withSpinner(color="#0dc5c1"),hr(),
                            h3("Computation to Cross-Correlation"),
                            plotlyOutput("Cross_Corr_AR"),
                            tableOutput("Cross_Corr.AR"))
      )
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output,session) {
    
  ##########
  #### Links
  ##########
  
  url <- a("Univariate Direct Filter Approach", href="http://signalextractionandforecasting.blogspot.com/2018/02/what-is-direct-filter-approach.html")
  output$tab <- renderUI({
    tagList("", url)})
  
  ###############################
  #### Update Input-Possibilities
  ###############################
  observeEvent(input$Reload.Button, session$reload())
  Chosen.Date <- reactive(input$dateRange[2])
  Chosen.Date1 <- reactive(input$dateRange[1])
  
  observeEvent(input$Stock.Choice, if(input$Stock.Choice){
    shinyjs::show("Stock1")
    }else{hide("Stock1")})
  
  observeEvent(input$Reset.Date,
               updateDateRangeInput(session,"dateRange",label = 'Periode:',
                                    start = Chosen.Date1(), end = Sys.Date()))
  
  
  
  observeEvent(input$dateRange,
    updateDateInput(session,"date.OOS",label = "Start of Out of Sample Period:", 
                    value = mean(c(as.Date(input$dateRange[1]),as.Date(input$dateRange[2])),na.rm = T),
                         min = input$dateRange[1], max = input$dateRange[2]))
  
  observeEvent(input$Freq,if(input$Freq == "3-Daily"){
    updateNumericInput(session,'k',"Trend-Length in 2 Weeks:",value = 4, min = 1, step = 1);
    updateNumericInput(session,'k1',"Trend-Length in 2 Weeks:",value = 4, min = 1, step = 1);
    updateNumericInput(session,'L',"Filter-Length:",value = 100, min = 3, max = 200)})
  
  observeEvent(input$Freq,if(input$Freq == "2-Daily"){
    updateNumericInput(session,'k',"Trend-Length in Weeks:",value = 6, min = 1, step = 1);
    updateNumericInput(session,'k1',"Trend-Length in Weeks:",value = 6, min = 1, step = 1);
    updateNumericInput(session,'L',"Filter-Length:",value = 100, min = 3, max = 200)})
  
  observeEvent(input$Freq,if(input$Freq == "Daily"){
    updateNumericInput(session,'k',"Trend-Length in Days:",value = 12, min = 1, step = 1);
    updateNumericInput(session,'k1',"Trend-Length in Days:",value = 12, min = 1, step = 1);
    updateNumericInput(session,'L',"Filter-Length:",value = 100, min = 3, max = 200)})
  
  observeEvent(input$Freq,if(input$Freq == "4-Daily"){
    updateNumericInput(session,'k',"Trend-Length in Months:",value = 2, min = 1, step = 1);
    updateNumericInput(session,'k1',"Trend-Length in Months:",value = 2, min = 1, step = 1);
    updateNumericInput(session,'L',"Filter-Length:",value = 100, min = 3, max = 200)})
  
  ########
  ### Data
  ########
  
  # Stock
  Stock.Symb <- reactive(if(input$Stock.Choice){
    input$Stock1
  } else {
    input$Stock
  })
  
  # Load
  data.Stock <- reactive(Load.Daten(Stock.Symb(),input$dateRange[1],input$dateRange[2]))
  # Manipulate
  data <- reactive(Choice(data.Stock(),input$Freq))
  data.1 <- reactive(data()[[1]])
  data.2 <- reactive(data()[[2]])
  
  data.OOS <- reactive(Choice(data.Stock(),input$Freq,OOS = TRUE,date2 = input$date.OOS))
  data.OOS1 <- reactive(data.OOS()[[1]])
  data.OOS2 <- reactive(data.OOS()[[2]])
  
  ##########
  ### Charts & Periodogram
  ##########
  
  output$Stock_plot <- renderPlotly(plot_ly(x = data.1()[,1], y = data.1()[,2],
                                            type = "scatter", mode = "lines"))
  
  output$Log_plot <- renderPlotly(plot_ly(x = data.1()[-1,1], y = Skalieren(data.2()[,1]),
                                          type = "scatter", mode = "lines") %>% 
                                    layout(yaxis = list(
                                      title = "",
                                      showticklabels = FALSE
                                    )))
    
  ### Periodogram
  
  output$LogP_plot <- renderPlotly(periodogram.Plot(data.2()[,1],b = TRUE))
  
  
  ######
  ## Inputs
  #####
  # Data
  data.3 <- reactive({
    validate(
      need(input$dateRange[2]-input$dateRange[1] > 6*365, "Bitte wähle die Periode genügend gross (mindestens 6 Jahre)")
    )
    Univariate.Filter.cust(data.2()[,1],0,input$L,1/input$k,input$lambda,input$eta,pi/input$k,input$Boolean)
  })
  
  data.4 <- reactive({
    validate(
      need(input$dateRange[2]-input$dateRange[1] > 6*365, "Bitte wähle die Periode genügend gross (mindestens 6 Jahre)")
    )
    list1()[[1]]
  })
  
  list1 <- reactive({
    validate(
      need(input$dateRange[2]-input$dateRange[1] > 6*365, "Bitte wähle die Periode genügend gross (mindestens 6 Jahre)")
    )
    if(!input$DF){data.3()} else {
    Univariate.Filter.cust.conv(data.2()[,1],0,input$L,1/input$k,input$lambda,input$eta,pi/input$k,input$Boolean,
                                0,input$L1,1/input$k1,input$lambda1,input$eta1,pi/input$k1,input$Boolean1)}
    })
  
  # Error
  output$Error_Table <- renderTable({Table_Error <- as.data.frame(MS_decomp_total(list1()[[5]],
                                                    list1()[[4]],
                                                    list1()[[6]],
                                                    list1()[[7]]))
    Table_Error$Accuracy <- formatC(Table_Error$Accuracy,format = "e",digits = 2)
    Table_Error$Smoothness <- formatC(Table_Error$Smoothness,format = "e",digits = 2)
    Table_Error$Timeliness <- formatC(Table_Error$Timeliness,format = "e",digits = 2)
    Table_Error$`MS-Error` <- formatC(Table_Error$`MS-Error`,format = "e",digits = 2)
    return(Table_Error)})
  
  
  # Filter Plot and Trend Segmente
  output$LogDFA_plot <- renderPlotly({if(input$DF){
    plot_ly(x = data.1()[-1,1], y = data.2()[,1],
            type = "scatter", mode = "lines", name = "Underlying Data") %>%
      layout(yaxis = list(
        title = "",
        showticklabels = FALSE
      )) %>%
      add_lines(y = data.3()[[1]], name = "Direct Filter") %>%
      add_lines(y = Skalieren(data.4(),max(data.3()[[1]],na.rm = T)), name = "Direct Double Filter")
  } else {
    plot_ly(x = data.1()[-1,1], y = data.2()[,1],
            type = "scatter", mode = "lines", name = "Underlying Data") %>%
      layout(yaxis = list(
        title = "",
        showticklabels = FALSE
      )) %>%
      add_lines(y = data.3()[[1]], name = "Direct Filter")}})
  
  output$Trends <- renderPlotly(Segmente.Trend(data.1()[,1],data.2()[,1],data.1()[,2],0,list1()[[8]],FALSE)[[1]])
  ######
  ### Analysis
  ######
  Gamma_k <- reactive(ifelse(!input$DF,1/input$k,1/input$k1))
  
  # Amplitude and Shift
  output$LogA_plot <- renderPlotly(Amplitude(list1()[[2]],Gamma_k()))
  output$LogS_plot <- renderPlotly(Shift(list1()[[3]]))
  
  # Periodogram
  output$LogPDFA_plot <- renderPlotly({if(!input$DF){
    periodogram.Plot(list1()[[1]][(input$L +1):length(data.2()[,1])],data.2()[-input$L,1],TRUE)
  } else {
    periodogram.Plot(list1()[[1]][(input$L1+input$L):length(data.2()[,1])], data.2()[-(input$L1+input$L-1),1],TRUE)
    }})
  
  # Coefficients
  output$Coefficient_DFA <- renderPlotly(plot_ly(x = seq(0,length(list1()[[8]])-1)) %>%
    add_trace(y = list1()[[8]], name = "Coefficients" , type = 'scatter', mode = "lines") %>%
      layout(xaxis = list(title = 'Lag'), yaxis = list(title = '')))
  #####
  ## Out of Sample Performance
  #####
  # 
  
  data.OOS3 <- reactive({
    Univariate.Filter.cust(data.OOS2()[,1],0,input$L,1/input$k,input$lambda,input$eta,pi/input$k,input$Boolean)
  })
  
  data.OOS4 <- reactive({
    listOOS1()[[1]]
  })
  
  listOOS1 <- reactive({
    if(!input$DF){data.OOS3()} else {
      Univariate.Filter.cust.conv(data.OOS2()[,1],0,input$L,1/input$k,input$lambda,input$eta,pi/input$k,input$Boolean,
                                  0,input$L1,1/input$k1,input$lambda1,input$eta1,pi/input$k1,input$Boolean1)}
  })
  
  
  # Filter and Data
  output$Out_of_Sample_p <- renderPlotly(Out_of_Sample_Chart(data()[[1]][-1,1],
                                                             data()[[2]][,1],input$date.OOS,listOOS1()[[8]],max(data.OOS3()[[1]],na.rm = T)))
  
  data.Var <- reactive(Choice(data.Stock(),"4-Daily"))
  Variance <- reactive(sd(diff(data.Var()[[1]][seq(1, length(data.Var()[[2]][,1]),3),2]),na.rm = T))
  # Trend Segmente
  Trend.OS <- reactive({
    Segmente.Trend(data()[[1]][,1],
                   data()[[2]][,1],data()[[1]][,2],
                   input$date.OOS,listOOS1()[[8]],input$OS.Perf,Variance(),input$Fee,input$Fee.num)
  })
  
  output$TrendsM_Out <- renderPlotly(Trend.OS()[[1]])
  
  # Cross-Correlation
  Cross.Correlation <- reactive({
    Cross.Corr(data()[[1]][-1,1],data()[[2]][,1],
               listOOS1()[[8]],input$date.OOS,input$OS.Cross,input$Freq)
  })
  
  output$Cross_Corr <- renderPlotly(Cross.Correlation()[[1]])
  output$Cross_Corr.Text <- renderPrint(Cross.Correlation()[[2]])
  
  # Hit-Ratio
  output$Performance_Chart <- renderPlotly(if(input$OS.Perf){Trend.OS()[[7]]} else {Trend.OS()[[6]]})
  output$OutPerformance_Chart <- renderPlotly(if(input$OS.Perf){Trend.OS()[[9]]} else {Trend.OS()[[8]]})
  
  ######################
  ### Aktuelle Resultate
  ######################
  
  # Daten mit Frequenzen
  data.D <- reactive(Choice(data.Stock(),"Daily"))
  data.W <- reactive(Choice(data.Stock(),"2-Daily"))
  data.2W <- reactive(Choice(data.Stock(),"3-Daily"))
  data.M <- reactive(Choice(data.Stock(),"4-Daily"))
  
  ############
  # Univariate
  ############
  
  # Daily
  data.1.D <- reactive({Univariate.Filter.cust(data.D()[[2]][,1],0,LangeT.1,1/trendT.1,lambdaT.1,etaT.1,pi/trendT.1,ConstraintT.1)[[1]]})
  
  listD <- reactive(if(!input$DF){data.1.D()} else {
    Univariate.Filter.cust.conv(data.D()[[2]][,1],0,LangeT.1,1/trendT.1,lambdaT.1,etaT.1,pi/trendT.1,ConstraintT.1,
                                0,LangeT.2,1/trendT.2,lambdaT.2,etaT.2,pi/trendT.2,ConstraintT.2)[[1]]})
  
  # Weekly
  data.1.W <- reactive({Univariate.Filter.cust(data.W()[[2]][,1],0,LangeW.1,1/trendW.1,lambdaW.1,etaW.1,pi/trendW.1,ConstraintW.1)[[1]]})
  
  listW <- reactive(if(!input$DF){data.1.W()} else {
    Univariate.Filter.cust.conv(data.W()[[2]][,1],0,LangeW.1,1/trendW.1,lambdaW.1,etaW.1,pi/trendW.1,ConstraintW.1,
                                0,LangeW.2,1/trendW.2,lambdaW.2,etaW.2,pi/trendW.2,ConstraintW.2)[[1]]})
  
  # 2-Weekly
  data.1.2W <- reactive({Univariate.Filter.cust(data.2W()[[2]][,1],0,Lange2W.1,1/trend2W.1,lambda2W.1,eta2W.1,pi/trend2W.1,Constraint2W.1)[[1]]})
  
  list2W <- reactive(if(!input$DF){data.1.2W()} else {
    Univariate.Filter.cust.conv(data.2W()[[2]][,1],0,Lange2W.1,1/trend2W.1,lambda2W.1,eta2W.1,pi/trend2W.1,Constraint2W.1,
                                0,Lange2W.2,1/trend2W.2,lambda2W.2,eta2W.2,pi/trend2W.2,Constraint2W.2)[[1]]})
  
  # Monthly
  data.1.M <- reactive({Univariate.Filter.cust(data.M()[[2]][,1],0,LangeM.1,1/trendM.1,lambdaM.1,etaM.1,pi/trendM.1,ConstraintM.1)[[1]]})
  
  listM <- reactive(if(!input$DF){data.1.M()} else {
    Univariate.Filter.cust.conv(data.M()[[2]][,1],0,LangeM.1,1/trendM.1,lambdaM.1,etaM.1,pi/trendM.1,ConstraintM.1,
                                0,LangeM.2,1/trendM.2,lambdaM.2,etaM.2,pi/trendM.2,ConstraintM.2)[[1]]})
  
  # Durchschnitt
  Mean.Vect <- reactive(Mean.of.Vectors(data.D()[[1]][-1,1],data.W()[[1]][-1,1],
                                        data.2W()[[1]][-1,1],data.M()[[1]][-1,1],
                                        Skalieren(listD()),Skalieren(listW()),
                                        Skalieren(list2W()),Skalieren(listM())))
  
  # Filter-Plot
  output$Filter_Combo <- renderPlotly(if(input$Zoom){
    plot_ly(x = data.D()[[1]][-1,1][which(as.character(data.D()[[1]][-1,1])>=Sys.Date()-730)], 
            y = Skalieren(listD())[which(as.character(data.D()[[1]][-1,1])>=Sys.Date()-730)],
            type = "scatter", mode = "lines", name = "Täglich") %>%
      layout(yaxis = list(
        title = "",
        showticklabels = FALSE
      )) %>%
      add_lines(x = data.W()[[1]][-1,1][which(as.character(data.W()[[1]][-1,1])>=Sys.Date()-730)],
                y = Skalieren(listW())[which(as.character(data.W()[[1]][-1,1])>=Sys.Date()-730)], name = "Wöchentlich") %>%
      add_lines(x = data.2W()[[1]][-1,1][which(as.character(data.2W()[[1]][-1,1])>=Sys.Date()-730)],
                y = Skalieren(list2W())[which(as.character(data.2W()[[1]][-1,1])>=Sys.Date()-730)], name = "2-Wöchentlich") %>%
      add_lines(x = data.M()[[1]][-1,1][which(as.character(data.M()[[1]][-1,1])>=Sys.Date()-730)],
                y = Skalieren(listM())[which(as.character(data.M()[[1]][-1,1])>=Sys.Date()-730)], name = "Monatlich") %>%
      add_lines(x = data.D()[[1]][-1,1][which(as.character(data.D()[[1]][-1,1])>=Sys.Date()-730)],
                y = Mean.Vect()[[1]][which(as.character(data.D()[[1]][-1,1])>=Sys.Date()-730)], name = "Durchschnitt") %>%
      add_lines(x = data.D()[[1]][-1,1][which(as.character(data.D()[[1]][-1,1])>=Sys.Date()-730)],
                y = Skalieren(data.D()[[2]][,1],3)[which(as.character(data.D()[[1]][-1,1])>=Sys.Date()-730)], name = "Underlying Data",opacity = 0.3)
  } else {
    plot_ly(x = data.D()[[1]][-1,1], y = Skalieren(listD()),
            type = "scatter", mode = "lines", name = "Täglich") %>%
      layout(yaxis = list(
        title = "",
        showticklabels = FALSE
      )) %>%
      add_lines(x = data.W()[[1]][-1,1],y = Skalieren(listW()), name = "Wöchentlich") %>%
      add_lines(x = data.2W()[[1]][-1,1],y = Skalieren(list2W()), name = "2-Wöchentlich") %>%
      add_lines(x = data.M()[[1]][-1,1],y = Skalieren(listM()), name = "Monatlich") %>%
      add_lines(x = data.D()[[1]][-1,1],
                y = Mean.Vect()[[1]], name = "Durchschnitt") %>%
      add_lines(x = data.D()[[1]][-1,1],
                y = Skalieren(data.D()[[2]][,1],3), name = "Underlying Data",opacity = 0.3)})
  
  # Segmente Plot
  output$Trend_Durchschnitt <- renderPlotly(Only.Segmente(data.D()[[1]][,1],Mean.Vect(),data.D()[[1]][,2],input$Freq1))
  
  # Cross Correlation
  CCAR <- reactive(Cross.Corr.akt.R(data.D()[[2]][,1],Mean.Vect()[[2]],Mean.Vect()[[3]],Mean.Vect()[[4]],Mean.Vect()[[5]],Mean.Vect()[[1]]))
  
  output$Cross_Corr_AR <- renderPlotly(CCAR()[[1]])
  output$Cross_Corr.AR <- renderTable(CCAR()[[2]],rownames = T)
  
  # Tabelle
  output$Resultate <- renderTable(Result(data.D()[[2]][,1],Mean.Vect()[[2]],Mean.Vect()[[3]],Mean.Vect()[[4]],Mean.Vect()[[5]],Mean.Vect()[[1]],data.D()[[1]][,1]),rownames = T)
  
  
}
shinyApp(ui = ui, server = server)