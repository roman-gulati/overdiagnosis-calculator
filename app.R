##################################################
# Shiny calculator to explore disease incidence
# and overdiagnosis in trial and population settings
##################################################
library(plyr)
library(reshape)
library(overdiag)
library(ggplot2)
library(shiny)
library(shinydashboard)

Header <- dashboardHeader(title='Overdiagnosis calculator',
                          titleWidth=250)

Sidebar <- dashboardSidebar(width=250,
             sidebarMenu(
               menuItem('Home', tabName='home'),
               menuItem('Trial setting', tabName='trial'),
               menuItem('Population setting', tabName='population'),
               menuItem('Documentation', tabName='documentation')))

Body <- dashboardBody(
           tags$head(tags$style(HTML('
                                     .skin-green .main-sidebar {
                                         background-color: #ffffff;
                                     }
                                     .skin-green .main-sidebar .sidebar .sidebar-menu .active a{
                                         background-color: #008d4c;
                                     }
                                     .skin-green .main-sidebar .sidebar .sidebar-menu a{
                                         background-color: #00a65a;
                                     }
                                     .content-wrapper,
                                     .right-side {
                                         background-color: #ffffff;
                                     }'))),
           tabItems(
           # Home page
           tabItem(tabName='home',
                   fluidRow(
                     box(width=12,
                         solidHeader=TRUE,
                         h1('Welcome!'),
                         br(),
                         h4('The overdiagnosis calculator provides an interface to a deterministic model of the effects of screening on disease incidence in randomized trials and population studies.'),
                         br(),
                         h4('To get started, select the', strong('Trial setting'), 'or', strong('Population setting'), 'from the tabs on the left. There you can specify input parameters that control multiple aspects of the chosen disease setting and the effects of screening. In either setting, a figure illustrates disease incidence and responds to the input parameters in real time. When available, the first time point at which empirical estimates are unbiased will be shown in the figure. Additional information is available under the', strong('Documentation'), 'tab.'),
                         hr(),
                         h5('Questions or comments? Please email', a(href='mailto:rgulati@fredhutch.org', 'rgulati@fredhutch.org'))))),
           # Trial setting
           tabItem(tabName='trial',
                   fluidRow(
                            box(width=12,
                                h2('Trial setting'),
                                solidHeader=TRUE),
                            column(width=4,
                                   box(width=NULL,
                                       title='Input parameters',
                                       solidHeader=TRUE,
                                       sliderInput('arm.size',
                                                   'Arm size:',
                                                   min=1000,
                                                   max=100000,
                                                   value=10000,
                                                   step=1000),
                                       sliderInput('onset.rate',
                                                   'Annual rate of onset:',
                                                   min=0.0001,
                                                   max=0.1,
                                                   value=0.001,
                                                   step=0.0001),
                                       sliderInput('sojourn.time',
                                                   'Range of preclinical durations:',
                                                   min=0,
                                                   max=20,
                                                   value=c(0, 6),
                                                   step=1),
                                       sliderInput('sensitivity',
                                                   'Episode sensitivity:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.5,
                                                   step=0.01),
                                       sliderInput('overdiag.rate',
                                                   'Overdiagnosis rate:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.25,
                                                   step=0.01),
                                       sliderInput('stop.year',
                                                   'Stop screen year:',
                                                   min=2,
                                                   max=13,
                                                   value=13,
                                                   step=1))),
                            column(width=8,
                                   box(width=NULL,
                                       title='Disease incidence',
                                       solidHeader=TRUE,
                                       plotOutput('trial.plot'))))),
           # Population setting
           tabItem(tabName='population',
                   fluidRow(
                            box(width=12,
                                h2('Population setting'),
                                solidHeader=TRUE),
                            column(width=4,
                                   box(width=NULL,
                                       title='Input parameters',
                                       solidHeader=TRUE,
                                       sliderInput('pop.size',
                                                   'Population size:',
                                                   min=1000,
                                                   max=100000,
                                                   value=10000,
                                                   step=1000),
                                       sliderInput('onset.rate',
                                                   'Annual rate of onset:',
                                                   min=0.0001,
                                                   max=0.1,
                                                   value=0.001,
                                                   step=0.0001),
                                       sliderInput('sojourn.time',
                                                   'Range of preclinical durations:',
                                                   min=0,
                                                   max=20,
                                                   value=c(0, 6),
                                                   step=1),
                                       sliderInput('sensitivity',
                                                   'Episode sensitivity:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.5,
                                                   step=0.01),
                                       sliderInput('overdiag.rate',
                                                   'Overdiagnosis rate:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.25,
                                                   step=0.01),
                                       textInput('proportion',
                                                 'Proportion receiving screening:',
                                                 value='0.05, 0.1, 0.15'),
                                       textInput('start.year',
                                                 'Corresponding start year:',
                                                 value='2, 3, 4'))),
                            column(width=8,
                                   box(width=NULL,
                                       title='Disease incidence',
                                       solidHeader=TRUE,
                                       plotOutput('pop.plot'))))),
           # Documentation
           tabItem(tabName='documentation',
                   fluidRow(
                     box(width=12,
                         h2('Documentation'),
                         solidHeader=TRUE)))
                      ))

ui <- dashboardPage(Header, Sidebar, Body, skin='green')

server <- function(input, output){
    followup.years <- 30
    output$trial.plot <- renderPlot({
        tset <- trial_setting(input$pop.size,
                              onset.rate=input$onset.rate,
                              sojourn.min=min(input$sojourn.time),
                              sojourn.max=max(input$sojourn.time),
                              sensitivity=input$sensitivity,
                              overdiag.rate=input$overdiag.rate,
                              screen.stop.year=input$stop.year,
                              followup.years=followup.years)
        print(trial_plot(tset))
    })
    output$pop.plot <- renderPlot({
        proportion <- as.numeric(unlist(strsplit(input$proportion, ',')))
        start.year <- as.numeric(unlist(strsplit(input$start.year, ',')))
        if(sum(proportion) < 1){
            proportion <- c(proportion, 1-sum(proportion))
            start.year <- c(start.year, followup.years-2)
        }
        mpset <- multipopulation_setting(input$pop.size,
                                    onset.rate=input$onset.rate,
                                    sojourn.min=min(input$sojourn.time),
                                    sojourn.max=max(input$sojourn.time),
                                    sensitivity=input$sensitivity,
                                    overdiag.rate=input$overdiag.rate,
                                    proportion=proportion,
                                    start.year=start.year,
                                    followup.years=followup.years)
        print(multipopulation_plot(mpset))
    })
}

shinyApp(ui, server)

