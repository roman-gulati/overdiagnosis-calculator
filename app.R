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

Header <- dashboardHeader(title='Overdiagnosis explorer',
                          titleWidth=250)

Sidebar <- dashboardSidebar(width=250,
             sidebarMenu(
               menuItem('Home', tabName='home'),
               menuItem('Trial setting', tabName='trial'),
               menuItem('Population setting', tabName='population'),
               menuItem('Documentation', tabName='documentation')))

Body <- dashboardBody(tabItems(
           # Home page
           tabItem(tabName='home',
                   fluidRow(
                     box(width=12,
                         h1('Overdiagnosis explorer'),
                         br(),
                         h4('Welcome! The overdiagnosis explorer provides an interface to a model of the effects of screening on disease incidence in randomized trials and population studies. You can control multiple aspects of either setting and use the projections to evaluate bias in empirical estimates of overdiagnosis and to determine the first time point at which empirical estimates are unbiased.'),
                         br(),
                         h4('To get started, select the appropriate', strong('Trial setting'), 'or', strong('Population setting'), 'from the tabs on the left. Note that the', strong('Trial setting'), 'is divided into 3 design types:'),
                         p(strong('Stop screening:'), 'screen arm stops screening after a fixed number of rounds and control arm receives no screens.'),
                         p(strong('Continue screening:'), 'screen arm continues screening indefinitely and control arm receives no screens.'),
                         p(strong('Delay screening:'), 'the screen arm continues screening indefinitely and control arm starts screening after a fixed delay.'),
                         br(),
                         h4('Additional information is available under the', strong('Documentation'), 'tab. For other questions or comments, please email', a(href='mailto:rgulati@fredhutch.org', 'rgulati@fredhutch.org'))))),
           # Trial setting
           tabItem(tabName='trial',
                   fluidRow(
                            box(width=12,
                                h2('Trial setting')),
                            column(width=4,
                                   box(width=NULL,
                                       title='Inputs',
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
                                h2('Population setting')),
                            column(width=4,
                                   box(width=NULL,
                                       title='Inputs',
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
                         h2('Documentation'))))
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

print(shinyApp(ui, server))

