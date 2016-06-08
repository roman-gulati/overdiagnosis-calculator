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
               menuItem('Help', tabName='help')))

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
                         h4('The overdiagnosis calculator provides an interface to a deterministic model of the effects of screening on disease incidence in randomized trials and population studies.'),
                         h4('To get started, select the', strong('Trial setting'), 'or', strong('Population setting'), 'from the tabs on the left. There you can specify input parameters that control multiple aspects of the chosen disease setting and the effects of screening. In either setting, a figure illustrates disease incidence and reacts to changes to the input parameters in real time. When available, the first time point at which empirical estimates are unbiased will be shown in the figure. Additional information is available under the', strong('Help'), 'tab.'),
                         hr(),
                         h5('Questions or comments? Please email', a(href='mailto:rgulati@fredhutch.org', 'rgulati@fredhutch.org'))))),
           # Trial setting
           tabItem(tabName='trial',
                   fluidRow(
                            box(width=12,
                                solidHeader=TRUE,
                                p('In the',
                                  strong('Trial setting'),
                                  'equal numbers of individuals are randomized to a screen or control arm. Latent disease onset occurs at a constant rate. The preclinical duration follows a uniform distribution. The control arm receives no screen tests, so control arm incidence matches the rate of onset. The screen arm begins screening in year 1. The empirical difference between screen and control arm incidence can provide an unbiased estimate of the number of overdiagnoses cases under 2 conditions. (1) The difference is based on', em('cumulative incidence'), 'if the screen arm stops screening and on', em('annual incidence'), 'if the screen arm continues screening. (2) The difference is calculated after screening stabilizes plus the maximum preclinical duration.')),
                            column(width=4,
                                   box(width=NULL,
                                       title='Input parameters',
                                       solidHeader=TRUE,
                                       sliderInput('trial.size',
                                                   'Arm size:',
                                                   min=1000,
                                                   max=100000,
                                                   value=10000,
                                                   step=1000),
                                       sliderInput('trial.onset.rate',
                                                   'Annual rate of onset:',
                                                   min=0.0001,
                                                   max=0.1,
                                                   value=0.001,
                                                   step=0.0001),
                                       sliderInput('trial.sojourn.time',
                                                   'Range of preclinical durations:',
                                                   min=0,
                                                   max=20,
                                                   value=c(0, 6),
                                                   step=1),
                                       sliderInput('trial.sensitivity',
                                                   'Episode sensitivity:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.5,
                                                   step=0.01),
                                       sliderInput('trial.overdiag.rate',
                                                   'Overdiagnosis rate:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.25,
                                                   step=0.01),
                                       sliderInput('trial.stop.year',
                                                   'Stop screen year:',
                                                   min=2,
                                                   max=13,
                                                   value=13,
                                                   step=1))),
                            column(width=8,
                                   box(width=NULL,
                                       title='Disease incidence',
                                       solidHeader=TRUE,
                                       plotOutput('trial.plot')),
                                   box(width=NULL,
                                       title='Save projections',
                                       solidHeader=TRUE,
                                       downloadButton('trial.data', 'Download'))))),
           # Population setting
           tabItem(tabName='population',
                   fluidRow(
                            box(width=12,
                                solidHeader=TRUE,
                                p('In the',
                                  strong('Population setting'),
                                  'latent disease onset occurs at a constant rate and the preclinical duration follows a uniform distribution. Initially, before screening starts, the model projects disease incidence in steady state, so that diagnosis without screening matches the rate of latent onset. Annual screening begins in segments of the population at specified starting years. The empirical difference between annual incidence with and without screening provides an unbiased estimate of overdiagnosis after screening stabilizes plus the maximum preclinical duration.')),
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
                                       sliderInput('pop.onset.rate',
                                                   'Annual rate of onset:',
                                                   min=0.0001,
                                                   max=0.1,
                                                   value=0.001,
                                                   step=0.0001),
                                       sliderInput('pop.sojourn.time',
                                                   'Range of preclinical durations:',
                                                   min=0,
                                                   max=20,
                                                   value=c(0, 6),
                                                   step=1),
                                       sliderInput('pop.sensitivity',
                                                   'Episode sensitivity:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.5,
                                                   step=0.01),
                                       sliderInput('pop.overdiag.rate',
                                                   'Overdiagnosis rate:',
                                                   min=0.01,
                                                   max=1,
                                                   value=0.25,
                                                   step=0.01),
                                       textInput('pop.proportion',
                                                 'Proportion receiving screening:',
                                                 value='0.05, 0.1, 0.15'),
                                       textInput('pop.start.year',
                                                 'Corresponding start year:',
                                                 value='2, 3, 4'))),
                            column(width=8,
                                   box(width=NULL,
                                       title='Disease incidence',
                                       solidHeader=TRUE,
                                       plotOutput('pop.plot')),
                                   box(width=NULL,
                                       title='Save projections',
                                       solidHeader=TRUE,
                                       downloadButton('pop.data', 'Download'))))),
           # Help
           tabItem(tabName='help',
                   fluidRow(
                     box(width=12,
                         solidHeader=TRUE,
                         h4('Conceptual details of the model are described here:'),
                         p(em('Gulati R, Feuer EJ, Etzioni R. Conditions for unbiased empirical estimates of cancer overdiagnosis in randomized trials and population studies. Am J Epidemiol, in press.')),
                         h4('An R package implementation of the model is available here:'),
                         a(href='https://github.com/roman-gulati/overdiag', 'https://github.com/roman-gulati/overdiag'))))
                      ))

ui <- dashboardPage(Header, Sidebar, Body, skin='green')

server <- function(input, output){
    followup.years <- 30
    followup.output <- 20
    tset <- reactive({
        dset <- trial_setting(arm.size=input$trial.size,
                              onset.rate=input$trial.onset.rate,
                              sojourn.min=min(input$trial.sojourn.time),
                              sojourn.max=max(input$trial.sojourn.time),
                              sensitivity=input$trial.sensitivity,
                              overdiag.rate=input$trial.overdiag.rate,
                              screen.stop.year=input$trial.stop.year,
                              followup.years=followup.years)
        dset <- subset(dset, year <= followup.output)
    })
    mpset <- reactive({
        proportion <- as.numeric(unlist(strsplit(input$pop.proportion, ',')))
        start.year <- as.numeric(unlist(strsplit(input$pop.start.year, ',')))
        if(sum(proportion) < 1){
            proportion <- c(proportion, 1-sum(proportion))
            start.year <- c(start.year, followup.years-2)
        }
        dset <- multipopulation_setting(pop.size=input$pop.size,
                                        onset.rate=input$pop.onset.rate,
                                        sojourn.min=min(input$pop.sojourn.time),
                                        sojourn.max=max(input$pop.sojourn.time),
                                        sensitivity=input$pop.sensitivity,
                                        overdiag.rate=input$pop.overdiag.rate,
                                        proportion=proportion,
                                        start.year=start.year,
                                        followup.years=followup.years)
        dset <- subset(dset, year <= followup.output)
    })
    output$trial.plot <- renderPlot({print(trial_plot(tset()))})
    output$trial.data <- downloadHandler(filename <- 'trial_incidence.csv',
                                         content <- function(filename){
                                             write.csv(tset(),
                                                       file=filename,
                                                       quote=FALSE,
                                                       row.names=FALSE)
                                         },
                                         contentType='text/csv')
    output$pop.plot <- renderPlot({print(multipopulation_plot(mpset()))})
    output$pop.data <- downloadHandler(filename <- 'population_incidence.csv',
                                       content <- function(filename){
                                           write.csv(mpset(),
                                                     file=filename,
                                                     quote=FALSE,
                                                     row.names=FALSE)
                                       },
                                       contentType='text/csv')
}

shinyApp(ui, server)

