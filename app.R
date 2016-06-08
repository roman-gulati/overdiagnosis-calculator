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
                                  'equal numbers of individuals are randomized to a screen or control arm. Latent disease onset occurs at a constant rate; the control arm receives no screen tests, so that control arm incidence matches the rate of onset. Sojourn time, or time between onset and diagnosis, follows a discrete uniform distribution. The screen arm begins screening in year 1. If the screen arm stops screening, the empirical difference between', em('cumulative incidence'), 'in the screen and control arm provides an unbiased estimate of overdiagnosis after screening stabilizes plus the maximum sojourn time. If the screen arm continues screening, the empirical difference between', em('annual incidence'), 'provides an unbiased estimate of overdiagnosis after screening stabilizes plus the maximum sojourn time.')),
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
                                       plotOutput('trial.plot'))))),
           # Population setting
           tabItem(tabName='population',
                   fluidRow(
                            box(width=12,
                                solidHeader=TRUE,
                                p('In the',
                                  strong('Population setting'),
                                  'the population is infinite and latent disease onset occurs at a constant rate. Initially, before screening starts, the model projects disease incidence in steady state, so that diagnosis without screening matches the rate of onset. Sojourn time, or time between onset and diagnosis, follows a discrete uniform distribution. Annual screening begins in segments of the population at specified starting years. The empirical difference between annual incidence with and without screening provides an unbiased estimate of overdiagnosis after screening stabilizes plus the maximum sojourn time.')),
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
                                       plotOutput('pop.plot'))))),
           # Help
           tabItem(tabName='help',
                   fluidRow(
                     box(width=12,
                         solidHeader=TRUE,
                         h4('Conceptual details of the model are described here:'),
                         p(em('Gulati R, Feuer EJ, Etzioni R. Conditions for unbiased empirical estimates of cancer overdiagnosis in randomized trials and population studies. Am J Epidemiol, in press.')),
                         h4('An R package implementation of the model is available here:'),
                         a(href='http://github.com/roman-gulati/overdiag', 'http://github.com/roman-gulati/overdiag'))))
                      ))

ui <- dashboardPage(Header, Sidebar, Body, skin='green')

server <- function(input, output){
    followup.years <- 30
    output$trial.plot <- renderPlot({
        tset <- trial_setting(arm.size=input$arm.size,
                              onset.rate=input$trial.onset.rate,
                              sojourn.min=min(input$trial.sojourn.time),
                              sojourn.max=max(input$trial.sojourn.time),
                              sensitivity=input$trial.sensitivity,
                              overdiag.rate=input$trial.overdiag.rate,
                              screen.stop.year=input$trial.stop.year,
                              followup.years=followup.years)
        print(trial_plot(tset))
    })
    output$pop.plot <- renderPlot({
        proportion <- as.numeric(unlist(strsplit(input$pop.proportion, ',')))
        start.year <- as.numeric(unlist(strsplit(input$pop.start.year, ',')))
        if(sum(proportion) < 1){
            proportion <- c(proportion, 1-sum(proportion))
            start.year <- c(start.year, followup.years-2)
        }
        mpset <- multipopulation_setting(pop.size=input$pop.size,
                                    onset.rate=input$pop.onset.rate,
                                    sojourn.min=min(input$pop.sojourn.time),
                                    sojourn.max=max(input$pop.sojourn.time),
                                    sensitivity=input$pop.sensitivity,
                                    overdiag.rate=input$pop.overdiag.rate,
                                    proportion=proportion,
                                    start.year=start.year,
                                    followup.years=followup.years)
        print(multipopulation_plot(mpset))
    })
}

shinyApp(ui, server)

