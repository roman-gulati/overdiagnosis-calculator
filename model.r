##################################################
# Illustrate effects of sojourn time, overdiagnosis
# and screening test sensitivity on hypothetical
# population and trial incidence patterns
##################################################

library(plyr)
library(reshape)
library(grid)
library(ggplot2)
library(scales)

basepath <- '..'
datapath <- file.path(basepath, 'data')
plotpath <- file.path(basepath, 'plots')

set.seed(12345)

##################################################
# Generate incidence in the absence of screening
# using fixed or uniform sojourn time distribution
##################################################

generate_absence <- function(population.size,
                             followup.years,
                             onset.rate,
                             sojourn.time.min,
                             sojourn.time.max){
    stopifnot(0 <= sojourn.time.min & sojourn.time.min <= sojourn.time.max)
    cat('sojourn time: [', sojourn.time.min, ',', sojourn.time.max, ']\n', sep='')
    if(sojourn.time.min == sojourn.time.max){
        sojourn.denom <- 1
        sojourn.time <- sojourn.time.min
    } else {
        sojourn.denom <- sojourn.time.max-sojourn.time.min+1
        sojourn.time <- seq(sojourn.time.min, sojourn.time.max)
    }
    probabilities <- rep(1/sojourn.denom, sojourn.denom)
    stopifnot(sum(probabilities) == 1)
    pset <- data.frame(onset_year=seq(-sojourn.time.max, followup.years),
                       count_onset=population.size*onset.rate)
    pset <- rbind(pset,
                  data.frame(onset_year=followup.years+seq(sojourn.time.max),
                             count_onset=0))
    pset <- ddply(pset,
                  .(onset_year),
                  function(x)
                      with(x, 
                           data.frame(onset_year=unique(onset_year),
                                      count_onset=probabilities*unique(count_onset),
                                      sojourn=sojourn.time)))
    pset <- transform(pset, clinical_year=onset_year+sojourn)
    return(pset)
}

##################################################
# Calculate clinical incidence. The number of false
# negatives are bounded by (1) the later of onset
# and the start of screening and (2) the earlier of
# clinical presentation and the end of screening.
# Cancers present clinically when all screening
# tests are false negatives. Imperfect attendance
# is represented by removing the expected number
# of cancers in men who attend a sensitive test
# scaled by the number of tests offered.
##################################################

calculate_clinical <- function(dset, screen.start.year, screen.stop.year, attendance){
    dset <- transform(dset, lower_year=pmax(screen.start.year, onset_year))
    dset <- transform(dset, upper_year=pmin(screen.stop.year, clinical_year))
    dset <- transform(dset, tests_offered=pmax(upper_year-lower_year, 0))
    #if(with(dset, any(sojourn > 0)) & attendance < 1)
    #    browser()
    dset <- transform(dset, count_clinical=count_onset*(attendance*(1-sensitivity)+(1-attendance))^tests_offered)
    dset <- subset(dset, select=-c(lower_year, upper_year))
    return(dset)
}

##################################################
# Calculate screen incidence. Create new variables
# to record year of each screening round (exclusive
# of stopping year). The number of false negatives
# are bounded
# (1) below by 0 and screen round - 1 and
# (2) above by number of tests - 1 and sojourn time - 1
# Check that incidence of onset matches diagnoses.
# Reshape dataset to indicate number of screen
# detections in each screen year.
##################################################

calculate_screen <- function(dset, screen.start.year, screen.stop.year, attendance){
    sojourn_time <- with(dset, unique(sojourn))
    for(screen_year in seq(screen.start.year, screen.stop.year-1))
        dset[[paste('count_screen', screen_year, sep='_')]] <- 0
    if(sojourn_time > 0){
        for(screen_year in seq(screen.start.year, screen.stop.year-1)){
            prev_years <- seq(screen_year-sojourn_time+1, screen_year)
            latent_index <- with(dset, which(onset_year %in% prev_years))
            latent_dset <- dset[latent_index, ]
            latent_times <- with(latent_dset, seq(sojourn_time-1, 0))
            latent_tests_offered <- with(latent_dset,
                                         pmin(latent_times, tests_offered-1))
            latent_tests_offered_and_valid <- with(latent_dset,
                                           pmin(latent_tests_offered,
                                                screen_year-screen.start.year))
            latent_tests_offered_and_valid <- with(latent_dset,
                                       pmax(latent_tests_offered_and_valid, 0))
            prob_onetp <- with(latent_dset,
                               sensitivity*attendance*(attendance*(1-sensitivity)+(1-attendance))^latent_tests_offered_and_valid)
            cases <- with(latent_dset, prob_onetp*count_onset)
            screen_index <- paste('count_screen', screen_year, sep='_')
            dset[latent_index, screen_index] <- cases
        }
    }
    count_all <- grep('count_screen|count_clinical', names(dset), value=TRUE)
    #if(sojourn_time > 0) browser()
    #if(any(rowSums(dset[, count_all]) != dset[, 'count_onset']))
    #    browser()
    #index <- which(rowSums(dset[, count_all]) != dset[, 'count_onset'])
    #r <- rowSums(dset[, count_all])
    #o <- dset[, 'count_onset']
    stopifnot(isTRUE(all.equal(rowSums(dset[, count_all]), dset[, 'count_onset'])))
    count_screen <- grep('count_screen', names(dset), value=TRUE)
    melted <- melt(dset, measure.vars=count_screen)
    melted <- rename(melted, c('variable'='screen_year', 'value'='count_screen'))
    melted <- transform(melted, screen_year=as.integer(sub('count_screen_', '', screen_year)))
    return(melted)
}

##################################################
# Generate incidence in the presence of screening
##################################################

generate_presence <- function(pset,
                              followup.years,
                              screen.start.year,
                              screen.stop.year,
                              sensitivity,
                              sojourn.time.min,
                              sojourn.time.max,
                              attendance=1){
    stopifnot(0 <= sensitivity & sensitivity <= 1)
    stopifnot(0 <= screen.start.year & screen.start.year <= screen.stop.year)
    stopifnot(screen.stop.year <= followup.years)
    pset <- transform(pset, sensitivity=sensitivity)
    pset <- ddply(pset,
                  .(sojourn),
                  calculate_clinical,
                  screen.start.year=screen.start.year,
                  screen.stop.year=screen.stop.year,
                  attendance=attendance)
    pset <- ddply(pset,
                  .(sojourn),
                  calculate_screen,
                  screen.start.year=screen.start.year,
                  screen.stop.year=screen.stop.year,
                  attendance=attendance)
    pset <- subset(pset, select=-tests_offered)
    return(pset)
}

##################################################
# Generate overdiagnosis
##################################################

generate_overdiag <- function(pset, overdiag.rate){
    stopifnot(0 <= overdiag.rate & overdiag.rate <= 1)
    pset <- transform(pset, count_overdiag=overdiag.rate*count_screen/(1-overdiag.rate))
    return(pset)
}

##################################################
# Control high-level graphics theme
##################################################

gg_theme <- function(...){
    theme_set(theme_bw())
    theme_update(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.border=element_blank(),
                 panel.margin=unit(0.025, 'npc'),
                 axis.ticks=element_line(size=0.25),
                 axis.text=element_text(size=14),
                 axis.title=element_text(angle=0, size=22),
                 axis.line=element_line(colour='black', size=0.25),
                 strip.background=element_rect(fill=NA, colour=NA),
                 strip.text=element_text(angle=0, size=14),
                 legend.position='none',
                 ...)
}

##################################################
# Examine population setting
##################################################

population_setting <- function(dset,
                               population.size=1e5,
                               screen.start.year=4,
                               screen.stop.year=30,
                               onset.rate=0.001,
                               followup.years=30,
                               sojourn.distribution=FALSE,
                               verbose=FALSE){
    with(dset, cat(paste(rep('-', 40), collapse=''),
                   '\nsensitivity:', as.character(unique(sensitivity)),
                   '\noverdiag.rate:', as.character(unique(overdiag.rate)),
                   '\n'))
    # generate population of population.size individuals and
    # record year of clinical diagnosis for batches of relevant
    # cancers that develop in each year with a given sojourn time
    pset <- with(dset, generate_absence(population.size,
                                        followup.years,
                                        onset.rate,
                                        sojourn.time.min,
                                        sojourn.time.max))
    # screen the population under assumed sensitivity by
    # counting screen diagnoses in each year of screening
    # for batches of relevant cancers that develop in each
    # year with a given sojourn time
    pset <- with(dset, generate_presence(pset,
                                         followup.years,
                                         screen.start.year,
                                         screen.stop.year,
                                         sensitivity,
                                         sojourn.time.min,
                                         sojourn.time.max))
    # append overdiagnoses as a constant fraction of screen
    # diagnoses in each year of screening
    pset <- with(dset, generate_overdiag(pset, overdiag.rate))
    # isolate clinical diagnoses by year and sojourn time
    clinical_sojourn <- ddply(pset,
                              .(clinical_year, sojourn),
                              summarize,
                              count_clinical=unique(count_clinical))
    # isolate screen diagnoses and overdiagnoses and sojourn time
    screened_sojourn <- ddply(pset,
                              .(screen_year, sojourn),
                              summarize,
                              count_screen=sum(count_screen),
                              count_overdiag=sum(count_overdiag))
    if(sojourn.distribution){
        merged_sojourn <- merge(rename(screened_sojourn, c('screen_year'='year')),
                                rename(clinical_sojourn, c('clinical_year'='year')),
                                by=c('year', 'sojourn'),
                                all=TRUE)
        merged_sojourn <- ddply(merged_sojourn,
                                .(year, sojourn),
                                summarize,
                                count_relevant=sum(count_screen, na.rm=TRUE)+unique(count_clinical))
        return(merged_sojourn)
    }
    # count clinical diagnoses in each year
    clinical <- ddply(clinical_sojourn,
                      .(clinical_year),
                      summarize,
                      count_clinical=sum(count_clinical))
    # count screen diagnoses and overdiagnoses in each year
    screened <- ddply(screened_sojourn,
                      .(screen_year),
                      summarize,
                      count_screen=sum(count_screen),
                      count_overdiag=sum(count_overdiag))
    if(verbose)
        if(with(dset, sojourn.time.min != sojourn.time.max & overdiag.rate == 0)){
            leadtime <- with(dset,
                             sum(sapply(seq(sojourn.time.max-1),
                                        function(x)
                                            sensitivity*(1-sensitivity)^(sojourn.time.max-x-1)*x)))
            cat('sensitivity:', with(dset, sensitivity), '\n')
            cat('lead time:', round(leadtime, 2), '\n')
        }
    # pad screen diagnoses and overdiagnoses with 0
    # in all pre-screening years
    screened <- rbind(data.frame(screen_year=seq(0, screen.start.year-1),
                                 count_screen=0,
                                 count_overdiag=0),
                      screened)
    # merge screen diagnoses, overdiagnoses, and clinical diagnoses
    # from pre-screening through all years of screening
    merged <- merge(rename(screened, c('screen_year'='year')),
                    rename(clinical, c('clinical_year'='year')))
    return(merged)
}

population_incidence <- function(population.size=1e5, screen.start.year=4){
    pset <- data.frame(sojourn.time.min=c(2, 0, 4, 2),
                       sojourn.time.max=c(2, 4, 4, 6))
    pset <- rbind(transform(pset, overdiag.rate=0),
                  transform(pset, overdiag.rate=0.25))
    pset <- rbind(transform(pset, sensitivity=1),
                  transform(pset, sensitivity=0.5))
    pset <- ddply(pset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    sensitivity,
                    overdiag.rate),
                  population_setting,
                  population.size=population.size,
                  screen.start.year=screen.start.year)
    return(pset)
}
#pset <- population_incidence()

##################################################
# Implement dissemination: call population_incidence
# for specified dissemination pattern then calculate
# incidence rates across years
##################################################

multipopulation_incidence <- function(population.size=1e5,
                                      span=5,
                                      proportion=rep(1/span, length=span),
                                      start.year=seq(4, length=span)){
    dissemination <- data.frame(population.size=population.size,
                                proportion=proportion,
                                start.year=start.year)
    stopifnot(with(dissemination, sum(proportion)) == 1)
    mpset <- ddply(dissemination,
                   .(proportion, start.year),
                   function(x){
                       cat('proportion:', with(x, proportion), '\n')
                       cat('start.year:', with(x, start.year), '\n')
                       with(x,
                            population_incidence(population.size=proportion*population.size,
                                                 screen.start.year=start.year))
                   })
    mpset <- ddply(mpset,
                   .(sojourn.time.min, sojourn.time.max, sensitivity, overdiag.rate, year),
                   summarize,
                   count_screen=sum(count_screen),
                   count_clinical=sum(count_clinical),
                   count_overdiag=sum(count_overdiag))
    return(mpset)
}
#pset <- multipopulation_incidence(span=1)
#mpset_fast <- multipopulation_incidence(span=2)
#mpset_slow <- multipopulation_incidence(span=4)

##################################################
# Illustrate population setting
##################################################

population_plot <- function(dset,
                            Sensitivity=1,
                            Overdiagnosis=0,
                            Dissemination=NA,
                            minyear=0,
                            maxyear=12,
                            minyear.text=4,
                            line.size=0.1,
                            text.size=3,
                            text.offset=20,
                            saveit=FALSE){
    dset <- droplevels(subset(dset,
                              minyear <= year & year <= maxyear &
                              sensitivity == Sensitivity &
                              overdiag.rate == Overdiagnosis,
                              select=-c(sensitivity, overdiag.rate)))
    dset <- transform(dset, sojourn.time=paste0('[',
                                                sojourn.time.min,
                                                ', ',
                                                sojourn.time.max,
                                                ']'))
    dset <- ddply(dset,
                  .(sojourn.time),
                  transform,
                  sojourn.mean=paste('Mean:',
                                     mean(seq(unique(sojourn.time.min),
                                              unique(sojourn.time.max))),
                                     'years'),
                  sojourn.range=paste('Range:\n',
                                      unique(sojourn.time.max-sojourn.time.min),
                                      'years'))
    dset <- transform(dset,
                      sojourn.time=factor(sojourn.time,
                                          levels=c('[2, 2]',
                                                   '[0, 4]',
                                                   '[4, 4]',
                                                   '[2, 6]')),
                      sojourn.mean=factor(sojourn.mean),
                      sojourn.range=relevel(factor(sojourn.range),
                                            ref=c('Range:\n 0 years')))
    dset.text <- subset(dset, year %in% seq(minyear.text, maxyear-1))
    dset.text <- transform(dset.text,
                           text_size=text.size,
                           text_offset=text.offset)
    gg_theme()
    gg <- ggplot(dset)
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=0,
                             ymax=count_clinical),
                         #fill='gray90')
                         fill=alpha('skyblue', 0.75))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical,
                             ymax=count_clinical+count_screen),
                         #fill='gray60')
                         fill=alpha('orange', 0.75))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical+count_screen,
                             ymax=count_clinical+count_screen+count_overdiag),
                         #fill='gray30')
                         fill=alpha('seagreen', 0.75))
    gg <- gg+geom_line(aes(x=year, y=count_clinical), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen+count_overdiag), size=line.size)
    gg <- gg+geom_hline(yintercept=0, size=0.25)
    gg <- gg+geom_vline(xintercept=minyear, size=0.25)
    gg <- gg+facet_grid(sojourn.range~sojourn.mean)
    gg <- gg+scale_x_continuous(name='\nYear',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Incidence rate per 100,000\n',
                                limits=c(0, 550),
                                expand=c(0, 0))
    if(Overdiagnosis == 0){
        gg <- gg+geom_text(data=dset.text,
                           aes(x=year,
                               y=count_clinical+count_screen+text_offset,
                               label=paste0('+', round(count_screen)),
                               size=text_size))
    } else {
        gg <- gg+geom_text(data=dset.text,
                           aes(x=year,
                               y=count_clinical+count_screen-text_offset,
                               label=paste0('+', round(count_screen)),
                               size=text_size))
        gg <- gg+geom_text(data=dset.text,
                           aes(x=year,
                               y=count_clinical+count_screen+count_overdiag+text_offset,
                               label=paste0('+', round(count_overdiag)),
                               size=text_size))
    }
    if(nrow(subset(dset.text, count_clinical <= text_offset)) > 0)
        gg <- gg+geom_text(data=subset(dset.text, count_clinical <= text_offset),
                           aes(x=year,
                               y=count_clinical+text_offset,
                               label=paste0('-', 100-round(count_clinical)),
                               size=text_size))
    gg <- gg+geom_text(data=subset(dset.text, count_clinical > text_offset),
                       aes(x=year,
                           y=count_clinical-text_offset,
                           label=paste0('-', 100-round(count_clinical)),
                           size=text_size))
    #print(gg)
    if(saveit){
        filename <- paste('population_prospective',
                          round(100*Sensitivity),
                          round(100*Overdiagnosis),
                          sep='_')
        if(!is.na(Dissemination))
            filename <- paste(filename, Dissemination, sep='_')
        filename <- paste(filename, 'pdf', sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=10,
               height=7)
    }
}
#population_plot(pset, Sensitivity=1, Overdiagnosis=0, saveit=TRUE)
#population_plot(pset, Sensitivity=0.5, Overdiagnosis=0, saveit=TRUE)
#population_plot(pset, Sensitivity=0.5, Overdiagnosis=0.25, saveit=TRUE)
#population_plot(mpset_fast,
#                Sensitivity=0.5,
#                Overdiagnosis=0.25,
#                Dissemination='fast',
#                saveit=TRUE)
#population_plot(mpset_slow,
#                Sensitivity=0.5,
#                Overdiagnosis=0.25,
#                Dissemination='slow',
#                saveit=TRUE)

##################################################
# Revised settings to show more realistic setting
# as base case and selected variants
##################################################

population_incidence_revised <- function(population.size=1e5,
                                         screen.start.year=4){
    pset <- data.frame(sojourn.time.min=c(0, 0),
                       sojourn.time.max=c(12, 6))
    pset <- rbind(transform(pset, overdiag.rate=0),
                  transform(pset, overdiag.rate=0.25))
    pset <- rbind(transform(pset, sensitivity=1),
                  transform(pset, sensitivity=0.5))
    pset <- ddply(pset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    sensitivity,
                    overdiag.rate),
                  population_setting,
                  population.size=population.size,
                  screen.start.year=screen.start.year)
    return(pset)
}

population_incidence_3ways <- function(population.size=1e5,
                                       screen.start.year=3,
                                       sojourn.distribution=FALSE){
    dset <- data.frame(sojourn.time.min=c(0, 0),
                       sojourn.time.max=c(12, 6),
                       sensitivity=0.5,
                       overdiag.rate=0.25)
    pset <- ddply(dset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    sensitivity,
                    overdiag.rate),
                  population_setting,
                  population.size=population.size,
                  screen.start.year=screen.start.year,
                  sojourn.distribution=sojourn.distribution)
    return(pset)
}

multipopulation_incidence_revised <- function(population.size=1e5,
                                              dissemination.pattern='default'){
    if(dissemination.pattern == 'default'){
        #proportion <- c(0.05, 0.1, 0.15, 0.15, 0.05, 0.5)
        #start.year <- c(4, 5, 6, 7, 8, 28)
        #dissemination.name <- 'Dissemination: 5%, 15%, 30%, 45%, 50%'
        proportion <- 1
        start.year <- 4
        dissemination.name <- 'Dissemination: 100%'
    } else {
        proportion <- c(0.1, 0.3, 0.1, 0.5)
        start.year <- c(4, 5, 6, 28)
        dissemination.name <- 'Dissemination: 10%, 40%, 50%'
    }
    dissemination <- data.frame(population.size=population.size,
                                proportion=proportion,
                                start.year=start.year)
    stopifnot(with(dissemination, sum(proportion)) == 1)
    mpset <- ddply(dissemination,
                   .(proportion, start.year),
                   function(x){
                       cat('proportion:', with(x, proportion), '\n')
                       cat('start.year:', with(x, start.year), '\n')
                       with(x,
        population_incidence_revised(population.size=proportion*population.size,
                             screen.start.year=start.year))
                   })
    mpset <- ddply(mpset,
                   .(sojourn.time.min,
                     sojourn.time.max,
                     sensitivity,
                     overdiag.rate,
                     year),
                   summarize,
                   count_screen=sum(count_screen),
                   count_clinical=sum(count_clinical),
                   count_overdiag=sum(count_overdiag))
    mpset <- transform(mpset, dissemination=dissemination.name)
    return(mpset)
}
#mpset_default <- multipopulation_incidence_revised(dissemination.pattern='default')
#mpset_variant <- multipopulation_incidence_revised(dissemination.pattern='variant')
#mpset <- rbind(mpset_default, mpset_variant)

multipopulation_incidence_3ways <- function(population.size=1e5,
                                            dissemination.pattern='default',
                                            sojourn.distribution=FALSE){
    if(dissemination.pattern == 'default'){
        proportion <- c(0.05, 0.1, 0.15, 0.15, 0.05, 0.5)
        start.year <- c(2, 3, 4, 5, 6, 28)
        dissemination.name <- 'Cumulative uptake (years 2-6): 5%,15%,30%,45%,50%'
    } else {
        proportion <- c(0.1, 0.3, 0.1, 0.5)
        start.year <- c(2, 3, 4, 28)
        dissemination.name <- 'Cumulative uptake (years 2-4): 10%,40%,50%'
    }
    dissemination <- data.frame(population.size=population.size,
                                proportion=proportion,
                                start.year=start.year)
    stopifnot(with(dissemination, sum(proportion)) == 1)
    mpset <- ddply(dissemination,
                   .(proportion, start.year),
                   function(x){
                       cat('proportion:', with(x, proportion), '\n')
                       cat('start.year:', with(x, start.year), '\n')
                       with(x,
        population_incidence_3ways(population.size=proportion*population.size,
                                   screen.start.year=start.year,
                                   sojourn.distribution=sojourn.distribution))
                   })
    if(sojourn.distribution){
        mpset <- ddply(mpset,
                       .(sojourn.time.min,
                         sojourn.time.max,
                         sensitivity,
                         overdiag.rate,
                         year,
                         sojourn),
                       summarize,
                       count_relevant=sum(count_relevant, na.rm=TRUE))
        mpset <- ddply(mpset,
                       .(sojourn.time.min,
                         sojourn.time.max,
                         sensitivity,
                         overdiag.rate,
                         year),
                       summarize,
                       sojourn=unique(sojourn),
                       count_relevant=count_relevant,
                       prop_relevant=count_relevant/sum(count_relevant, na.rm=TRUE))
    } else {
        mpset <- ddply(mpset,
                       .(sojourn.time.min,
                         sojourn.time.max,
                         sensitivity,
                         overdiag.rate,
                         year),
                       summarize,
                       count_screen=sum(count_screen),
                       count_clinical=sum(count_clinical),
                       count_overdiag=sum(count_overdiag))
    }
    mpset <- transform(mpset, dissemination=dissemination.name)
    return(mpset)
}
#mpset_default <- multipopulation_incidence_3ways(dissemination.pattern='default')
#mpset_variant <- multipopulation_incidence_3ways(dissemination.pattern='variant')
#mpset_3ways <- rbind(mpset_default, mpset_variant)

#msset_default <- multipopulation_incidence_3ways(dissemination.pattern='default',
#                                                 sojourn.distribution=TRUE)
#msset_variant <- multipopulation_incidence_3ways(dissemination.pattern='variant',
#                                                 sojourn.distribution=TRUE)
#msset_3ways <- rbind(msset_default, msset_variant)

##################################################
# Illustrate population setting
##################################################

population_plot_revised <- function(dset,
                                    Overdiagnosis=0,
                                    minyear=0,
                                    maxyear=22,
                                    minyear.text=4,
                                    line.size=0.1,
                                    text.size=3,
                                    text.offset=15,
                                    stable.arrows=FALSE,
                                    make.table=FALSE,
                                    saveit=FALSE){
    dset <- droplevels(subset(dset,
                              minyear <= year & year <= maxyear &
                              overdiag.rate == Overdiagnosis,
                              select=-overdiag.rate))
    dset <- transform(dset,
                      sensitivity=factor(paste0('Effective sensitivity: ',
                                                round(100*sensitivity),
                                                '%')),
                      sojourn.time=factor(paste('Max preclinical detectable period:',
                                                 sojourn.time.max)))
    dset <- transform(dset, label=apply(data.frame(sensitivity,
                                                   dissemination,
                                                   sojourn.time),
                                        1,
                                        paste,
                                        collapse='\n'))
    facets <- with(dset, as.data.frame(table(sensitivity,
                                             dissemination,
                                             sojourn.time)))
    # rows 1, 2, 3, 5 specify desired base case and variants
    #facets <- droplevels(subset(facets[c(5, 6, 7, 1), ], select=-Freq))
    facets <- droplevels(subset(facets[c(1, 2, 5, 6), ], select=-Freq))
    facets <- summarize(facets, label=apply(facets, 1, paste, collapse='\n'))
    dset <- merge(facets, dset, by='label')
    dset <- transform(dset, label=factor(label, levels=facets$label))
    dset <- sort_df(dset, vars=c('label', 'year'))
    dset <- ddply(dset,
                  .(label),
                  function(x){
                      baseline <- x[with(x, year == 0), 'count_clinical']
                      odx_count <- with(x, count_overdiag)
                      sdx_count <- with(x, count_screen)
                      cdx_count <- with(x, count_clinical)
                      year_offset <- sdx_count > 0
                      rough <- year_offset & odx_count/(cdx_count+sdx_count+odx_count-baseline) > 0.9
                      rough_year <- min(which(rough))
                      exact <- year_offset & abs(cdx_count+sdx_count-baseline) < 0.000001
                      exact_year <- min(which(exact))
                      print(with(x, as.character(unique(label))))
                      cat('Rough:', rough_year, '\n')
                      cat('Exact:', exact_year, '\n')
                      x <- transform(x,
                                     rough=year == rough_year,
                                     rough_tag='90%',
                                     exact=year == exact_year,
                                     exact_tag='100%')
                      return(x)
                  })
    if(make.table){
        ests <- ddply(dset,
                      .(label),
                      function(x){
                          baseline <- x[with(x, year == 0), 'count_clinical']
                          beg_year <- x[with(x, min(which(count_screen > 0))), 'year']
                          end_year <- x[with(x, exact), 'year']
                          xset <- subset(x, beg_year <= year & year <= end_year)
                          xset <- summarize(xset,
                                            followup=end_year,
                                            end_baseline=baseline,
                                            end_total=(count_clinical+count_screen+count_overdiag)[year == end_year],
                                            end_actual=count_overdiag[year == end_year],
                                            end_scases=count_screen[year == end_year],
                                            total_baseline=baseline*nrow(xset),
                                            total_incidence=sum(count_clinical+count_screen+count_overdiag),
                                            total_actual=sum(count_overdiag),
                                            total_excess=sum(count_clinical+count_screen+count_overdiag-baseline))
                          return(xset)
                      })
        ests <- transform(ests,
                          end_excess=end_total-end_baseline,
                          end_ratio=end_total/end_baseline,
                          abs_error=total_excess-total_actual,
                          rel_error=total_excess/total_actual)
        write.csv(ests,
                  file=file.path(datapath, 'population_estimates.csv'),
                  quote=FALSE,
                  row.names=FALSE)
    }
    dset.text <- subset(dset, year %in% seq(minyear.text, maxyear-1))
    dset.text <- transform(dset.text,
                           text_size=text.size,
                           text_offset=text.offset)
    gg_theme(panel.margin=unit(0.025, 'npc'),
             strip.text.x=element_text(size=16, hjust=0))
    gg <- ggplot(dset)
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=0,
                             ymax=count_clinical),
                         #fill='gray90')
                         fill=alpha('skyblue', 0.75))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical,
                             ymax=count_clinical+count_screen),
                         #fill='gray60')
                         fill=alpha('orange', 0.75))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical+count_screen,
                             ymax=count_clinical+count_screen+count_overdiag),
                         #fill='gray30')
                         fill=alpha('seagreen', 0.75))
    gg <- gg+geom_line(aes(x=year, y=count_clinical), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen+count_overdiag), size=line.size)
    gg <- gg+geom_hline(yintercept=0, size=0.25)
    gg <- gg+geom_vline(xintercept=minyear, size=0.25)
    gg <- gg+facet_wrap(~label, nrow=2)
    gg <- gg+scale_x_continuous(name='\nYear',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Incidence rate per 100,000\n',
                                #limits=c(0, 250),
                                #breaks=seq(0, 250, by=50),
                                expand=c(0, 0))
    if(Overdiagnosis == 0){
        gg <- gg+geom_text(data=dset.text,
                           aes(x=year,
                               y=count_clinical+count_screen+text_offset,
                               label=paste0('+', round(count_screen)),
                               size=text_size))
    } else {
        gg <- gg+geom_text(data=dset.text,
                           aes(x=year,
                               y=count_clinical+count_screen-text_offset,
                               label=paste0('+', round(count_screen)),
                               size=text_size))
        gg <- gg+geom_text(data=dset.text,
                           aes(x=year,
                               y=count_clinical+count_screen+count_overdiag+text_offset,
                               label=paste0('+', round(count_overdiag)),
                               size=text_size))
        if(stable.arrows){
            gg <- gg+geom_segment(data=subset(dset, exact),
                                  aes(x=year,
                                      xend=year,
                                      y=200,
                                      yend=150),
                                  colour='gray30',
                                  arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
            gg <- gg+geom_text(data=subset(dset, exact),
                               aes(x=year,
                                   y=220,
                                   label=exact_tag),
                               colour='gray30',
                               vjust=0,
                               size=4)
            gg <- gg+geom_segment(data=subset(dset, rough),
                                  aes(x=year,
                                      xend=year,
                                      y=200,
                                      yend=150),
                                  colour='gray60',
                                  arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
            gg <- gg+geom_text(data=subset(dset, rough),
                               aes(x=year,
                                   y=220,
                                   label=rough_tag),
                               colour='gray60',
                               vjust=0,
                               size=4)
        }
    }
    if(nrow(subset(dset.text, count_clinical <= text_offset)) > 0)
        gg <- gg+geom_text(data=subset(dset.text, count_clinical <= text_offset),
                           aes(x=year,
                               y=count_clinical+text_offset,
                               label=paste0('-', 100-round(count_clinical)),
                               size=text_size))
    gg <- gg+geom_text(data=subset(dset.text, count_clinical > text_offset),
                       aes(x=year,
                           y=count_clinical-text_offset,
                           label=paste0('-', 100-round(count_clinical)),
                           size=text_size))
    print(gg)
    if(saveit){
        filename <- paste('population_revised', round(100*Overdiagnosis), sep='_')
        if(stable.arrows)
            filename <- paste(filename, 'arrows', sep='_')
        filename <- paste(filename, '2015-02-11', sep='_')
        filename <- paste(filename, 'png', sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=12,
               height=7)
    }
}
#population_plot_revised(mpset, Overdiagnosis=0, saveit=TRUE)
#population_plot_revised(mpset, Overdiagnosis=0.25, saveit=TRUE)

#population_plot_revised(mpset,
#                        Overdiagnosis=0.25,
#                        stable.arrows=TRUE,
#                        make.table=FALSE,
#                        saveit=TRUE)

##################################################
# Illustrate population setting varying only
# dissemination and the max preclinical period
##################################################

population_plot_3ways <- function(dset,
                                  minyear=0,
                                  maxyear=20,
                                  minyear.text=2,
                                  line.size=0.1,
                                  text.size=3,
                                  text.offset=12,
                                  text.angle=30,
                                  make.table=FALSE,
                                  stable.arrows=FALSE,
                                  ext='png',
                                  saveit=FALSE){
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    dset <- transform(dset, sojourn.time=factor(paste('Max preclinical detectable period:',
                                                      sojourn.time.max,
                                                      'years')))
    dset <- transform(dset, type=apply(data.frame(dissemination, sojourn.time),
                                       1,
                                       paste,
                                       collapse='\n'))
    facets <- with(dset, as.data.frame(table(dissemination, sojourn.time)))
    facets <- droplevels(subset(facets[c(3, 4, 1), ], select=-Freq))
    facets <- summarize(facets, type=apply(facets, 1, paste, collapse='\n'))
    dset <- merge(facets, dset, by='type')
    dset <- transform(dset, type=factor(type, levels=facets$type))
    dset <- sort_df(dset, vars=c('type', 'year'))
    dset <- ddply(dset,
                  .(type),
                  function(x){
                      with(x, cat(paste(rep('-', 40), collapse=''),
                                  '\nsetting:', as.character(unique(type)),
                                  '\n'))
                      # extract background incidence
                      baseline <- x[with(x, year == 0), 'count_clinical']
                      # extract incidence under screening
                      overdiag_inc <- with(x, count_overdiag)
                      screen_inc <- with(x, count_screen)
                      clinical_inc <- with(x, count_clinical)
                      # calculate empirical excess
                      inc_total <- clinical_inc+screen_inc+overdiag_inc
                      excess_total <- inc_total-baseline
                      # specify valid years for evaluating empirical excess
                      year_offset <- screen_inc > 0
                      # calculate first year excess matches true overdiagnosis
                      exact <- year_offset & abs(overdiag_inc/excess_total-1) < 0.000001
                      rough <- year_offset & abs(overdiag_inc/excess_total-1) < 0.1
                      exact_year <- min(which(exact))
                      rough_year <- min(which(rough))
                      cat('exact:', with(x, year[exact_year]), '\n')
                      cat('rough:', with(x, year[rough_year]), '\n')
                      # append exact and rough approximation years to stratum
                      x <- transform(x,
                                     rough=year == year[rough_year],
                                     rough_tag='90%',
                                     exact=year == year[exact_year],
                                     exact_tag='Unbiased')
                      return(x)
                  })
    if(make.table){
        ests <- ddply(dset,
                      .(type),
                      function(x){
                          baseline <- x[with(x, year == 0), 'count_clinical']
                          beg_year <- x[with(x, min(which(count_screen > 0))), 'year']
                          end_year <- x[with(x, exact), 'year']
                          xset <- subset(x, beg_year <= year & year <= end_year)
                          xset <- summarize(xset,
                                            followup=end_year,
                                            end_baseline=baseline,
                                            end_total=(count_clinical+count_screen+count_overdiag)[year == end_year],
                                            end_actual=count_overdiag[year == end_year],
                                            end_scases=count_screen[year == end_year],
                                            total_baseline=baseline*nrow(xset),
                                            total_incidence=sum(count_clinical+count_screen+count_overdiag),
                                            total_actual=sum(count_overdiag),
                                            total_excess=sum(count_clinical+count_screen+count_overdiag-baseline))
                          return(xset)
                      })
        ests <- transform(ests,
                          end_excess=end_total-end_baseline,
                          end_ratio=end_total/end_baseline,
                          abs_error=total_excess-total_actual,
                          rel_error=total_excess/total_actual)
        write.csv(ests,
                  file=file.path(datapath, 'population_estimates.csv'),
                  quote=FALSE,
                  row.names=FALSE)
    }
    dset.text <- subset(dset, year %in% seq(minyear.text, maxyear-1))
    dset.text <- transform(dset.text,
                           text_size=text.size,
                           text_offset=text.offset,
                           text_angle=text.angle)
    gg_theme(panel.margin=unit(0.025, 'npc'),
             panel.border=element_rect(fill='transparent'),
             strip.text.x=element_text(size=12, hjust=0),
             axis.title=element_text(size=18))
    gg <- ggplot(dset)
    gg <- gg+geom_blank(aes(y=220))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=0,
                             ymax=count_clinical),
                         #fill='gray90')
                         fill=alpha('skyblue', 0.75))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical,
                             ymax=count_clinical+count_screen),
                         #fill='gray60')
                         fill=alpha('orange', 0.75))
    gg <- gg+geom_ribbon(aes(x=year,
                             ymin=count_clinical+count_screen,
                             ymax=count_clinical+count_screen+count_overdiag),
                         #fill='gray30')
                         fill=alpha('seagreen', 0.75))
    gg <- gg+geom_line(aes(x=year, y=count_clinical), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen), size=line.size)
    gg <- gg+geom_line(aes(x=year, y=count_clinical+count_screen+count_overdiag), size=line.size)
    gg <- gg+geom_hline(yintercept=0, size=0.25)
    gg <- gg+geom_vline(xintercept=minyear, size=0.25)
    gg <- gg+facet_grid(.~type)
    gg <- gg+scale_x_continuous(name='\nYear',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Annual incidence\nper 100,000 individuals\n',
                                #limits=c(0, 250),
                                #breaks=seq(0, 250, by=50),
                                expand=c(0, 0))
    gg <- gg+geom_text(data=dset.text,
                       aes(x=year,
                           y=count_clinical+count_screen-text_offset,
                           label=paste0('+', round(count_screen)),
                           size=text_size,
                           angle=text_angle))
    gg <- gg+geom_text(data=dset.text,
                       aes(x=year,
                           y=count_clinical+count_screen+count_overdiag+text_offset,
                           label=paste0('+', round(count_overdiag)),
                           size=text_size,
                           angle=text_angle))
    if(nrow(subset(dset.text, count_clinical <= text_offset)) > 0)
        gg <- gg+geom_text(data=subset(dset.text, count_clinical <= text_offset),
                           aes(x=year,
                               y=count_clinical+text_offset,
                               label=paste0('-', 100-round(count_clinical)),
                               size=text_size,
                               angle=text_angle))
    gg <- gg+geom_text(data=subset(dset.text, count_clinical > text_offset),
                       aes(x=year,
                           y=count_clinical-text_offset,
                           label=paste0('-', 100-round(count_clinical)),
                           size=text_size,
                           angle=text_angle))
    if(stable.arrows){
        gg <- gg+geom_segment(data=subset(dset, exact),
                              aes(x=year,
                                  xend=year,
                                  y=175,
                                  yend=140),
                              colour='gray30',
                              arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
        gg <- gg+geom_text(data=subset(dset, exact),
                           aes(x=year,
                               y=185,
                               label=exact_tag),
                           colour='gray30',
                           vjust=0,
                           size=4)
        #gg <- gg+geom_segment(data=subset(dset, rough),
        #                      aes(x=year,
        #                          xend=year,
        #                          y=175,
        #                          yend=140),
        #                      colour='gray60',
        #                      arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
        #gg <- gg+geom_text(data=subset(dset, rough),
        #                   aes(x=year,
        #                       y=185,
        #                       label=rough_tag),
        #                   colour='gray60',
        #                   vjust=0,
        #                   size=4)
    }
    print(gg)
    if(saveit){
        filename <- 'population_3ways'
        datestamp <- '2015-04-13'
        if(stable.arrows)
            filename <- paste(filename, 'arrows', sep='_')
        filename <- paste(filename, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=15,
               height=5)
    }
}
#population_plot_3ways(mpset_3ways,
#                      make.table=TRUE,
#                      stable.arrows=TRUE,
#                      ext='svg',
#                      saveit=FALSE)

##################################################
# Illustrate relevant sojourn time distribution
##################################################

sojourn_plot_3ways <- function(dset,
                               minyear=0,
                               maxyear=20,
                               minyear.text=3,
                               line.size=0.1,
                               text.size=3,
                               text.offset=15,
                               ext='png',
                               saveit=FALSE){
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    dset <- transform(dset, sojourn.time=factor(paste('Max preclinical detectable period:',
                                                      sojourn.time.max,
                                                      'years')))
    dset <- transform(dset, type=apply(data.frame(dissemination, sojourn.time),
                                       1,
                                       paste,
                                       collapse='\n'))
    facets <- with(dset, as.data.frame(table(dissemination, sojourn.time)))
    facets <- droplevels(subset(facets[c(3, 4, 1), ], select=-Freq))
    facets <- summarize(facets, type=apply(facets, 1, paste, collapse='\n'))
    dset <- merge(facets, dset, by='type')
    dset <- transform(dset, type=factor(type, levels=facets$type))
    dset <- sort_df(dset, vars=c('type', 'year', 'sojourn'))
    gg_theme(panel.margin=unit(0.025, 'npc'),
             strip.text.x=element_text(size=12, hjust=0),
             axis.title=element_text(size=18))
    gg <- ggplot(dset)
    gg <- gg+geom_bar(aes(x=year,
                          #y=prop_relevant,
                          y=count_relevant,
                          fill=log(sojourn+1)),
                      stat='identity',
                      position='stack')
    gg <- gg+facet_grid(.~type)
    gg <- gg+scale_x_continuous(name='\nYear',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear-0.5, maxyear+0.5),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(#name='Percent of diagnoses\nwith each preclinical period\n',
                                name='Annual incidence\nin each preclinical period\nper 100,000 individuals\n',
                                #limits=c(0, 1),
                                #breaks=seq(0, 1, by=1/7),
                                breaks=seq(0, 200, by=50),
                                #labels=percent_format(),
                                expand=c(0, 0))
    gg <- gg+scale_fill_gradient(low='brown', high='papayawhip')
    print(gg)
    if(saveit){
        filename <- 'sojourn_3ways'
        datestamp <- '2015-02-15'
        filename <- paste(filename, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=15,
               height=5)
    }
}
#sojourn_plot_3ways(msset_3ways, saveit=TRUE)

##################################################
# Examine trial setting
##################################################

trial_setting <- function(dset,
                          arm.size=1e5,
                          screen.start.year=1,
                          screen.stop.year=30,
                          onset.rate=0.001,
                          followup.years=30){
    tset <- generate_absence(arm.size,
                             followup.years,
                             onset.rate,
                             with(dset, sojourn.time.min),
                             with(dset, sojourn.time.max))
    cset <- transform(tset, count_clinical=count_onset)
    control_arm <- subset(cset,
                          0 <= clinical_year &
                          clinical_year <= followup.years,
                          select=c(clinical_year, count_clinical))
    control_arm <- ddply(control_arm,
                         .(clinical_year),
                         summarize,
                         arm='control',
                         count_clinical=sum(count_clinical),
                         count_screen=0,
                         count_overdiag=0)
    control_arm <- rename(control_arm, c('clinical_year'='year'))
    sset <- generate_presence(tset,
                              followup.years,
                              screen.start.year,
                              screen.stop.year,
                              with(dset, sensitivity),
                              with(dset, sojourn.time.min),
                              with(dset, sojourn.time.max))
    sset <- generate_overdiag(sset, with(dset, overdiag.rate))
    screen_screen <- ddply(sset,
                           .(screen_year),
                           summarize,
                           arm='screen',
                           count_screen=sum(count_screen),
                           count_overdiag=sum(count_overdiag))
    screen_clinical <- ddply(sset,
                             .(clinical_year, sojourn),
                             summarize,
                             count_clinical=unique(count_clinical))
    screen_clinical <- ddply(screen_clinical,
                             .(clinical_year),
                             summarize,
                             arm='screen',
                             count_clinical=sum(count_clinical))
    screen_arm <- merge(rename(screen_clinical, c('clinical_year'='year')),
                        rename(screen_screen, c('screen_year'='year')),
                        all=TRUE)
    screen_arm[is.na(screen_arm)] <- 0
    screen_arm <- subset(screen_arm, 0 <= year & year <= followup.years)
    merged <- rbind(control_arm, screen_arm)
    return(merged)
}

trial_incidence <- function(screen.stop.year=30){
    tset <- data.frame(sojourn.time.min=c(1, 2, 0),
                       sojourn.time.max=c(1, 2, 4))
    tset <- rbind(transform(tset, overdiag.rate=0),
                  transform(tset, overdiag.rate=0.25))
    tset <- rbind(transform(tset, sensitivity=1),
                  transform(tset, sensitivity=0.5))
    tset <- transform(tset, screen.stop.year=screen.stop.year)
    #tset <- data.frame(sojourn.time.min=0,
    #                   sojourn.time.max=4,
    #                   overdiag.rate=0,
    #                   sensitivity=1)
    #tset <- data.frame(sojourn.time.min=5,
    #                   sojourn.time.max=5,
    #                   overdiag.rate=0,
    #                   sensitivity=0.5)
    tset <- ddply(tset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    overdiag.rate,
                    sensitivity),
                  trial_setting,
                  screen.stop.year=screen.stop.year)
    return(tset)
}
#tset <- trial_incidence(screen.stop.year=30)
#tset <- trial_incidence(screen.stop.year=10)

##################################################
# Illustrate trial setting
##################################################

trial_plot <- function(dset, Sensitivity=1, earlyend=FALSE, saveit=FALSE){
    dset <- droplevels(subset(dset, sensitivity == Sensitivity, select=-sensitivity))
    dset <- transform(dset,
                      sojourn.time=paste0('[',
                                          sojourn.time.min,
                                          ', ',
                                          sojourn.time.max,
                                          ']'),
                      overdiag=paste0('Overdiagnosis rate: ',
                                      round(100*overdiag.rate), '%'))
    dset <- transform(dset, sojourn.time=factor(sojourn.time,
                                            levels=c('[1, 1]',
                                                     '[2, 2]',
                                                     '[0, 4]')))
    dset <- sort_df(dset, vars=c('arm', 'year'))
    dset <- ddply(dset,
                  .(sojourn.time, overdiag, arm),
                  transform,
                  cuminc_nonoverdiag=cumsum(count_clinical+count_screen),
                  cuminc_total=cumsum(count_clinical+count_screen+count_overdiag))
    gg_theme()
    gg <- ggplot(dset)
    gg <- gg+geom_ribbon(data=subset(dset, arm == 'screen'),
                         aes(x=year,
                             ymin=cuminc_nonoverdiag,
                             ymax=cuminc_total),
                         fill='gray90')
    gg <- gg+geom_line(aes(x=year, y=cuminc_total, group=arm, colour=arm), size=0.75)
    gg <- gg+geom_line(aes(x=year, y=cuminc_nonoverdiag, group=arm, colour=arm), size=0.75)
    gg <- gg+geom_hline(yintercept=0, size=0.25)
    gg <- gg+geom_vline(xintercept=0, size=0.25)
    gg <- gg+facet_grid(sojourn.time~overdiag)
    gg <- gg+scale_x_continuous(name='\nYear',
                                breaks=seq(0, 20, by=5),
                                limits=c(0, 20),
                                expand=c(0, 0))
    gg <- gg+scale_y_continuous(name='Cumulative incidence\n',
                                breaks=seq(0, 4000, by=1000),
                                limits=c(0, 4000),
                                expand=c(0, 0))
    gg <- gg+scale_colour_manual(values=c('control'='black',
                                          'screen'='darkgray'))
    #print(gg)
    if(saveit){
        suffix <- ifelse(earlyend, 'earlyend', 'continues')
        filename <- paste('trial', 'prospective', suffix, round(100*Sensitivity), sep='_')
        filename <- paste(filename, 'pdf', sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=10,
               height=7)
    }
}
#trial_plot(tset, Sensitivity=1, earlyend=FALSE, saveit=TRUE)
#trial_plot(tset, Sensitivity=0.5, earlyend=FALSE, saveit=TRUE)
#trial_plot(tset, Sensitivity=1, earlyend=TRUE, saveit=TRUE)
#trial_plot(tset, Sensitivity=0.5, earlyend=TRUE, saveit=TRUE)

##################################################
# Examine simple trial setting
##################################################

trial_setting_simple <- function(dset,
                                 arm.size=50000,
                                 onset.rate=0.001,
                                 followup.years=30){
    # generate trial population of arm.size individuals and
    # record year of clinical diagnosis for batches of relevant
    # cancers that develop in each year with a given sojourn time
    # to serve as a basis for either arm
    tset <- with(dset, generate_absence(arm.size,
                                        followup.years,
                                        onset.rate,
                                        sojourn.time.min,
                                        sojourn.time.max))
    # construct the control arm by "screening" the trial population
    # under 0 sensitivity and counting screen diagnoses in each year
    # of screening for batches of relevant cancers that develop in each
    # year with a given sojourn time
    cset <- with(dset, generate_presence(tset,
                                         followup.years,
                                         screen.start.year,
                                         screen.stop.year,
                                         0,
                                         sojourn.time.min,
                                         sojourn.time.max))
    # append overdiagnoses as a constant fraction of screen
    # diagnoses in each year of screening
    cset <- generate_overdiag(cset, with(dset, overdiag.rate))
    # count control arm screen diagnoses in each year of screening
    control_screen <- ddply(cset,
                           .(screen_year),
                           summarize,
                           arm='control',
                           count_screen=sum(count_screen),
                           count_overdiag=sum(count_overdiag))
    # count clinical diagnoses in each year after first collapsing
    # across screening years
    control_clinical <- ddply(cset,
                             .(clinical_year, sojourn),
                             summarize,
                             count_clinical=unique(count_clinical))
    control_clinical <- ddply(control_clinical,
                             .(clinical_year),
                             summarize,
                             arm='control',
                             count_clinical=sum(count_clinical))
    # merge control arm screen and clinical diagnoses
    control_arm <- merge(rename(control_clinical, c('clinical_year'='year')),
                         rename(control_screen, c('screen_year'='year')),
                         all=TRUE)
    # coerce missing control arm screen diagnoses and overdiagnoses to 0
    control_arm[is.na(control_arm)] <- 0
    # process the trial population under assumed sensitivity
    # to construct the screen arm
    # construct the screen arm by screening the trial population under
    # given sensitivity and counting screen diagnoses in each year of
    # screening for batches of relevant cancers that develop in each
    # year with a given sojourn time
    sset <- with(dset, generate_presence(tset,
                                         followup.years,
                                         screen.start.year,
                                         screen.stop.year,
                                         sensitivity,
                                         sojourn.time.min,
                                         sojourn.time.max,
                                         attendance))
    # append overdiagnoses as a constant fraction of screen
    # diagnoses in each year of screening
    sset <- generate_overdiag(sset, with(dset, overdiag.rate))
    # count screen arm screen diagnoses in each year of screening
    screen_screen <- ddply(sset,
                           .(screen_year),
                           summarize,
                           arm='screen',
                           count_screen=sum(count_screen),
                           count_overdiag=sum(count_overdiag))
    # count clinical diagnoses in each year after first collapsing
    # across screening years
    screen_clinical <- ddply(sset,
                             .(clinical_year, sojourn),
                             summarize,
                             count_clinical=unique(count_clinical))
    screen_clinical <- ddply(screen_clinical,
                             .(clinical_year),
                             summarize,
                             arm='screen',
                             count_clinical=sum(count_clinical))
    # merge screen arm screen and clinical diagnoses
    screen_arm <- merge(rename(screen_clinical, c('clinical_year'='year')),
                        rename(screen_screen, c('screen_year'='year')),
                        all=TRUE)
    # coerce missing screen arm screen diagnoses and overdiagnoses to 0
    screen_arm[is.na(screen_arm)] <- 0
    # merge screen and control arms
    merged <- rbind(control_arm, screen_arm)
    # restrict to given years of follow-up
    merged <- subset(merged, 0 <= year & year <= followup.years)
    return(merged)
}

trial_incidence_simple <- function(screen.stop.year=4){
    tset <- data.frame(sojourn.time.min=c(0, 0),
                       sojourn.time.max=c(12, 6))
    tset <- rbind(transform(tset, overdiag.rate=0),
                  transform(tset, overdiag.rate=0.25))
    tset <- rbind(transform(tset, attendance=0.6),
                  transform(tset, attendance=1))
    tset <- rbind(transform(tset, sensitivity=0.5),
                  transform(tset, sensitivity=1))
    tset <- transform(tset, screen.stop.year=screen.stop.year)
    tset <- ddply(tset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    sensitivity,
                    attendance,
                    overdiag.rate),
                  trial_setting_simple,
                  screen.stop.year=screen.stop.year)
    tset <- transform(tset, count_clinical=ifelse(year == 0, 0, count_clinical))
    return(tset)
}
#tset_stop <- trial_incidence_simple(screen.stop.year=11)

trial_incidence_3types <- function(){
    tset <- data.frame(sojourn.time.min=0,
                       sojourn.time.max=6,
                       sensitivity=0.5,
                       attendance=0.8,
                       overdiag.rate=0.25)
    tset <- rbind(transform(tset, screen.start.year=1,
                                  screen.stop.year=5,
                                  type='screen-stop'),
                  transform(tset, screen.start.year=1,
                                  screen.stop.year=30,
                                  type='screen-continue'),
                  transform(tset, screen.start.year=5,
                                  screen.stop.year=30,
                                  type='screen-continue-delay'))
    tset <- ddply(tset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    sensitivity,
                    attendance,
                    overdiag.rate,
                    screen.start.year,
                    screen.stop.year,
                    type),
                  trial_setting_simple)
    tset <- transform(tset, count_clinical=ifelse(year == 0, 0, count_clinical))
    newscreen <- subset(tset, arm == 'screen' & type == 'screen-continue')
    newcontrol <- subset(tset, arm == 'screen' & type == 'screen-continue-delay')
    newscreen <- transform(newscreen, type='screen-continue-control-start')
    newcontrol <- transform(newcontrol, arm='control', type='screen-continue-control-start')
    tset <- subset(tset, type != 'screen-continue-delay')
    tset <- rbind(tset, newscreen, newcontrol)
    return(tset)
}
#tset_3types <- trial_incidence_3types()

##################################################
# Illustrate ERSPC setting
##################################################

trial_incidence_erspc <- function(){
    tset <- data.frame(sojourn.time.min=0,
                       sensitivity=0.48, # Hakama et al. (2007)
                       attendance=0.76)  # Otto et al. (2010)
    tset <- rbind(transform(tset, screen.start.year=1,
                                  screen.stop.year=5))
    tset <- rbind(transform(tset,
                            sojourn.time.max=3,
                            overdiag.rate=0.4,
                            panel='short'),
                  transform(tset,
                            sojourn.time.max=4,
                            overdiag.rate=0.25,
                            panel='medium'),
                  transform(tset,
                            sojourn.time.max=5,
                            overdiag.rate=0.1,
                            panel='long'))
    screen_arm <- ddply(tset,
                        .(sojourn.time.min,
                          sojourn.time.max,
                          sensitivity,
                          attendance,
                          overdiag.rate,
                          screen.start.year,
                          screen.stop.year,
                          panel),
                        trial_setting_simple,
                        arm.size=72891,
                        onset.rate=0.022,
                        followup.years=5)
    screen_arm <- subset(screen_arm, arm == 'screen')
    screen_arm <- transform(screen_arm,
                            count_clinical=1000*count_clinical/72891,
                            count_screen=1000*count_screen/72891,
                            count_overdiag=1000*count_overdiag/72891)
    control_arm <- ddply(tset,
                         .(sojourn.time.min,
                           sojourn.time.max,
                           sensitivity,
                           attendance,
                           overdiag.rate,
                           screen.start.year,
                           screen.stop.year,
                           panel),
                         trial_setting_simple,
                         arm.size=89352,
                         onset.rate=0.022,
                         followup.years=5)
    control_arm <- subset(control_arm, arm == 'control')
    control_arm <- transform(control_arm,
                             count_clinical=1000*count_clinical/89352,
                             count_screen=1000*count_screen/89352,
                             count_overdiag=1000*count_overdiag/89352)
    tset <- rbind(screen_arm, control_arm)
    tset <- transform(tset, count_clinical=ifelse(year == 0, 0, count_clinical))
    return(tset)
}
#tset_erspc <- trial_incidence_erspc()

##################################################
# Illustrate simple trial setting
##################################################

trial_plot_simple <- function(dset,
                              Sensitivity=1,
                              Overdiagnosis=0,
                              Cumulative=FALSE,
                              minyear=0,
                              maxyear=8,
                              line.size=0.15,
                              point.size=0.2,
                              point.shape=19,
                              cont=FALSE,
                              saveit=FALSE){
    dset <- droplevels(subset(dset,
                              sensitivity == Sensitivity &
                              overdiag.rate == Overdiagnosis,
                              select=-c(sensitivity, overdiag.rate)))
    dset <- transform(dset, sojourn.time=paste0('[',
                                                sojourn.time.min,
                                                ', ',
                                                sojourn.time.max,
                                                ']'))
    dset <- ddply(dset,
                  .(sojourn.time),
                  transform,
                  sojourn.mean=paste('Mean:',
                                     mean(seq(unique(sojourn.time.min),
                                              unique(sojourn.time.max))),
                                     'years'),
                  sojourn.range=paste('Range:\n',
                                      unique(sojourn.time.max-sojourn.time.min),
                                      'years'))
    dset <- transform(dset,
                      line.size=line.size,
                      point.size=point.size,
                      point.shape=point.shape,
                      sojourn.time=factor(sojourn.time,
                                          levels=c('[2, 2]',
                                                   '[0, 4]',
                                                   '[4, 4]',
                                                   '[2, 6]')),
                      sojourn.mean=factor(sojourn.mean),
                      sojourn.range=relevel(factor(sojourn.range),
                                            ref=c('Range:\n 0 years')),
                      arm=relevel(factor(arm), ref='screen'))
    dset <- ddply(dset,
                  .(sojourn.time, arm),
                  transform,
                  total=count_clinical+count_screen+count_overdiag)
    ylabel <- 'Total number of cancers\n'
    if(Cumulative){
        dset <- ddply(dset,
                      .(sojourn.time, arm),
                      transform,
                      nonoverdiag=cumsum(count_clinical+count_screen),
                      total=cumsum(total))
        ylabel <- 'Cumulative number of cancers\n'
    }
    gg_theme(panel.border=element_rect(fill='transparent'))
    gg <- ggplot(dset)
    gg <- gg+geom_vline(xintercept=seq(3), linetype='dashed', colour='gray')
    #gg <- gg+geom_vline(xintercept=-0.5, size=0.5)
    #gg <- gg+geom_hline(yintercept=-50, size=0.5)
    if(Overdiagnosis > 0){
        dset.screen <- subset(dset, arm == 'screen')
        gg <- gg+geom_ribbon(data=dset.screen,
                             aes(x=year,
                                 ymin=nonoverdiag,
                                 ymax=total,
                                 size=line.size),
                             colour='lightgreen',
                             fill=alpha('lightgreen', 0.4))
        gg <- gg+geom_line(data=dset.screen,
                           aes(x=year,
                               y=nonoverdiag,
                               size=line.size),
                           colour='lightgreen')
        gg <- gg+geom_point(data=dset.screen,
                            aes(x=year,
                                y=nonoverdiag,
                                size=point.size,
                                shape=as.character(point.shape)),
                            colour='lightgreen')
    }
    gg <- gg+geom_line(aes(x=year,
                           y=total,
                           group=arm,
                           colour=arm,
                           size=line.size))
    gg <- gg+geom_point(aes(x=year,
                            y=total,
                            group=arm,
                            colour=arm,
                            size=point.size,
                            shape=as.character(point.shape)))
    gg <- gg+facet_grid(sojourn.range~sojourn.mean)
    gg <- gg+scale_x_continuous(name='\nYears of follow-up',
                                breaks=seq(minyear, maxyear),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0.5))
    gg <- gg+scale_y_continuous(name=ylabel,
                                breaks=seq(0, 500, by=100),
                                limits=c(0, 500),
                                expand=c(0, 50))
    gg <- gg+scale_colour_manual(values=c('control'=muted('aliceblue'),
                                          'screen'='hotpink'))
    print(gg)
    if(saveit){
        suffix <- ifelse(cont, 'continue', 'stop')
        filename <- paste('trial',
                          'prospective',
                          suffix,
                          round(100*Sensitivity),
                          round(100*Overdiagnosis),
                          sep='_')
        filename <- paste(filename, 'svg', sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=10,
               height=7)
    }
}
#trial_plot_simple(tset_stop,
#                  Sensitivity=1,
#                  Overdiagnosis=0,
#                  Cumulative=TRUE,
#                  cont=FALSE,
#                  saveit=FALSE)
#trial_plot_simple(tset_stop,
#                  Sensitivity=0.5,
#                  Overdiagnosis=0,
#                  Cumulative=TRUE,
#                  cont=FALSE,
#                  saveit=FALSE)
#trial_plot_simple(tset_stop,
#                  Sensitivity=0.5,
#                  Overdiagnosis=0.25,
#                  Cumulative=TRUE,
#                  cont=FALSE,
#                  saveit=FALSE)

##################################################
# Illustrate trial setting
##################################################

trial_plot_revised <- function(dset,
                               Overdiagnosis=0,
                               Cumulative=FALSE,
                               minyear=0,
                               maxyear=21,
                               minyear.text=1,
                               text.size=3,
                               text.offset=75,
                               stable.arrows=FALSE,
                               make.table=FALSE,
                               saveit=FALSE){
    dset <- droplevels(subset(dset,
                              minyear <= year & year <= maxyear &
                              overdiag.rate == Overdiagnosis,
                              select=-overdiag.rate))
    dset <- transform(dset,
                      sensitivity=factor(paste0('Effective sensitivity: ',
                                                round(100*sensitivity),
                                                '%')),
                      attendance=factor(paste0('Screening attendance: ',
                                                round(100*attendance),
                                                '%')),
                      sojourn.time=factor(paste('Max preclinical detectable period:',
                                                 sojourn.time.max)),
                      arm=relevel(factor(arm), ref='screen'))
    dset <- ddply(dset,
                  .(sensitivity, attendance, sojourn.time, arm),
                  transform,
                  total=count_clinical+count_screen+count_overdiag)
    if(Cumulative){
        dset <- ddply(dset,
                      .(sensitivity, attendance, sojourn.time, arm),
                      transform,
                      count_clinical=cumsum(count_clinical),
                      count_screen=cumsum(count_screen),
                      nonoverdiag=cumsum(count_clinical+count_screen),
                      total=cumsum(total))
        ylabel <- 'Cumulative number of cancers\n'
        ymax <- 1200
        ystep <- 200
        yoffset <- 50
    } else {
        dset <- ddply(dset,
                      .(sensitivity, attendance, sojourn.time, arm),
                      transform,
                      nonoverdiag=count_clinical+count_screen)
        ylabel <- 'Number of cancers\n'
        ymax <- 200
        ystep <- 50
        yoffset <- 10
    }
    dset <- transform(dset, label=apply(data.frame(sensitivity,
                                                   attendance,
                                                   sojourn.time),
                                        1,
                                        paste,
                                        collapse='\n'))
    facets <- with(dset, as.data.frame(table(sensitivity,
                                             attendance,
                                             sojourn.time)))
    # rows 7, 8, 5, 3 specify desired base case and variants
    #facets <- droplevels(subset(facets[c(7, 8, 5, 3), ], select=-Freq))
    facets <- droplevels(subset(facets[c(1, 2, 5, 6), ], select=-Freq))
    facets <- summarize(facets, label=apply(facets, 1, paste, collapse='\n'))
    dset <- merge(facets, dset, by='label')
    dset <- transform(dset, label=factor(label, levels=facets$label))
    dset <- sort_df(dset, vars=c('label', 'year'))
    if(Cumulative){
        dset <- ddply(dset,
                      .(label),
                      function(x){
                          tot_count <- with(x, total[arm == 'screen'])
                          ndx_count <- with(x, nonoverdiag[arm == 'screen'])
                          cdx_count <- with(x, count_clinical[arm == 'control'])
                          screen_year <- with(x, year[arm == 'screen'])
                          control_year <- with(x, year[arm == 'control'])
                          year_offset <- cdx_count > 0
                          rough <- year_offset & (tot_count-ndx_count)/(tot_count-cdx_count) > 0.9
                          rough_year <- min(which(rough))
                          exact <- year_offset & abs(ndx_count-cdx_count) < 0.1
                          exact_year <- min(which(exact))
                          print(with(x, as.character(unique(label))))
                          cat('rough:', rough_year, '\n')
                          cat('exact:', exact_year, '\n')
                          stopifnot(isTRUE(all.equal(screen_year[exact_year],
                                                     control_year[exact_year])))
                          stopifnot(isTRUE(all.equal(screen_year[rough_year],
                                                     control_year[rough_year])))
                          x <- transform(x,
                                         rough=arm == 'screen' & year == screen_year[rough_year],
                                         rough_tag='90%',
                                         exact=arm == 'screen' & year == screen_year[exact_year],
                                         exact_tag='100%')
                          return(x)
                      })
        if(make.table){
            ests <- ddply(dset,
                          .(label),
                          function(x){
                              end_year <- x[with(x, exact), 'year']
                              xset <- subset(x, year <= end_year)
                              xset <- summarize(xset,
            followup=end_year,
            control=total[arm == 'control' & year == end_year],
            screen=total[arm == 'screen' & year == end_year],
            actual=sum(count_overdiag),
            scases=sum(count_screen[arm == 'screen' & year == end_year]))
                              return(xset)
                          })
            ests <- transform(ests,
                              excess=screen-control,
                              ratio=screen/control)
            write.csv(ests,
                      file=file.path(datapath, 'trial_estimates.csv'),
                      quote=TRUE,
                      row.names=FALSE)
        }
        dset.text <- subset(dset, year %in% seq(minyear.text, maxyear))
        dset.text <- transform(dset.text,
                               text_size=text.size,
                               text_offset=text.offset)
    }
    gg_theme(panel.margin=unit(0.03, 'npc'),
             panel.border=element_rect(fill='transparent'),
             strip.text.x=element_text(size=16, hjust=0))
    gg <- ggplot(dset)
    gg <- gg+geom_vline(data=subset(dset, arm == 'screen'),
                        xintercept=seq(10),
                        size=1,
                        linetype='dashed',
                        colour='lightgray')
    if(Overdiagnosis > 0){
        dset.screen <- subset(dset, arm == 'screen')
        gg <- gg+geom_ribbon(data=dset.screen,
                             aes(x=year,
                                 ymin=nonoverdiag,
                                 ymax=total),
                             size=0.5,
                             colour='lightgreen',
                             fill=alpha('lightgreen', 0.4))
        gg <- gg+geom_line(data=dset.screen,
                           aes(x=year, y=nonoverdiag),
                           size=1,
                           colour='lightgreen')
        gg <- gg+geom_point(data=dset.screen,
                            aes(x=year, y=nonoverdiag),
                            size=3,
                            shape=19,
                            colour='lightgreen')
    }
    gg <- gg+geom_line(aes(x=year, y=total, group=arm, colour=arm),
                       size=1)
    gg <- gg+geom_point(aes(x=year,
                            y=total,
                            group=arm,
                            colour=arm),
                        size=3,
                        shape=19)
    gg <- gg+facet_wrap(~label, nrow=2)
    gg <- gg+scale_x_continuous(name='\nYears of follow-up',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0.5))
    gg <- gg+scale_y_continuous(name=ylabel,
                                #breaks=seq(0, ymax, by=ystep),
                                #limits=c(0, ymax),
                                expand=c(0, yoffset))
    gg <- gg+scale_size(range=c(2, 2))
    gg <- gg+scale_colour_manual(values=c('control'=muted('skyblue'),
                                          'screen'='lightgreen'))
    if(Cumulative){
        if(Overdiagnosis == 0){
            gg <- gg+geom_text(data=dset.text,
                               aes(x=year,
                                   y=nonoverdiag+text_offset,
                                   label=paste0('+', round(count_screen)),
                                   size=text_size))
        } else {
            gg <- gg+geom_text(data=subset(dset.text, arm == 'screen'),
                               aes(x=year,
                                   y=total+text_offset,
                                   label=round(total),
                                   size=text_size))
            if(stable.arrows){
                gg <- gg+geom_segment(data=subset(dset, exact),
                                      aes(x=year,
                                          xend=year,
                                          y=total+300,
                                          yend=total+150),
                                      colour='gray30',
                                      arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
                gg <- gg+geom_text(data=subset(dset, exact),
                                   aes(x=year,
                                       y=total+350,
                                       label=exact_tag),
                                   colour='gray30',
                                   vjust=0,
                                   size=4)
                #gg <- gg+geom_segment(data=subset(dset, rough),
                #                      aes(x=year,
                #                          xend=year,
                #                          y=total+300,
                #                          yend=total+150),
                #                      colour='gray60',
                #                      arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
                #gg <- gg+geom_text(data=subset(dset, rough),
                #                   aes(x=year,
                #                       y=total+350,
                #                       label=rough_tag),
                #                   colour='gray60',
                #                   vjust=0,
                #                   size=4)
            }
        }
        gg <- gg+geom_text(data=subset(dset.text, arm == 'control'),
                           aes(x=year,
                               y=count_clinical-text_offset,
                               label=round(count_clinical),
                               size=text_size))
    }
    print(gg)
    if(saveit){
        filename <- paste('trial_revised', round(100*Overdiagnosis), sep='_')
        if(Cumulative)
            filename <- paste(filename, 'cumulative', sep='_')
        if(stable.arrows)
            filename <- paste(filename, 'arrows', sep='_')
        filename <- paste(filename, '2015-02-11', sep='_')
        filename <- paste(filename, 'png', sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=12,
               height=7)
    }
}

#trial_plot_revised(tset_stop,
#                   Overdiagnosis=0.25,
#                   Cumulative=TRUE,
#                   stable.arrows=TRUE,
#                   make.table=FALSE,
#                   saveit=FALSE)
#trial_plot_revised(tset_stop,
#                   Overdiagnosis=0.25,
#                   Cumulative=FALSE,
#                   stable.arrows=TRUE,
#                   make.table=FALSE,
#                   saveit=FALSE)

##################################################
# Illustrate 3 types of trial settings:
# 1. both arms stop screening
# 2. screen arm continues screening
# 3. both arms continue screening
##################################################

trial_plot_3types <- function(dset,
                              minyear=0,
                              maxyear=12,
                              text.offset=75,
                              make.table=FALSE,
                              stable.arrows=FALSE,
                              ext='png',
                              saveit=FALSE){
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    dset <- ddply(dset,
                  .(arm, type),
                  transform,
                  count_nonoverdiag=count_clinical+count_screen,
                  count_total=count_clinical+count_screen+count_overdiag)
    dset <- transform(dset, measure='annual')
    cset <- transform(dset, measure='cumulative')
    cset <- ddply(cset,
                  .(arm, type),
                  transform,
                  count_clinical=cumsum(count_clinical),
                  count_screen=cumsum(count_screen),
                  count_overdiag=cumsum(count_overdiag),
                  count_nonoverdiag=cumsum(count_nonoverdiag),
                  count_total=cumsum(count_total))
    dset <- rbind(dset, cset)
    dset <- ddply(dset,
                  .(type, measure),
                  function(x){
                      with(x, cat(paste(rep('-', 40), collapse=''),
                                  '\nsetting:', as.character(unique(type)),
                                  '\nmeasure:', as.character(unique(measure)),
                                  '\n'))
                      # extract years for each arm
                      screen_year <- with(x, year[arm == 'screen'])
                      control_year <- with(x, year[arm == 'control'])
                      # extract screen diagnoses for each arm
                      screen_screen <- with(x, count_screen[arm == 'screen'])
                      control_screen <- with(x, count_screen[arm == 'control'])
                      # extract overdiagnoses for each arm
                      screen_overdiag <- with(x, count_overdiag[arm == 'screen'])
                      control_overdiag <- with(x, count_overdiag[arm == 'control'])
                      # extract total diagnoses for each arm
                      screen_total <- with(x, count_total[arm == 'screen'])
                      control_total <- with(x, count_total[arm == 'control'])
                      # calculate empirical excess and true overdiagnosis
                      excess_overdiag <- screen_overdiag-control_overdiag
                      excess_total <- screen_total-control_total
                      # specify valid years for evaluating empirical excess
                      if(any(control_screen > 0))
                          year_offset <- FALSE#control_screen > 0
                      else
                          year_offset <- screen_screen > 0
                      # calculate first year excess matches true overdiagnosis
                      exact <- year_offset & abs(excess_overdiag/excess_total-1) < 0.000001
                      rough <- year_offset & abs(excess_overdiag/excess_total-1) < 0.1
                      exact_year <- min(which(exact))
                      rough_year <- min(which(rough))
                      cat('exact:', screen_year[exact_year], '\n')
                      cat('rough:', screen_year[rough_year], '\n')
                      stopifnot(isTRUE(all.equal(screen_year[exact_year],
                                                 control_year[exact_year])))
                      stopifnot(isTRUE(all.equal(screen_year[rough_year],
                                                 control_year[rough_year])))
                      #if(with(x, unique(type) == 'screen-stop'))
                      #    browser()
                      #if(with(x, unique(type) == 'screen-continue'))
                      #    browser()
                      #if(with(x, unique(type) == 'screen-continue-control-start'))
                      #    browser()
                      # append exact and rough approximation years to stratum
                      x <- transform(x,
                                     rough=arm == 'screen' & year == screen_year[rough_year],
                                     rough_tag='90%',
                                     exact=arm == 'screen' & year == screen_year[exact_year],
                                     exact_tag='Unbiased')
                      return(x)
                  })
    dset <- transform(dset,
                      measure=sapply(as.character(measure), switch,
                                     'annual'='\tAnnual\n\tincidence',
                                     'cumulative'='\tCumulative\n\tincidence'),
                      type=sapply(as.character(type), switch,
                                'screen-stop'=paste('Screen arm receives tests during first 4 years',
                                                    'Control arm receives no tests',
                                                    sep='\n'),
                                'screen-continue'=paste('Screen arm receives tests in all years',
                                                        'Control arm receives no tests',
                                                        sep='\n'),
                                'screen-continue-control-start'=paste('Screen arm receives tests in all years',
                                                              'Control arm receives tests after first 4 years',
                                                              sep='\n')))
    dset <- transform(dset, type=factor(type, levels=unique(type)))
    if(make.table){
        ests <- ddply(dset,
                      .(type),
                      function(x){
                          end_year <- x[with(x, exact), 'year']
                          xset <- subset(x, year <= end_year)
                          xset <- summarize(xset,
                                followup=end_year,
                                control=count_total[arm == 'control' & year == end_year],
                                screen=count_total[arm == 'screen' & year == end_year],
                                actual=sum(count_overdiag),
                                scases=sum(count_screen[arm == 'screen' & year == end_year]))
                          return(xset)
                      })
        ests <- transform(ests,
                          excess=screen-control,
                          ratio=screen/control)
        write.csv(ests,
                  file=file.path(datapath, 'trial_estimates.csv'),
                  quote=TRUE,
                  row.names=FALSE)
    }
    gg_theme(panel.margin=unit(0.02, 'npc'),
             panel.border=element_rect(fill='transparent'),
             strip.text.x=element_text(size=12, hjust=0),
             strip.text.y=element_text(size=14, hjust=0),
             axis.title=element_text(size=18))
    gg <- ggplot(dset)
    # highlight overdiagnosis in the screen arm
    {
        dset.screen <- subset(dset, arm == 'screen')
        gg <- gg+geom_ribbon(data=dset.screen,
                             aes(x=year,
                                 ymin=count_nonoverdiag,
                                 ymax=count_total),
                             size=0.5,
                             fill=alpha('darkgreen', 0.4),
                             colour=alpha('darkgreen', 0.4))
    }
    # highlight overdiagnosis in the control arm
    {
        dset.control <- subset(dset, arm == 'control' & type == 'Screen arm receives tests in all years\nControl arm receives tests after first 4 years')
        gg <- gg+geom_ribbon(data=dset.control,
                             aes(x=year,
                                 ymin=count_nonoverdiag,
                                 ymax=count_total),
                             size=0.5,
                             fill=alpha('purple', 0.4),
                             colour=alpha('purple', 0.4))
    }
    gg <- gg+geom_line(aes(x=year, y=count_total, colour=arm),
                       alpha=0.75,
                       size=1)
    gg <- gg+geom_point(aes(x=year,
                            y=count_total,
                            group=arm,
                            colour=arm),
                        alpha=1,
                        size=2,
                        shape=19)
    gg <- gg+facet_grid(measure~type, scales='free_y')
    gg <- gg+scale_x_continuous(name='\nYears of follow-up',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0.5))
    gg <- gg+scale_y_continuous(name='Number of cancers\n',
                                #breaks=seq(0, ymax, by=ystep),
                                #limits=c(0, ymax),
                                expand=c(0.05, 0))
    gg <- gg+scale_size(range=c(2, 2))
    gg <- gg+geom_blank(data=subset(dset, measure == '\tAnnual\n\tincidence'), aes(y=150))
    gg <- gg+scale_colour_manual(values=c('control'=alpha('purple', 0.4),
                                          'screen'=alpha('darkgreen', 0.4)))
    gg <- gg+scale_fill_manual(values=c('control'=alpha('purple', 0.4),
                                        'screen'=alpha('darkgreen', 0.4)))
    {# highlight screen events in each arm
        gg <- gg+geom_segment(data=subset(dset, arm == 'screen' &
                                                count_screen > 0 &
                                                measure=='\tAnnual\n\tincidence'),
                              aes(x=year,
                                  xend=year,
                                  y=-Inf,
                                  yend=0),
                              size=1,
                              colour=alpha('darkgreen', 0.4))
        gg <- gg+geom_segment(data=subset(dset, arm == 'control' &
                                                count_screen > 0 &
                                                measure=='\tAnnual\n\tincidence'),
                              aes(x=year,
                                  xend=year,
                                  y=0,
                                  yend=10),
                              size=1,
                              colour=alpha('purple', 0.4))
    }
    if(stable.arrows){
        gg <- gg+geom_segment(data=subset(dset, exact),
                              aes(x=year,
                                  xend=year,
                                  y=1.75*count_total,
                                  yend=1.25*count_total),
                              colour='gray30',
                              arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
        gg <- gg+geom_text(data=subset(dset, exact),
                           aes(x=year,
                               y=2*count_total,
                               label=exact_tag),
                           colour='gray30',
                           vjust=0,
                           size=4)
        #gg <- gg+geom_segment(data=subset(dset, rough),
        #                      aes(x=year,
        #                          xend=year,
        #                          y=1.75*count_total,
        #                          yend=1.25*count_total),
        #                      colour='gray60',
        #                      arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
        #gg <- gg+geom_text(data=subset(dset, rough),
        #                   aes(x=year,
        #                       y=2*count_total,
        #                       label=rough_tag),
        #                   colour='gray60',
        #                   vjust=0,
        #                   size=4)
    }
    print(gg)
    if(saveit){
        filename <- 'trial_3types'
        datestamp <- '2015-04-13'
        if(stable.arrows)
            filename <- paste(filename, 'arrows', sep='_')
        filename <- paste(filename, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=15,
               height=5)
    }
}
#trial_plot_3types(tset_3types,
#                  make.table=FALSE,
#                  stable.arrows=TRUE,
#                  ext='svg',
#                  saveit=FALSE)

##################################################
# Illustrate 2 types of trial settings as above
# but without overdiagnosis
##################################################

trial_incidence_2types <- function(){
    tset <- data.frame(sojourn.time.min=0,
                       sojourn.time.max=6,
                       sensitivity=0.5,
                       attendance=0.8,
                       overdiag.rate=0)
    tset <- rbind(transform(tset, screen.start.year=1,
                                  screen.stop.year=5,
                                  type='screen-stop'),
                  transform(tset, screen.start.year=1,
                                  screen.stop.year=30,
                                  type='screen-continue'))
    tset <- ddply(tset,
                  .(sojourn.time.min,
                    sojourn.time.max,
                    sensitivity,
                    attendance,
                    overdiag.rate,
                    screen.start.year,
                    screen.stop.year,
                    type),
                  trial_setting_simple)
    tset <- transform(tset, count_clinical=ifelse(year == 0, 0, count_clinical))
    return(tset)
}
#tset_2types <- trial_incidence_2types()

trial_plot_2types <- function(dset,
                              minyear=0,
                              maxyear=12,
                              text.offset=75,
                              stable.arrows=FALSE,
                              ext='png',
                              saveit=FALSE){
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    dset <- ddply(dset,
                  .(arm, type),
                  transform,
                  count_nonoverdiag=count_clinical+count_screen,
                  count_total=count_clinical+count_screen+count_overdiag)
    dset <- transform(dset, measure='annual')
    cset <- transform(dset, measure='cumulative')
    cset <- ddply(cset,
                  .(arm, type),
                  transform,
                  count_clinical=cumsum(count_clinical),
                  count_screen=cumsum(count_screen),
                  count_overdiag=cumsum(count_overdiag),
                  count_nonoverdiag=cumsum(count_nonoverdiag),
                  count_total=cumsum(count_total))
    dset <- rbind(dset, cset)
    dset <- ddply(dset,
                  .(type, measure),
                  function(x){
                      with(x, cat(paste(rep('-', 40), collapse=''),
                                  '\nsetting:', as.character(unique(type)),
                                  '\nmeasure:', as.character(unique(measure)),
                                  '\n'))
                      # extract years for each arm
                      screen_year <- with(x, year[arm == 'screen'])
                      control_year <- with(x, year[arm == 'control'])
                      # extract screen diagnoses for each arm
                      screen_screen <- with(x, count_screen[arm == 'screen'])
                      control_screen <- with(x, count_screen[arm == 'control'])
                      # extract overdiagnoses for each arm
                      screen_overdiag <- with(x, count_overdiag[arm == 'screen'])
                      control_overdiag <- with(x, count_overdiag[arm == 'control'])
                      # extract total diagnoses for each arm
                      screen_total <- with(x, count_total[arm == 'screen'])
                      control_total <- with(x, count_total[arm == 'control'])
                      # calculate empirical excess and true overdiagnosis
                      excess_overdiag <- screen_overdiag-control_overdiag
                      excess_total <- screen_total-control_total
                      # specify valid years for evaluating empirical excess
                      year_offset <- screen_total > 0
                      # calculate first year excess matches true overdiagnosis
                      exact <- year_offset & abs(excess_total) < 0.000001
                      rough <- year_offset & abs(excess_overdiag/excess_total-1) < 0.1
                      exact_year <- min(which(exact))
                      cat('exact:', screen_year[exact_year], '\n')
                      stopifnot(isTRUE(all.equal(screen_year[exact_year],
                                                 control_year[exact_year])))
                      x <- transform(x,
                                     exact=arm == 'screen' & year == screen_year[exact_year],
                                     exact_tag='Unbiased')
                      return(x)
                  })
    dset <- transform(dset, measure=sapply(as.character(measure),
                                           switch,
                                           'annual'='\tAnnual\n\tincidence',
                                           'cumulative'='\tCumulative\n\tincidence'))
    dset <- transform(dset, type=factor(type, levels=unique(type)))
    gg_theme(panel.margin=unit(0.02, 'npc'),
             panel.border=element_rect(fill='transparent'),
             strip.text.x=element_text(size=12, hjust=0),
             strip.text.y=element_text(size=14, hjust=0),
             axis.title=element_text(size=18))
    gg <- ggplot(dset)
    gg <- gg+geom_line(aes(x=year, y=count_total, colour=arm),
                       alpha=0.75,
                       size=1)
    gg <- gg+geom_point(aes(x=year,
                            y=count_total,
                            group=arm,
                            colour=arm),
                        alpha=1,
                        size=2,
                        shape=19)
    gg <- gg+facet_grid(measure~type, scales='free_y')
    gg <- gg+scale_x_continuous(name='\nYears of follow-up',
                                breaks=seq(minyear, maxyear, by=2),
                                limits=c(minyear, maxyear),
                                expand=c(0, 0.5))
    gg <- gg+scale_y_continuous(name='Number of cancers\n',
                                expand=c(0.05, 0))
    gg <- gg+geom_blank(data=subset(dset, measure == '\tAnnual\n\tincidence'), aes(y=120))
    gg <- gg+scale_colour_manual(values=c('control'=alpha('purple', 0.4),
                                          'screen'=alpha('darkgreen', 0.4)))
    {# highlight screen events in each arm
        gg <- gg+geom_segment(data=subset(dset, arm == 'screen' &
                                                count_screen > 0 &
                                                measure=='\tAnnual\n\tincidence'),
                              aes(x=year,
                                  xend=year,
                                  y=-Inf,
                                  yend=0),
                              size=1,
                              colour=alpha('darkgreen', 0.4))
    }
    if(stable.arrows){
        gg <- gg+geom_segment(data=subset(dset, exact),
                              aes(x=year,
                                  xend=year,
                                  y=1.75*count_total,
                                  yend=1.25*count_total),
                              colour='gray30',
                              arrow=arrow(length=unit(0.3, 'cm'), type='closed'))
        gg <- gg+geom_text(data=subset(dset, exact),
                           aes(x=year,
                               y=2*count_total,
                               label=exact_tag),
                           colour='gray30',
                           vjust=0,
                           size=4)
    }
    print(gg)
    if(saveit){
        filename <- 'trial_2types'
        datestamp <- '2015-10-23'
        if(stable.arrows)
            filename <- paste(filename, 'arrows', sep='_')
        filename <- paste(filename, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=10,
               height=5)
    }
}
#trial_plot_2types(tset_2types,
#                  stable.arrows=TRUE,
#                  ext='svg',
#                  saveit=FALSE)

##################################################
# Observed ERSPC incidence (with rescaled years)
##################################################

obs_erspc <- data.frame(year=c(9, 11, 13)/4,
                        screen=1000*c(5990, 6963, 7408)/72891,
                        control=1000*c(4307, 5396, 6107)/89352)
obs_erspc <- melt(obs_erspc, id.vars='year')
obs_erspc <- rename(obs_erspc, c(variable='arm', value='count_total'))

nnd_erspc <- ddply(obs_erspc, .(year), summarize, inc_diff=-diff(count_total))
nnd_erspc <- transform(nnd_erspc, mort_diff=c(0.71, 1.07, 1.28))
nnd_erspc <- transform(nnd_erspc, nnd=inc_diff/mort_diff)

##################################################
# Illustrate ERSPC setting
##################################################

trial_plot_erspc <- function(dset,
                             erspc,
                             minyear=0,
                             maxyear=4,
                             ext='png',
                             saveit=FALSE){
    years <- seq(minyear, maxyear, by=0.25)
    dset <- droplevels(subset(dset, minyear <= year & year <= maxyear))
    dset <- ddply(dset,
                  .(arm, panel),
                  transform,
                  count_nonoverdiag=count_clinical+count_screen,
                  count_total=count_clinical+count_screen+count_overdiag)
    dset <- ddply(dset,
                  .(arm, panel),
                  transform,
                  count_clinical=cumsum(count_clinical),
                  count_screen=cumsum(count_screen),
                  count_overdiag=cumsum(count_overdiag),
                  count_nonoverdiag=cumsum(count_nonoverdiag),
                  count_total=cumsum(count_total))
    dset <- transform(dset,
                      sojourn.time=factor(paste('Max preclinical detectable period:',
                                                4*sojourn.time.max,
                                                'years')),
                      overdiagnosis=factor(paste0('True overdiagnosis rate: ',
                                                  round(100*overdiag.rate),
                                                  '%')))
    dset <- transform(dset, label=apply(data.frame(sojourn.time,
                                                   overdiagnosis),
                                        1,
                                        paste,
                                        collapse='\n'))
    dset <- subset(dset, select=-c(sojourn.time, overdiagnosis))
    calculations <- ddply(dset,
                          .(panel),
                          function(x, years=c(9, 11, 13)){
                              with(x, cat(paste(rep('-', 40), collapse=''),
                                          '\npanel:', as.character(unique(panel)),
                                          '\n'))
                              screen_arm <- subset(x,
                                                   arm == 'screen',
                                                   select=c(year, count_overdiag))
                              screen_arm <- transform(screen_arm, year=year*4)
                              estimates <- data.frame(year=years,
                                                      estimate=with(screen_arm,
                                                                    approx(year,
                                                                           count_overdiag,
                                                                           xout=years))$y,
                                                      panel=with(x, unique(panel)))
                              return(estimates)
                          })
    calculations <- melt(calculations, id.vars=c('year', 'panel'))
    calculations <- subset(calculations, select=-variable)
    casted <- cast(calculations, year~panel)
    #round(casted)
    gg_theme(panel.margin=unit(0.02, 'npc'),
             panel.border=element_rect(fill='transparent'),
             strip.text.x=element_text(size=12, hjust=0),
             strip.text.y=element_text(size=14, hjust=0),
             axis.title=element_text(size=18))
    gg <- ggplot(dset)
    # illustrate screening tests
    gg <- gg+geom_vline(xintercept=seq(0, 3), size=0.15, linetype='dashed')
    # highlight overdiagnosis in the screen arm
    {
        dset.screen <- subset(dset, arm == 'screen')
        gg <- gg+geom_ribbon(data=dset.screen,
                             aes(x=year,
                                 ymin=count_nonoverdiag,
                                 ymax=count_total),
                             size=0.5,
                             fill=alpha('darkgreen', 0.4),
                             colour=alpha('darkgreen', 0.4))
    }
    gg <- gg+geom_line(aes(x=year, y=count_total, colour=arm),
                       alpha=0.75,
                       size=1)
    gg <- gg+geom_point(aes(x=year,
                            y=count_total,
                            group=arm,
                            colour=arm),
                        alpha=1,
                        size=2,
                        shape=19)
    gg <- gg+geom_point(data=erspc,
                        aes(x=year,
                            y=count_total,
                            group=arm,
                            colour=arm),
                        colour='black',
                        size=3,
                        shape=19)
    gg <- gg+facet_grid(.~label)
    gg <- gg+scale_x_continuous(name='\nYears of follow-up',
                                breaks=years,
                                labels=4*years,
                                limits=c(minyear, maxyear),
                                expand=c(0, 0.15))
    gg <- gg+scale_y_continuous(name='Cumulative prostate cancers\nper 1000 men\n',
                                breaks=seq(0, 150, by=25),
                                expand=c(0.025, 0))
    gg <- gg+scale_colour_manual(values=c('control'=alpha('purple', 0.4),
                                          'screen'=alpha('darkgreen', 0.4)))
    gg <- gg+scale_fill_manual(values=c('control'=alpha('purple', 0.4),
                                        'screen'=alpha('darkgreen', 0.4)))
    print(gg)
    if(saveit){
        filename <- 'trial_erspc'
        datestamp <- '2015-08-26'
        filename <- paste(filename, datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(file.path(plotpath, filename),
               plot=gg,
               width=15,
               height=5)
    }
}
#trial_plot_erspc(tset_erspc,
#                 obs_erspc,
#                 ext='svg',
#                 saveit=FALSE)

