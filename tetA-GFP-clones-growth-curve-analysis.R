## tetA-GFP-clones-growth-curve-analysis.R by Rohan Maddamsetti.

library(tidyverse)
library(cowplot)
library(broom)

## TODO: POLISH PLOTS, AND SUMMARIZE DATA/RESULTS!

#############################################################
## Global constants for calculating max growth rate and time lags based on an OD600 interval of interest.

SUBTRACTED.OD600.INTERVAL.MIN <- 0.02
SUBTRACTED.OD600.INTERVAL.MAX <- 0.14
SUBTRACTED.OD600.INTERVAL.MIDPOINT <- (SUBTRACTED.OD600.INTERVAL.MIN + SUBTRACTED.OD600.INTERVAL.MAX)/2


RAW.OD600.INTERVAL.MIN <- 0.3
RAW.OD600.INTERVAL.MAX <- 0.4
RAW.OD600.INTERVAL.MIDPOINT <- (RAW.OD600.INTERVAL.MIN + RAW.OD600.INTERVAL.MAX)/2

#############################################################
## Functions

well.to.Tet <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the tetracycline concentration.
    rowletter <- substring(well,1,1)
    LB.rows <- c("A", "B", "C")
    LB.Tet50.rows <- c("F", "G", "H")
    Blank.rows <- c("D", "E")
    if (rowletter %in% LB.Tet50.rows) return(50)
    if (rowletter %in% LB.rows) return(0)
    if (rowletter %in% Blank.rows) return(0)
    return(NA)
}


well.to.Plasmid <- function(well) {
    ## This helper function maps the row of the well in
    ## the 96-well plate to the plasmid treatment.
    rowletter <- substring(well,1,1)
    no.plasmid.rows <- c("A", "F") 
    p15A.plasmid.rows <- c("B", "G")
    pUC.plasmid.rows <- c("C", "H")
    Blank.rows <- c("D", "E")
    if (rowletter %in% no.plasmid.rows) return("No plasmid")
    if (rowletter %in% p15A.plasmid.rows) return("p15A")
    if (rowletter %in% pUC.plasmid.rows) return("pUC")
    if (rowletter %in% Blank.rows) return("LB Blank")
    return(NA)
}


well.to.population <- function(well) {
    column.num <- as.numeric(substring(well,2))
    ## subtract one since the first and last columns of the plate were blanks.
    population <- column.num - 1
    return(population)
}


tidy.my.tetA.GFP.clone.data <- function(long.format.data, plate.reader.measurement.type) {
    ## only OD600, GFP55, and GFP100 data allowed here.
    stopifnot(plate.reader.measurement.type %in% c("OD600", "GFP55", "GFP100"))

    tidy.data <- long.format.data %>%
        pivot_longer(cols=V1, values_to="Well") %>%
        select(-name) %>% ## drop this useless column
        pivot_longer(!Well,
                     names_to = "Time",
                     values_to = plate.reader.measurement.type,
                     names_pattern = "V(.+)") %>%
        mutate(Time = as.numeric(Time) - 1) %>%
        mutate(minutes = 10*Time) %>%
        mutate(hours = minutes/60) %>%
        mutate(Tet = as.factor(sapply(Well, well.to.Tet))) %>%
        mutate(Plasmid = sapply(Well, well.to.Plasmid)) %>%
        mutate(Population = as.factor(sapply(Well, well.to.population))) %>%
        ## use purrr::map2 black magic to make a column using values in two other columns:
        ## https://github.com/rstudio/cheatsheets/blob/main/purrr.pdf
        ## https://dcl-prog.stanford.edu/purrr-parallel.html
        ## populations 6 and 10 from the pUC treatment lost the pUC plasmid.
        mutate(lost_pUC = map2_lgl(.x=Plasmid, .y=Population,  .f= ~ ifelse(.x=="pUC" && .y %in% c(6,10), TRUE, FALSE))) %>%
        ## we want to encode the lost_pUC state separately.
        mutate(EvolvedPlasmid = ifelse(lost_pUC, "lost pUC", Plasmid)) %>%
        ## odd place to make Plasmid a factor, but need to do it after generating the EvolvedPlasmid column,
        ## so that it gets string values rather than integer factor values.
        mutate(Plasmid = as.factor(Plasmid)) %>%
        mutate(EvolvedPlasmid = as.factor(EvolvedPlasmid)) %>%
    
    return(tidy.data)
}


subtract.LB.blanks <- function(OD600.GFP.tidy.data) {
    ## 1) average blank measurement at every time point.
    ## 2) subtract average blank measurement from each time point.
    blanks <- filter(OD600.GFP.tidy.data, Plasmid == 'LB Blank') %>% group_by(Date) %>%
        summarise(
            blank.OD600.avg = mean(OD600),
            blank.GFP55.avg = mean(GFP55),
            blank.GFP100.avg = mean(GFP100))

    subtracted.OD600.GFP.tidy.data <- OD600.GFP.tidy.data %>%
        left_join(blanks) %>%
        mutate(rawOD600 = OD600) %>%
        mutate(rawGFP55 = GFP55) %>%
        mutate(rawGFP100 = GFP100) %>%
        mutate(OD600 = rawOD600 - blank.OD600.avg) %>%
        mutate(GFP55 = rawGFP55 - blank.GFP55.avg) %>%
        mutate(GFP100 = rawGFP100 - blank.GFP100.avg) %>%
        ## remove the blanks from the polished dataset.
        filter(Plasmid != 'LB Blank') %>%
        ungroup()
    
    return(subtracted.OD600.GFP.tidy.data)
}


tidy.and.merge.OD600.and.GFP.data <- function(
                                              OD600.tetA.GFP.clones.long.format.data,
                                              GFP55.tetA.GFP.clones.long.format.data,
                                              GFP100.tetA.GFP.clones.long.format.data,
                                              datestring,
                                              subtract.LB.blanks=FALSE) {


    ## tidy the data.
    OD600.tidy.data <- tidy.my.tetA.GFP.clone.data(OD600.tetA.GFP.clones.long.format.data, "OD600")
    GFP55.tidy.data <- tidy.my.tetA.GFP.clone.data(GFP55.tetA.GFP.clones.long.format.data, "GFP55")
    GFP100.tidy.data <- tidy.my.tetA.GFP.clone.data(GFP100.tetA.GFP.clones.long.format.data, "GFP100")
    
    ## now merge them together, and add additional columns.
    OD600.GFP.tidy.data <- OD600.tidy.data %>%
        full_join(GFP55.tidy.data) %>%
        full_join(GFP100.tidy.data) %>%
        ## encode the date of data collection in the dataframe.
        ## IMPORTANT: I use the Date column in the subtract.LB.blanks() function!
        mutate(Date = datestring)


    ## subtract LB blanks if set to TRUE.
    if (subtract.LB.blanks) {
        OD600.GFP.tidy.data <- OD600.GFP.tidy.data %>%
            subtract.LB.blanks()        
    }

    ## regardless, remove the LB Blank wells from the data,
    ## and calculate additional columns.
    OD600.GFP.tidy.data <- OD600.GFP.tidy.data %>%
        filter(Plasmid != "LB Blank") %>%
        mutate(normalized.GFP55 = GFP55/OD600) %>%
        mutate(normalized.GFP100 = GFP100/OD600) %>%
        mutate(log.normalized.GFP55 = log(GFP55/OD600)) %>%
        mutate(log.normalized.GFP100 = log(GFP100/OD600))
    
    tidycsvfilename <- paste0(datestring, "-tetA-GFP-evolved-clones-OD600.csv")
    tidycsvpath <- paste0("../results/tetA-GFP-growth-results/", tidycsvfilename)
    ## and write to file.
    write.csv(OD600.GFP.tidy.data, tidycsvpath, row.names=FALSE, quote=FALSE)
    
    return(OD600.GFP.tidy.data)
}


make.timeseries.plots <- function(OD600.GFP.tidy.data) {
    ## make all the timeseries plots.

    datestring <- unique(OD600.GFP.tidy.data$Date)
    
    ## make plots of unnormalized OD600 and GFP.
    OD600.plot <- ggplot(OD600.GFP.tidy.data,
                         aes(x = hours,
                             y = OD600,
                             color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
    OD600.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                              paste0(datestring, "-OD600.pdf"))
    ggsave(filename=OD600.plot.name, OD600.plot, width=10)
    

    GFP55.plot <- ggplot(OD600.GFP.tidy.data,
                         aes(x = hours,
                             y = GFP55,
                             color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
    GFP55.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                              paste0(datestring, "-GFP55.pdf"))
    ggsave(filename=GFP55.plot.name, GFP55.plot, width=10)


    GFP100.plot <- ggplot(OD600.GFP.tidy.data,
                          aes(x = hours,
                              y = GFP100,
                              color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
    GFP100.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                               paste0(datestring, "-GFP100.pdf"))
    ggsave(filename=GFP100.plot.name, GFP100.plot, width=10)


    ## make plots of unnormalized OD600 and GFP on log-scale.
    log.OD600.plot <- ggplot(OD600.GFP.tidy.data,
                             aes(x = hours,
                                 y = log(OD600),
                                 color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
    log.OD600.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                  paste0(datestring, "-log-OD600.pdf"))
    ggsave(filename=log.OD600.plot.name, log.OD600.plot, width=10)


    log.GFP55.plot <- ggplot(OD600.GFP.tidy.data,
                             aes(x = hours,
                                 y = log(GFP55),
                                 color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
    log.GFP55.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                              paste0(datestring, "-log-GFP55.pdf"))
    ggsave(filename=log.GFP55.plot.name, log.GFP55.plot, width=10)


    log.GFP100.plot <- ggplot(OD600.GFP.tidy.data,
                              aes(x = hours,
                                  y = log(GFP100),
                                  color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
     log.GFP100.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                              paste0(datestring, "-log-GFP100.pdf"))
    ggsave(filename=log.GFP55.plot.name, log.GFP100.plot, width=10)

    
    ## make plots of normalized GFP/OD600.
    normalized.GFP55.plot <- ggplot(OD600.GFP.tidy.data,
                                    aes(x = hours,
                                        y = normalized.GFP55,
                                        color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
      normalized.GFP55.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                              paste0(datestring, "-normalized-GFP55.pdf"))
    ggsave(filename=normalized.GFP55.plot.name, normalized.GFP55.plot, width=10)


normalized.GFP100.plot <- ggplot(OD600.GFP.tidy.data,
                             aes(x = hours,
                                 y = normalized.GFP100,
                                 color = EvolvedPlasmid)) +
    geom_point(size=0.5) +
    facet_wrap(Tet~.)
        ## save the plot.
    normalized.GFP100.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                         paste0(datestring, "-normalized-GFP100.pdf"))
    ggsave(filename=normalized.GFP100.plot.name, normalized.GFP100.plot, width=10)
    

## make plots of log-normalized GFP/OD600.
    log.normalized.GFP55.plot <- ggplot(OD600.GFP.tidy.data,
                                        aes(x = hours,
                                            y = log.normalized.GFP55,
                                            color = EvolvedPlasmid)) +
        geom_point(size=0.5) +
        facet_wrap(Tet~.)
    ## save the plot.
    log.normalized.GFP55.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                          paste0(datestring, "-log-normalized-GFP55.pdf"))
    ggsave(filename=log.normalized.GFP55.plot.name, log.normalized.GFP55.plot, width=10)


log.normalized.GFP100.plot <- ggplot(OD600.GFP.tidy.data,
                             aes(x = hours,
                                 y = log.normalized.GFP100,
                                 color = EvolvedPlasmid)) +
    geom_point(size=0.5) +
    facet_wrap(Tet~.)
## save the plot.
    log.normalized.GFP100.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                          paste0(datestring, "-log-normalized-GFP100.pdf"))
    ggsave(filename=log.normalized.GFP100.plot.name, log.normalized.GFP100.plot, width=10)


    ## let's plot the derivatives of log(OD600) over time.
    delta.log2.OD600.data <- OD600.GFP.tidy.data %>%
        group_by(Well, Tet, EvolvedPlasmid, Population) %>%
        mutate(log2.OD600 = log2(OD600)) %>%
        mutate(delta.log2.OD600 = log2.OD600 - lag(log2.OD600))

    delta.log2.OD600.plot <- delta.log2.OD600.data %>%
        group_by(Well) %>%
        filter(hours > 1) %>%
        ggplot(
            aes(x = hours,
                y = delta.log2.OD600,
                color = Well)) + geom_line() + theme_classic() +
        facet_grid(EvolvedPlasmid~Tet) +
        ylab("Delta log2(OD600)") +
        xlab("Time (hours)") +
        guides(color = "none")
    ## save the plot.
    delta.log2.OD600.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                              paste0(datestring, "-delta-OD600.pdf"))
    ggsave(delta.log2.OD600.plot.name,
           delta.log2.OD600.plot, width=10)

    ## return nothing, so that we don't implicitly return the last object from this function.
    return()
}


helper.to.calculate.time.lag <- function(well.df, were.blanks.subtracted = FALSE) {
    ## to compare time lags, measure point when OD600 hits the midpoint of the interval of interest.
    if (were.blanks.subtracted)
        lag.min.index <- max(min(which(well.df$OD600 >= SUBTRACTED.OD600.INTERVAL.MIDPOINT)), 1)
    else
        lag.min.index <- max(min(which(well.df$OD600 >= RAW.OD600.INTERVAL.MIDPOINT)), 1)
        
    t.OD.hit.lag.min <- well.df[lag.min.index,]$hours
    ret.df <- mutate(well.df, time.lag = t.OD.hit.lag.min)
    return(ret.df)
}


calculate.time.lag <- function(OD600.GFP.tidy.data, blanks.subtracted = FALSE) {

    ## generate a time lag calculation function depending on whether blanks were subtracted or not.
    time.lag.helper.function <- partial(helper.to.calculate.time.lag, were.blanks.subtracted = blanks.subtracted)
    
    ## calculate the time to reach the midpoint of the OD600 interval of interest.
    OD600.GFP.tidy.data %>%
        split(.$Well) %>%
        map_dfr(.f=time.lag.helper.function) %>%
        group_by(Population, EvolvedPlasmid) %>%
        ## clunky hack to so that time.lag is not shown for every single time point!
        ## TODO: come up with a simpler solution.
        summarize(mean.time.lag = mean(time.lag))
}


calculate.max.growth.rate <- function(OD600.GFP.tidy.data, blanks.subtracted = FALSE) {
    ## calculate the max growth rate (look at points in the OD600 interval of interest).
    if (blanks.subtracted) {
        tetA.GFP.evolved.clones.max.growth.rate.data <- OD600.GFP.tidy.data %>%
            filter(OD600 >= SUBTRACTED.OD600.INTERVAL.MIN) %>%
            filter(OD600 <= SUBTRACTED.OD600.INTERVAL.MAX)
    } else {
        tetA.GFP.evolved.clones.max.growth.rate.data <- OD600.GFP.tidy.data %>%
            filter(OD600 >= RAW.OD600.INTERVAL.MIN) %>%
            filter(OD600 <= RAW.OD600.INTERVAL.MAX)
    }

    tetA.GFP.evolved.clones.max.growth.rate.data <- tetA.GFP.evolved.clones.max.growth.rate.data %>%
        group_by(Well, Tet, EvolvedPlasmid, Population) %>%
        nest() %>% ## nest the data by well, then calculate the log2(OD) slope using lm().
        mutate(fit = map(data, ~lm(log2(OD600) ~ hours, data = .))) %>%
        ## get the fit parameters.
        ## see tutorial at:
        ## https://www.kaylinpavlik.com/linear-regression-with-nested-data/
        mutate(tidied = map(fit, tidy)) %>%
        select(Well, Tet, EvolvedPlasmid, Population, tidied) %>%
        unnest(tidied) %>%
        filter(term == 'hours') ## we only care about the slope parameter.

    return(tetA.GFP.evolved.clones.max.growth.rate.data)
}


make.growth.parameter.plots <- function(OD600.GFP.tidy.data) {

    datestring <- unique(OD600.GFP.tidy.data$Date)
    
    ## let's plot the distribution of OD at 24h.
    OD600.GFP.tidy.24h.data <- OD600.GFP.tidy.data %>%
        filter(hours == 24)

    tetA.GFP.evolved.clones.24h.OD600.plot <- ggplot(
        OD600.GFP.tidy.24h.data,
        aes(x = EvolvedPlasmid,
            y = OD600,
            color = Population,
            shape=EvolvedPlasmid)) +
        facet_wrap(Tet~.) +
        geom_point() + theme_classic() +
        ylab("OD600 at 24h timepoint") +
        xlab("EvolvedPlasmid presence")
    ##save the plot.
    OD600.at.24h.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                     paste0(datestring, "-evolved-clones-24h-OD600.pdf"))
    ggsave(OD600.at.24h.plot.name,
           tetA.GFP.evolved.clones.24h.OD600.plot, width=10)


    ## now, calculate the max growth rate.
    tetA.GFP.evolved.clones.max.growth.rate.data <- calculate.max.growth.rate(OD600.GFP.tidy.data)
    
    tetA.GFP.evolved.clones.24h.max.growth.rate.plot <- ggplot(
        tetA.GFP.evolved.clones.max.growth.rate.data,
        aes(x = EvolvedPlasmid,
            y = estimate,
            color = Population,
            shape = EvolvedPlasmid)) +
        geom_point() + theme_classic() +
        facet_wrap(Tet~.) +
        ylab("max growth rate") +
        xlab("Presence/absence of target plasmid")
    ## save the plot.
    max.growth.rate.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                        paste0(datestring, "-evolved-clones-growth-rate.pdf"))
    ggsave(max.growth.rate.plot.name,
           tetA.GFP.evolved.clones.24h.max.growth.rate.plot, width=10)


    ## replot these figures as histograms, per Lingchong's request.
    tetA.GFP.evolved.clones.24h.OD600.histogram <- ggplot(
        OD600.GFP.tidy.data,
        aes(fill = EvolvedPlasmid,
            x = OD600)) + geom_histogram(alpha=0.5) + theme_classic() +
        facet_grid(EvolvedPlasmid~Tet) +
        ylab("OD600 at 24h timepoint") +
        labs(fill = "EvolvedPlasmid type") +
        theme(legend.position="top")
    ## save the plot.
    OD600.at.24h.histogram.name <- paste0("../results/tetA-GFP-growth-results/",
                                        paste0(datestring, "-evolved-clones-24h-OD600-histogram.pdf"))
    ggsave(OD600.at.24h.histogram.name,
           tetA.GFP.evolved.clones.24h.OD600.histogram, width=10)

    
    tetA.GFP.evolved.clones.24h.max.growth.rate.histogram <- ggplot(
        tetA.GFP.evolved.clones.max.growth.rate.data,
        aes(fill = EvolvedPlasmid,
            x = estimate)) + geom_histogram(alpha=0.5) + theme_classic() +
        facet_grid(EvolvedPlasmid~Tet) +
        xlab("max growth rate") +
        labs(fill = "EvolvedPlasmid") +
        theme(legend.position="top")
    ## save the plot.
    max.growth.rate.histogram.name <- paste0("../results/tetA-GFP-growth-results/",
                                        paste0(datestring, "-evolved-clones-growth-rate-histogram.pdf"))
    ggsave(max.growth.rate.histogram.name,
           tetA.GFP.evolved.clones.24h.max.growth.rate.histogram, width=10)


    ## now, calculate the time lag.
    tetA.GFP.evolved.clones.tidy.lagtime.data <- calculate.time.lag(OD600.GFP.tidy.data)
    
    tetA.GFP.evolved.clones.timelag.plot <- ggplot(
        tetA.GFP.evolved.clones.tidy.lagtime.data,
        aes(x = EvolvedPlasmid,
            y = mean.time.lag,
            color = Population)) + geom_point() + theme_classic() +
        ylab("Time lag") +
        xlab("EvolvedPlasmid")
    ## save the plot.
    timelag.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                paste0(datestring, "-evolved-clones-OD600-time-lag.pdf"))
    ggsave(timelag.plot.name,
           tetA.GFP.evolved.clones.timelag.plot, width=10)

    ## make a combined plot.
    combined.tetA.GFP.evolved.clones.growth.parameter.plot <- plot_grid(
        tetA.GFP.evolved.clones.timelag.plot,
        tetA.GFP.evolved.clones.24h.max.growth.rate.plot,
        tetA.GFP.evolved.clones.24h.OD600.plot,
        labels=c('A','B','C'), nrow = 1)
    ## save the plot.
    combined.plot.name <- paste0("../results/tetA-GFP-growth-results/",
                                paste0(datestring, "-evolved-clones-combined-growth-plot.pdf"))
    ggsave(combined.plot.name,
       combined.tetA.GFP.evolved.clones.growth.parameter.plot, width = 10)

    ## return nothing, so that we don't implicitly return the last object from this function.
    return()    
}


run.and.print.statistics <- function(OD600.GFP.tidy.data) {
    ## let's calculate some statistics.

    tetA.GFP.evolved.clones.24h.OD600.summary <- OD600.GFP.tidy.data %>%
        group_by(EvolvedPlasmid, Tet) %>%
        summarize(mean_OD600 = mean(OD600), var_OD600 = var(OD600))

    print("Running Wilcox test to compare 24h OD600 between No plasmid and pUC treatments")
    print(wilcox.test(x=filter(tetA.GFP.evolved.clones.24h.OD600.summary,
                         EvolvedPlasmid == 'No plasmid')$mean_OD600,
                y=filter(tetA.GFP.evolved.clones.24h.OD600.summary,
                     EvolvedPlasmid == 'pUC')$mean_OD600))

    print("Running Wilcox test to compare variance in OD600 between No plasmid and pUC treatments")
    print(wilcox.test(x=filter(tetA.GFP.evolved.clones.24h.OD600.summary,
                     EvolvedPlasmid == 'No plasmid')$var_OD600,
            y=filter(tetA.GFP.evolved.clones.24h.OD600.summary,
                     EvolvedPlasmid == 'pUC')$var_OD600))

    tetA.GFP.evolved.clones.max.growth.rate.data <- calculate.max.growth.rate(OD600.GFP.tidy.data)

    tetA.GFP.evolved.clones.max.growth.rate.summary <- tetA.GFP.evolved.clones.max.growth.rate.data %>%
    group_by(Tet, EvolvedPlasmid) %>%
    summarize(mean_growth_rate = mean(estimate), var_growth_rate = var(estimate))

    print("Running Wilcox test to compare max growth rate between No plasmid and pUC treatments")
    print(wilcox.test(x=filter(tetA.GFP.evolved.clones.max.growth.rate.summary,
                     EvolvedPlasmid == 'No plasmid')$mean_growth_rate,
            y=filter(tetA.GFP.evolved.clones.max.growth.rate.summary,
                     EvolvedPlasmid == 'pUC')$mean_growth_rate))

        print("Running Wilcox test to compare variance in max growth rate between No plasmid and pUC treatments")
    print(wilcox.test(x=filter(tetA.GFP.evolved.clones.max.growth.rate.summary,
                     EvolvedPlasmid == 'No plasmid')$var_growth_rate,
            y=filter(tetA.GFP.evolved.clones.max.growth.rate.summary,
                     EvolvedPlasmid == 'pUC')$var_growth_rate))

    tetA.GFP.evolved.clones.tidy.lagtime.data <- calculate.time.lag(OD600.GFP.tidy.data)
    
    no.plasmid.timelag.data <- filter(tetA.GFP.evolved.clones.tidy.lagtime.data,
                                      EvolvedPlasmid == "No plasmid")
    pUC.plasmid.timelag.data <- filter(tetA.GFP.evolved.clones.tidy.lagtime.data,
                                       EvolvedPlasmid == "pUC")
    
    print("Running Wilcox test to compare time lag between No plasmid and pUC treatments")
    print(wilcox.test(no.plasmid.timelag.data$mean.time.lag,
                pUC.plasmid.timelag.data$mean.time.lag))

    ## return nothing, so that we don't implicitly return the last object from this function.
    return()    

}

############################################################################################
## Code to fit data to the logistic growth function.

sigmoid <- function(x) {
    ## See: https://en.wikipedia.org/wiki/Logistic_function
    ## sigmoid is the standard logistic function (L = 1, k = 1, x0 = 0)
    return(1 / (1 + exp(-x)))
}


logistic.growth.function <- function(t, k_max, tau, L) {
    ## In the fitted function, the growth rate is represented by k_max (the maximal exponential rate),
    ## lag time by tau (the time before maximum growth is reached),
    ## yield by L (the maximal optical density reached),
    return(L*sigmoid(k_max*(t-tau)))
}


fit.logistic.growth.function.on.one.well <- function(well.df) {

    ## get the metadata for the input growth data.
    well <- unique(well.df$Well)
    tet <- unique(well.df$Tet)
    plasmid <- unique(well.df$Plasmid)
    population <- unique(well.df$Population)
    lost.pUC <- unique(well.df$lost_pUC)
    evolved.plasmid <- unique(well.df$EvolvedPlasmid)
    date <- unique(well.df$Date)
    preconditioning.tet <- unique(well.df$PreconditioningTet)
    preconditioning.length <- unique(well.df$preconditioning_length)

    ## assert that well.df has a unique combination of these metadata-- only single values allowed!
    stopifnot(all(length(well) == length(tet), 
                  length(tet) == length(plasmid), 
                  length(plasmid) == length(population), 
                  length(population) == length(lost.pUC), 
                  length(lost.pUC) == length(evolved.plasmid), 
                  length(evolved.plasmid) == length(date), 
                  length(date) == length(preconditioning.tet),
                  length(preconditioning.tet) == length(preconditioning.length),
                  length(preconditioning.length) == 1))

    ## get the growth time series data for the well.
    t_data <- well.df$Time
    y_data <- well.df$OD600

    ## initial guesses for the parameters.
    initial_params_set <- list(k_max = 0.05, tau = 50, L = 1)

    ## For debugging nls fits that fail:
    cat("Data collection date:", unique(well.df$Date), "\n")
    cat("Well:", unique(well.df$Well), "\n")
    
    ## Fit the data to the logistic growth function
    fit <- nls(y_data ~ logistic.growth.function(t_data, k_max, tau, L),
               start = initial_params_set)
  
    ## Calculate R-squared
    fitted_values <- predict(fit)
    residuals <- y_data - fitted_values
    r_squared <- 1 - sum(residuals^2) / sum((y_data - mean(y_data))^2)

    ## Extract the best-fitted parameters
    best_k_max_fit <- coef(fit)["k_max"]
    best_tau_fit <- coef(fit)["tau"]
    best_L_fit <- coef(fit)["L"]

    ## Print the best-fitted parameters and R-squared
    cat("Fitted k_max:", best_k_max_fit, "\n")
    cat("Fitted tau:", best_tau_fit, "\n")
    cat("Fitted L:", best_L_fit, "\n")
    cat("R-squared:", r_squared, "\n")

    ## return the results as a dataframe.
    well.results.df <- data.frame(
        Well = well,
        Tet = tet,
        Plasmid = plasmid,
        Population = population,
        lost_pUC = lost.pUC,
        EvolvedPlasmid = evolved.plasmid,
        Date = date,
        PreconditioningTet = preconditioning.tet,
        preconditioning_length = preconditioning.length,
        max_growth_rate = best_k_max_fit,
        lag_time = best_tau_fit,
        yield = best_L_fit,
        Rsquared = r_squared,
        row.names = NULL
    )

    return(well.results.df)
}


########################################################################
## Data analysis.

## I split the original data from Tecan J into three separate files, for OD600, GFP55, and GFP100 readings.
## each contains a 30h time series of the 30 tetA-GFP clones 

## LB preconditioned data, Block 1.
August8.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-08-LBandLBTet50-tetA-GFP-evolved-clones-OD600.csv", header=FALSE)
August8.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-08-LBandLBTet50-tetA-GFP-evolved-clones-GFP55.csv", header=FALSE)
August8.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-08-LBandLBTet50-tetA-GFP-evolved-clones-GFP100.csv", header=FALSE)

## LB preconditioned data, Block 2.
August16.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-16-LBandLBTet50-tetA-GFP-evolved-clones-OD600.csv", header=FALSE)
August16.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-16-LBandLBTet50-tetA-GFP-evolved-clones-GFP55.csv", header=FALSE)
August16.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-16-LBandLBTet50-tetA-GFP-evolved-clones-GFP100.csv", header=FALSE)

## LB preconditioned data, Block 3.
August19.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-19-LBandLBTet50-tetA-GFP-evolved-clones-OD600.csv", header=FALSE)
August19.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-19-LBandLBTet50-tetA-GFP-evolved-clones-GFP55.csv", header=FALSE)
August19.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-19-LBandLBTet50-tetA-GFP-evolved-clones-GFP100.csv", header=FALSE)

## LB+Tet50 preconditioned data, Block 1.
August21.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-21-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-OD600.csv", header=FALSE)
August21.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-21-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-GFP55.csv", header=FALSE)
August21.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-21-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-GFP100.csv", header=FALSE)

## LB+Tet50 preconditioned data, Block 2.
## CRITICAL NOTE: Not really planned--
## I let the cultures grow into stationary phase for an extended period-- I started
## growth curves ~10pm, instead of 5pm.
## So far, this is the only growth data in which I see long lags of the no plasmid clones.

## Question: does the extended lag seen for some no plasmid clones in these data relate to growth
## in stationary phase for an extended period of time?

## Is there some causal variable I can control affecting the appearance of extended lags in the no plasmid
## clones, or is this largely a stochastic factor??

## I can test this question by starting growth curves after exactly 24h, and after letting the cultures sit
## for a longer period, say 30, 36, or 48h.

August23.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-23-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-OD600.csv", header=FALSE)
August23.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-23-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-GFP55.csv", header=FALSE)
August23.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-23-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-GFP100.csv", header=FALSE)


## LB+Tet50 preconditioned data, Block 3.
August29.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-29-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-OD600.csv", header=FALSE)
August29.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-29-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-GFP55.csv", header=FALSE)
August29.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-29-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-GFP100.csv", header=FALSE)

#################################
## LB+Tet50 preconditioned, after ~48h of growth at 37C. 90/180 cycles completed.
## These are the same cultures as used for the batch 3 LB+Tet50 cultures (August 29),
## but after 48h growth, instead of ~24h growth.
## My prediction is that I will see extended lags in these data. This could be caused
## mechanistically by TetA degradation during stationary phase, such that the no plasmid
## strains are more likely to have not enough TetA to start growing immediately in fresh LB+Tet50 media.

## I see extended lags in wells F7, and F11, which are no plasmid LB+Tet50 pops 6 and 10.
## These are the same clones showing extended lags in the Batch 2 data.
## Maybe the right comparison is lags for the no plasmid clones and the p15A clones.

August30.OD600.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-30-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-after48h-growth-OD600.csv", header=FALSE)
August30.GFP55.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-30-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-after48h-growth-GFP55.csv", header=FALSE)
August30.GFP100.tetA.GFP.clones.long.format.data <- read.csv("../data/tetA-GFP-OD600-data/2023-08-30-LBandLBTet50-tetA-GFP-evolved-clones-Tet50-preconditioning-after48h-growth-GFP100.csv", header=FALSE)


## Now tidy and merge the OD600 and GFP datasets.
August8.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August8.OD600.tetA.GFP.clones.long.format.data,
    August8.GFP55.tetA.GFP.clones.long.format.data,
    August8.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-8") %>%
    mutate(PreconditioningTet = 0) %>%
    mutate(preconditioning_length = "24h")

August16.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August16.OD600.tetA.GFP.clones.long.format.data,
    August16.GFP55.tetA.GFP.clones.long.format.data,
    August16.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-16") %>%
    mutate(PreconditioningTet = 0) %>%
    mutate(preconditioning_length = "24h")

August19.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August19.OD600.tetA.GFP.clones.long.format.data,
    August19.GFP55.tetA.GFP.clones.long.format.data,
    August19.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-19") %>%
    mutate(PreconditioningTet = 0) %>%
    mutate(preconditioning_length = "24h")

August21.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August21.OD600.tetA.GFP.clones.long.format.data,
    August21.GFP55.tetA.GFP.clones.long.format.data,
    August21.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-21") %>%
    mutate(PreconditioningTet = 50) %>%
    mutate(preconditioning_length = "24h")

August23.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August23.OD600.tetA.GFP.clones.long.format.data,
    August23.GFP55.tetA.GFP.clones.long.format.data,
    August23.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-23") %>%
    mutate(PreconditioningTet = 50) %>%
    mutate(preconditioning_length = "24h")

August29.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August29.OD600.tetA.GFP.clones.long.format.data,
    August29.GFP55.tetA.GFP.clones.long.format.data,
    August29.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-29") %>%
    mutate(PreconditioningTet = 50) %>%
    mutate(preconditioning_length = "24h")

## 48h pre-conditioning growth data.
August30.OD600.GFP.tidy.data <- tidy.and.merge.OD600.and.GFP.data(
    August30.OD600.tetA.GFP.clones.long.format.data,
    August30.GFP55.tetA.GFP.clones.long.format.data,
    August30.GFP100.tetA.GFP.clones.long.format.data,
    "2023-8-30") %>%
    mutate(PreconditioningTet = 50) %>%
    mutate(preconditioning_length = "48h")

## join all these data.
all.tetA.GFP.clone.growth.data <- August8.OD600.GFP.tidy.data %>%
    full_join(August16.OD600.GFP.tidy.data) %>%
    full_join(August19.OD600.GFP.tidy.data) %>%
    full_join(August21.OD600.GFP.tidy.data) %>%
    full_join(August23.OD600.GFP.tidy.data) %>%
    full_join(August29.OD600.GFP.tidy.data) %>%
    full_join(August30.OD600.GFP.tidy.data)

## write these data to file.
write.csv(all.tetA.GFP.clone.growth.data, "../results/tetA-GFP-growth-results/August2023-tetA-GFP-clone-growth-curves.csv", row.names=FALSE, quote=FALSE)

##############################################
## Make timeseries plots.
all.tetA.GFP.clone.growth.data %>%
    split(.$Date) %>%
    map_dfr(.f = make.timeseries.plots)

## make growth parameter plots:
all.tetA.GFP.clone.growth.data %>%
    split(.$Date) %>%
    map_dfr(.f = make.growth.parameter.plots)

## Run some simple statistics comparing the no plasmid and pUC plasmid populations.
all.tetA.GFP.clone.growth.data %>%
    split(.$Date) %>%
    map_dfr(.f = run.and.print.statistics)

## Fit the logistic function to the data.
logistic.growth.parameters.df <- all.tetA.GFP.clone.growth.data %>%
    split(list(.$Date, .$Well)) %>%
    map_dfr(.f = fit.logistic.growth.function.on.one.well)


## let's plot the growth parameters.
max.growth.rate.plot <- ggplot(
    data = logistic.growth.parameters.df,
    aes(x = Well, y = max_growth_rate, color = EvolvedPlasmid)) +
    geom_point() +
    theme_classic() +
    facet_wrap(PreconditioningTet~Tet)


max.growth.rate.plot

lagtime.plot <- ggplot(
    data = logistic.growth.parameters.df,
    aes(x = Well, y = lag_time, color = EvolvedPlasmid)) +
    geom_point() +
    theme_classic() +
    facet_wrap(PreconditioningTet~Tet)


lagtime.plot

