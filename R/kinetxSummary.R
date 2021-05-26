#' Function to process
#'
#' @param inputFile, location and name of summary file to work on
#' @param outputLocation, where to put the results
#' @param meas, the measure to process e.g. FL1
#' @return create a new summary file V2 in outputLocation
kinetxSummary <-
  function(inputFile,
           outputLoc,
           meas) {

    library(psych)
    library(dplyr)
    library(ggplot2)
    library(tidyverse)

    df <-
      read.csv(inputFile,
               stringsAsFactors = FALSE,
               header = TRUE)

    df <- subset(df, measure == meas)

    df <- df[, c("donor", "output", "measure", "agonist", "RoC")]
    df <- subset(df, donor != "MET053A")
    df <- subset(df, donor != "MET052A")
    df <- subset(df, donor != "MET066A")
    df <- subset(df, donor != "MET020A")

    #This code presumes there is only 1 output!
    outputs <-
      df[!duplicated(df$output), c("output")]

    agonists <-
      df[!duplicated(df$agonist), c("agonist")]

    donors <-
      df[!duplicated(df$donor), c("donor")]

    df <- df %>% drop_na("RoC")

    #Create some summary stats for RoC
    df.aggregates.max <-
      aggregate(df$RoC,
                by = df["agonist"],
                FUN = max,
                na.rm = TRUE)
    df.aggregates.min <-
      aggregate(
        df$RoC,
        by = df["agonist"],
        FUN = min,
        na.rm = TRUE,
        na.action = "na.pass"
      )

    names(df.aggregates.max) <- c("agonist", "maxRoC")
    names(df.aggregates.min) <- c("agonist", "minRoC")

    df.aggregates.max$maxRoC <-
      as.numeric(df.aggregates.max$maxRoC)
    df.aggregates.min$minRoC <-
      as.numeric(df.aggregates.min$minRoC)
    df.new <- left_join(df, df.aggregates.max, by = "agonist")
    df.new <- left_join(df.new, df.aggregates.min, by = "agonist")

    #Scale over each agonist
    df.new$agonistMetric <- df.new$RoC / df.new$maxRoC

    #Sum the agonist RoC rescaled for each donor
    df.aggregates.kinetxMetric <-
      aggregate(df.new$agonistMetric, by = df["donor"], FUN = sum)
    names(df.aggregates.kinetxMetric) <-
      c("donor", "kinetxMetric")

    #Divide by the number of agonists
    df.aggregates.count <-
      aggregate(df.new$RoC, by = df["donor"], FUN = length)
    names(df.aggregates.count) <- c("donor", "length")

    df.aggregates.count$length <-
      as.numeric(df.aggregates.count$length)

    #Create a new simplified data frame with just donor, output, measure and kinetxMetric
    df.aggregates <-
      merge(df.aggregates.count, df.aggregates.kinetxMetric, by = "donor")

    df.aggregates$kinetxMetric <-
      df.aggregates$kinetxMetric / df.aggregates$length

    df.new <- unique(df.new[, 1:3])

    df.new <-
      merge(df.new, df.aggregates, by = "donor")


    write.csv(
      df.new,
      file = paste(outputLoc, "summary", meas, "FileV2.csv", sep = ""),
      row.names = FALSE
    )

    print("finished")

  }
