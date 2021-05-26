#' Function to process
#'
#' @param inputLocation, the full pathname of the data directory containing donor directories (named after donors)
#  The inputLoc directory should contain one subdirectory per donor, plus a WellData.txt and isotype.txt file.
#' @param outputLocation, pathname where generated figures will be created
#' @return create figures in outputLocation
kinetxProcess <-
  function(inputLoc,
           outputLoc,
           measure) {
    library(minpack.lm)
    library(psych)

    #Extract donor numbers
    donorDirList <-
      list.dirs(inputLoc,
                full.names = FALSE,
                recursive = FALSE)

    #Open the WellData file, this lists all of the files, with associated agonists and concentrations
    wellData <- NULL
    try(wellData <-
          read.csv(
            paste(inputLoc, "WellData.txt", sep = ""),
            stringsAsFactors = FALSE,
            encoding = "UTF-8"
          ),
        silent = TRUE)
    if (is.null(wellData)) {
      stop("no wellData.txt")
    }
    
    #function or calcium
    refs <-
      wellData[!duplicated(wellData$ref), c("ref")]
    #refs <- subset(refs,ref != "NULL")

    #Open the isotypeData file, this tells us which isytope to use, this was run using an antibody that binds to nothing.
    isoypeData <- NULL
    try(isotypeData <-
          read.csv(
            paste(inputLoc, "IsotypeData.txt", sep = ""),
            stringsAsFactors = FALSE,
            encoding = "UTF-8"
          ),
        silent = TRUE)
    if (is.null(isotypeData)) {
      stop("no isotypeData.txt")
    }

    isotypeList <-
      isotypeData$isotypes
    
    #create structure for the summary file
    summaryFile <-
      matrix(nrow = length(donorDirList) * 15,
             ncol = 39)
    summaryFileCount <- 1
    colnames(summaryFile) <-
      c(
        "donor",
        "output",
        "measure",
        "agonist",
        "conc",
        "drug",
        "0secs",
        "10secs",
        "20secs",
        "30secs",
        "40secs",
        "60secs",
        "90secs",
        "120secs",
        "180secs",
        "10secRoC",
        "20secRoC",
        "30secRoC",
        "40secRoC",
        "60secRoC",
        "90secRoC",
        "120secRoC",
        "180secRoC",
        "isotypeMedian",
        "zeroPoint",
        "maxRateChange",
        "maxRateChangePosition",
        "minRateChange",
        "minRateChangePosition",
        "maxRateAccel",
        "maxRateAccelPosition",
        "minRateAccel",
        "minRateAccelPosition",
        "Shape",
        "ShapeMetric",
        "RateLabel",
        "RoC",
        "EarlyRateLabel",
        "EarlyRateMetric"
      )

    #Loop through each donor's directory
    for (i in 1:length(donorDirList)) {
      carryon <- 0
      if (measure == "FL1") {
        if (file.exists(paste(inputLoc, donorDirList[i], "/", "D3.csv", sep = ""))) {
          dfIsotype <-
            read.csv(
              paste(inputLoc, donorDirList[i], "/", "D3.csv", sep = ""),
              stringsAsFactors = FALSE,
              header = TRUE
            )
          dfIsotype <- dfIsotype[, c("FL1.A", "Time")]
          names(dfIsotype) <- c("fl", "timePoints")
        } else {
          carryon <- 1
        }
      } else {
        if (file.exists(paste(inputLoc, donorDirList[i], "/", "D4.csv", sep = ""))) {
          dfIsotype <-
            read.csv(
              paste(inputLoc, donorDirList[i], "/", "D4.csv", sep = ""),
              stringsAsFactors = FALSE,
              header = TRUE
            )
          dfIsotype <- dfIsotype[, c("FL4.A", "Time")]
          names(dfIsotype) <- c("fl", "timePoints")
        } else {
          carryon <- 1
        }
      }

      if (carryon == 0) {
        #process the files for each donor specified in wellData.txt (there is one file per agonist/concentration)
        for (j in 1:length(refs)) {
          if (measure == "FL1" | (measure == "FL4" & refs[j] != "calcium")) {
            wellDatas <-
              wellData[which(wellData$ref == refs[j]), c("well", "agonist", "conc", "units")]
            ThirtyMinWellDatas <-
              wellData[which(wellData$ref == refs[j]), "endpointWell"]

            print(
              paste(
                "processing: ",
                inputLoc,
                donorDirList[i],
                ",",
                refs[j],
                ",",
                measure,
                ".csv",
                sep = ""
              )
            )

            #Set up output plot
            OutF <-
              paste(outputLoc,
                    donorDirList[i],
                    "_",
                    refs[j],
                    "_",
                    measure,
                    ".png",
                    sep = "")
            png(
              filename = OutF,
              width = 12,
              height = 11,
              units = 'in',
              res = 400
            )

            par(mfrow = c(5, 4))
            par(oma = c(0.2, 0, 1.5, 0.5))
            par(mgp = c(1.2, 0.3, 0))
            par(mar = c(5.4, 2, 1, 1.75))
            par(
              ps = 12,
              cex = 1,
              cex.main = 0.8,
              cex.axis = 0.8,
              cex.lab = 0.8,
              font.main = 1
            )

            #Process individual files for this donor, 1 by 1
            for (k in 1:nrow(wellDatas)) {
              destfile <-
                paste(inputLoc,
                      donorDirList[i],
                      "/",
                      wellDatas[k, "well"],
                      ".csv",
                      sep = "")

              if (file.exists(destfile)) {
                df <-
                  read.csv(destfile,
                           stringsAsFactors = FALSE,
                           header = TRUE)

                if (measure == "FL1") {
                  df.raw <- df[, c("FL1.A", "Time")]
                  names(df.raw) <- c("fl", "timePoints")
                }
                if (measure == "FL4") {
                  df.raw <- df[, c("FL4.A", "Time")]
                  names(df.raw) <- c("fl", "timePoints")
                }

                #Remove anything greater than 200000 and less than 1 - we presume anything out of this range is nonesense
                df.raw <- subset(df.raw, fl < 200000)
                df.raw <- subset(df.raw, fl > 1)

                #Exclude early timepoints, before an agonist was applied
                df.raw <- subset(df.raw, timePoints > 280)

                #Reset all data to be the distance above the Isotype
                df.raw$fl <- df.raw$fl - median(dfIsotype$fl)

                #Grab hold of stats for the raw data - will use the median in the following plot
                df.psych <-
                  data.frame(describeBy(df.raw$fl, df.raw$timePoints, mat = TRUE))

                #Local Regression (fitting), smooths the data using points in a neighbourhood (controlled by span) and weighted
                #by their distance from a point (x)
                loe <-
                  loess(df.raw$fl ~ df.raw$timePoints,
                        span = 0.1,
                        degree = 1)

                df.smooth <-
                  data.frame(seq(from = 281, to = 2000), predict(loe, seq(
                    from = 281, to = 2000
                  )))
                names(df.smooth) <- c("timePoints", "fl")

                #Now I'm generating an EXTRA smoothed line
                dfZero <-
                  subset(df.raw, timePoints > 280 &
                           timePoints < 320)
                if (nrow(dfZero) == 0) {
                  zeroPoint <- 0
                } else {
                  zeroPoint <-
                    sum(median(dfZero$fl)) / length(median(dfZero$fl))
                }
                df.raw.zero <- df.raw
                df.raw.zero$fl <- df.raw.zero$fl - zeroPoint

                loe <-
                  loess(
                    df.raw.zero$fl ~ df.raw.zero$timePoints,
                    span = 0.5,
                    degree = 1
                  )
                df.smooth.extra <-
                  data.frame(seq(from = 281, to = 2000), predict(loe, seq(
                    from = 281, to = 2000
                  )))
                names(df.smooth.extra) <- c("timePoints", "fl")

                #Calculate rate of change of the smoothed data
                smoothPoints <- 200 #
                dY <-
                  diff(df.smooth.extra$fl,
                       lag = smoothPoints)
                dX <- df.smooth.extra$timePoints + smoothPoints / 2
                chopdX <- length(dX) - smoothPoints #
                dX <- dX[1:chopdX]
                rateOfChange <-
                  round(dY / smoothPoints, digits = 2) #so x distance of differential = smoothPoints

                maxRateOfChange <- max(rateOfChange, na.rm = TRUE)
                whereMaxRateOfChange <- which.max(rateOfChange)
                minRateOfChange <- min(rateOfChange, na.rm = TRUE)
                whereMinRateOfChange <- which.min(rateOfChange)

                dYY <-
                  diff(df.smooth.extra$fl,
                       lag = smoothPoints,
                       differences = 2)
                dXX <- dX + smoothPoints / 2
                chopdXX <- length(dXX) - smoothPoints #
                dXX <- dXX[1:chopdXX]
                rateOfAccel <-
                  round(dYY / smoothPoints, digits = 2) #so x distance of differential = smoothPoints

                maxRateOfAccel <- max(rateOfAccel, na.rm = TRUE)
                whereMaxRateOfAccel <- which.max(rateOfAccel)
                minRateOfAccel <- min(rateOfAccel, na.rm = TRUE)
                whereMinRateOfAccel <- which.min(rateOfAccel)

                pred0 <- predict(loe, 281)
                pred10 <- predict(loe, 281 + 100)
                pred20 <- predict(loe, 281 + 200)
                pred30 <- predict(loe, 281 + 300)
                pred40 <- predict(loe, 281 + 400)
                pred60 <- predict(loe, 281 + 600)
                pred90 <- predict(loe, 281 + 900)
                pred120 <- predict(loe, 281 + 1200)
                if (is.na(pred120)) {
                  pred120 <- 0
                }
                pred180 <- predict(loe, 281 + 1800)
                RoC10 <- round((pred10 - pred0) / 10, 2)
                RoC20 <- round((pred20 - pred10) / 10, 2)
                RoC30 <- round((pred30 - pred20) / 10, 2)
                RoC40 <- round((pred40 - pred30) / 10, 2)
                RoC60 <- round((pred60 - pred40) / 20, 2)
                RoC90 <- round((pred90 - pred60) / 30, 2)
                RoC120 <- round((pred120 - pred90) / 30, 2)
                RoC180 <- round((pred180 - pred120) / 60, 2)
                RoC <- round((pred120 - pred0) / 120, 2)

                Label <- "NA"
                LabelMetric <- "NA"
                
                if (refs[j] != "calcium") {
                  
                  if (measure == "FL1" & pred120 < 4000) {
                    plotLabel <- "Y"
                    Label <- "low responder"
                    LabelMetric <- pred120
                  }
                  if (measure == "FL4" & pred120 < 1000) {
                    plotLabel <- "Y"
                    Label <- "low responder"
                    LabelMetric <- pred120
                  }
                  #If not low responder
                  if (Label == "NA") {
                    if (!is.na(RoC120) & !is.na(RoC10) & !is.na(RoC20)) {
                      LabelMetric <- RoC120 / ((RoC10 + RoC20) / 2)
                      if (RoC120 / ((RoC10 + RoC20) / 2) > 2) {
                        Label <- "increasing"
                      } else if (RoC120 / ((RoC10 + RoC20) / 2) < 0.9) {
                        Label <- "decreasing"
                      } else {
                        Label <- "linear"
                      }
                    }
                  }
                  LabelRate <- "NA"
                  if (!is.na(RoC)) {
                      if (measure == "FL1" & RoC > 80) {
                          LabelRate <- "fast"
                        }
                        if (measure == "FL4" & RoC > 10) {
                          LabelRate <- "fast"
                        }
                        if (measure == "FL1" & RoC <= 80) {
                          LabelRate <- "medium"
                        }
                        if (measure == "FL4" & RoC <= 10) {
                          LabelRate <- "medium"
                        }
                        if (measure == "FL1" & RoC <= 40) {
                          LabelRate <- "slow"
                        }
                        if (measure == "FL4" & RoC <= 5) {
                          LabelRate <- "slow"
                        }
                      }
                    
                  EarlyRate <- "NA"
                  EarlyRateMetric <- "NA"
                  if (!is.na(RoC10) & !is.na(RoC20)) {
                    EarlyRateMetric <- (RoC10 + RoC20) / 2
                    if (measure == "FL1" &
                        ((RoC10 + RoC20) / 2 > 80)) {
                      EarlyRate <- "fast"
                    }
                    if (measure == "FL4" &
                        ((RoC10 + RoC20) / 2 > 12)) {
                      EarlyRate <- "fast"
                    }
                    if (measure == "FL1" &
                        ((RoC10 + RoC20) / 2 <= 80)) {
                      EarlyRate <- "medium"
                    }
                    if (measure == "FL4" &
                        ((RoC10 + RoC20) / 2 <= 12)) {
                      EarlyRate <- "medium"
                    }
                    if (measure == "FL1" &
                        ((RoC10 + RoC20) / 2 <= 40)) {
                      EarlyRate <- "slow"
                    }
                    if (measure == "FL4" &
                        ((RoC10 + RoC20) / 2 <= 4)) {
                      EarlyRate <- "slow"
                    }
                  }

                  #Now lets plot it all
                  mainTitle = paste(wellDatas[k, "conc"], wellDatas[k, "units"], " ", wellDatas[k, "agonist"])

                  if (measure == "FL1") {
                    ylimit <- c(1, 800000)
                  } else {
                    ylimit <- c(1, 80000)
                  }
                  #Show the raw data on a log scale
                  plot(
                    df.raw$timePoints / 10,
                    df.raw$fl,
                    ylim = ylimit,
                    log = "y",
                    cex = 0.2,
                    col = "grey70",
                    xlab = "time (secs)",
                    ylab = "log(flurescence)",
                    main = mainTitle
                  )
                  #plot the median of the data
                  lines(
                    as.numeric(as.vector(df.psych$group1)) / 10,
                    df.psych$median,
                    log = "y",
                    col = "grey40",
                    lwd = 2
                  )

                  #Plot the smoothed line
                  lines(
                    df.smooth$timePoints / 10,
                    df.smooth$fl,
                    log = "y",
                    col = "green",
                    lwd = 2
                  )

                  #Now show the smoothed line and the rate of change (not on a log scale)
                  if (measure == "FL1") {
                    ylimit <- c(1, 40000)
                  } else {
                    ylimit <- c(1, 10000)
                  }
                  if (refs[j] == "calcium") {
                    xlimit <- c(0, 140)
                  } else {
                    xlimit <- c(0, 200)
                  }

                  #Show the extra smoothed line - not on a log scale
                  plot(
                    df.smooth.extra$timePoints / 10,
                    df.smooth.extra$fl,
                    col = "green",
                    lwd = 2,
                    ylim = ylimit,
                    xlim = xlimit,
                    cex = 0.2,
                    xlab = "time (secs)",
                    ylab = "flurescence",
                    main = mainTitle
                  )

                  if (measure == "FL1") {
                    text(100, 34000, Label,
                         cex = 1.5)
                    text(100, 26000, LabelRate,
                         cex = 1.5)
                  } else {
                    text(100, 8000, Label,
                         cex = 1.5)
                    text(100, 6800, LabelRate,
                         cex = 1.5)
                  }

                  #Add some text
                  if (measure == "FL1") {
                    mtext(
                      paste(donorDirList[i], refs[j], " fibrinogen"),
                      outer = TRUE,
                      cex = 1
                    )
                  } else {
                    mtext(
                      paste(donorDirList[i], refs[j], " P-selectin"),
                      outer = TRUE,
                      cex = 1
                    )
                  }

                  #Add a line to the summary file
                  summaryFile[summaryFileCount, 1] <-
                    donorDirList[i]
                  summaryFile[summaryFileCount, 2] <- refs[j]
                  summaryFile[summaryFileCount, 3] <- measure
                  summaryFile[summaryFileCount, 4] <-
                    wellDatas[k, "agonist"]
                  summaryFile[summaryFileCount, 5] <-
                    wellDatas[k, "conc"]
                  summaryFile[summaryFileCount, 6] <- "none"
                  summaryFile[summaryFileCount, 7] <- pred0
                  summaryFile[summaryFileCount, 8] <- pred10
                  summaryFile[summaryFileCount, 9] <- pred20
                  summaryFile[summaryFileCount, 10] <- pred30
                  summaryFile[summaryFileCount, 11] <- pred40
                  summaryFile[summaryFileCount, 12] <- pred60
                  summaryFile[summaryFileCount, 13] <- pred90
                  summaryFile[summaryFileCount, 14] <- pred120
                  summaryFile[summaryFileCount, 15] <- pred180
                  summaryFile[summaryFileCount, 16] <- RoC10
                  summaryFile[summaryFileCount, 17] <- RoC20
                  summaryFile[summaryFileCount, 18] <- RoC30
                  summaryFile[summaryFileCount, 19] <- RoC40
                  summaryFile[summaryFileCount, 20] <- RoC60
                  summaryFile[summaryFileCount, 21] <- RoC90
                  summaryFile[summaryFileCount, 22] <- RoC120
                  summaryFile[summaryFileCount, 23] <- RoC180
                  summaryFile[summaryFileCount, 24] <-
                    median(dfIsotype$fl)
                  summaryFile[summaryFileCount, 25] <-
                    zeroPoint
                  try(summaryFile[summaryFileCount, 26] <-
                        maxRateOfChange,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 27] <-
                        whereMaxRateOfChange / 10,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 28] <-
                        minRateOfChange,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 29] <-
                        whereMinRateOfChange / 10,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 30] <-
                        maxRateOfAccel,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 31] <-
                        whereMaxRateOfAccel / 10,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 32] <-
                        minRateOfAccel,
                      silent = TRUE)
                  try(summaryFile[summaryFileCount, 33] <-
                        whereMinRateOfAccel / 10,
                      silent = TRUE)
                  summaryFile[summaryFileCount, 34] <-
                    Label
                  summaryFile[summaryFileCount, 35] <-
                    LabelMetric
                  summaryFile[summaryFileCount, 36] <-
                    LabelRate
                  summaryFile[summaryFileCount, 37] <-
                    RoC
                  summaryFile[summaryFileCount, 38] <-
                    EarlyRate
                  summaryFile[summaryFileCount, 39] <-
                    EarlyRateMetric
                  summaryFileCount <- summaryFileCount + 1

                } else {
                  print(paste(destfile, " does not exist"))
                }

              }
            }

            dev.off()

          }
        }
      }
    }

    #Save the summary
    summaryFileName <-
      paste(outputLoc, "summary", measure, "FileV1.csv", sep = "")

    write.csv(summaryFile, file = summaryFileName, row.names = FALSE)

    print("finished")

  }

