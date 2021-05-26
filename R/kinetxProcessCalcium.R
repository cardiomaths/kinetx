#' Function to process
#'
#' @param inputLocation, the full pathname of the data directory containing donor directories (named after donors)
#  The inputLoc directory should contain one subdirectory per donor, plus a WellData.txt and isotype.txt file.
#' @param outputLocation, pathname where generated figures are created
#' @return create figures in outputLocation
kinetxProcessCalcium <-
  function(inputLoc,
           outputLoc,
           measure,
           donor='All') {
    
    library(minpack.lm)
    library(psych)

    #This provides a list of all the donors
    donorDirList <-
      list.dirs(inputLoc,
                full.names = FALSE,
                recursive = FALSE)

    if (donor != "ALL") {
      donorDirList <- donorDirList[grep(donor, donorDirList)]
    }

    #Open the WellData file, this lists all the donors files that we wish to process
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
    wellData <- subset(wellData, ref == "calcium")

    #function or calcium
    refs <-
      wellData[!duplicated(wellData$ref), c("ref")]

    #create structure for the summary file
    summaryFile <-
      matrix(nrow = length(donorDirList) * 15,
             ncol = 37)
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
        "Drop",
        "AbsoluteDrop",
        "maxRateChange",
        "maxRateChangePosition",
        "minRateChange",
        "minRateChangePosition",
        "maxRateAccel",
        "maxRateAccelPosition",
        "minRateAccel",
        "minRateAccelPosition",
        "category",
        "categoryMetric",
        "Rate",
        "RateMetric"
      )

    #Loop through each donor's directory
    for (i in 1:length(donorDirList)) {
      wellDatas <-
        wellData[, c("well", "agonist", "conc", "units")]

      print(
        paste(
          "processing: ",
          inputLoc,
          donorDirList[i],
          ",",
          "calcium",
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
              "calcium",
              "_",
              measure,
              ".png",
              sep = "")
      png(
        filename = OutF,
        width = 11,
        height = 8,
        units = 'in',
        res = 400
      )
      par(mfrow = c(3, 4))
      par(oma = c(0.2, 0, 1.5, 0.5)) #Sets outer region around plot (bottom, left, top, right)
      par(mgp = c(1.2, 0.3, 0)) #Sets space between data area and axis labels
      par(mar = c(5.4, 2, 1, 1.75)) #sets the space between plots (bottom, left, top, right)
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
        #print(paste("destfile: ", destfile))

        if (file.exists(destfile)) {
          df <-
            read.csv(destfile,
                     stringsAsFactors = FALSE,
                     header = TRUE)
          df.raw <- df[, c("FL1.A", "Time")]
          names(df.raw) <- c("fl", "timePoints")

          df.psych <-
            data.frame(describeBy(df.raw$fl, df.raw$timePoints, mat = TRUE))

          loe <-
            loess(df.raw$fl ~ df.raw$timePoints,
                  span = 0.1,
                  degree = 1)

          df.smooth <-
            data.frame(seq(from = 281, to = 2000), predict(loe, seq(from = 281, to = 2000)))
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
            data.frame(seq(from = 281, to = 2000), predict(loe, seq(from = 281, to = 2000)))
          names(df.smooth.extra) <- c("timePoints", "fl")

          maxPoint <- max(df.smooth.extra$fl, na.rm = TRUE)
          minPoint <- min(df.smooth.extra$fl, na.rm = TRUE)
          maxTime <- (which.max(df.smooth.extra$fl))/10 + 27

          #Trying to calculate rate of change of the smoothed data
          smoothPoints <- 200 #
          dY <-
            diff(df.smooth.extra$fl, lag = smoothPoints) # smoothing out the rate of change, basing diff on smoothPoints datapoints
          dX <- df.smooth.extra$timePoints + smoothPoints / 2
          chopdX <- length(dX) - smoothPoints #
          dX <- dX[1:chopdX]
          rateOfChange <-
            round(dY / smoothPoints, digits = 2) #so x distance of differential = smoothPoints

          maxRateOfChange <- max(rateOfChange)
          whereMaxRateOfChange <- which.max(rateOfChange)
          minRateOfChange <- min(rateOfChange)
          whereMinRateOfChange <- which.min(rateOfChange)

          dYY <-
            diff(df.smooth.extra$fl,
                 lag = smoothPoints,
                 differences = 2) # smoothing out the rate of change, basing diff on smoothPoints datapoints
          dXX <- dX + smoothPoints / 2
          chopdXX <- length(dXX) - smoothPoints #
          dXX <- dXX[1:chopdXX]
          rateOfAccel <-
            round(dYY / smoothPoints, digits = 2) #so x distance of differential = smoothPoints

          maxRateOfAccel <- max(rateOfAccel)
          whereMaxRateOfAccel <- which.max(rateOfAccel)
          minRateOfAccel <- min(rateOfAccel)
          whereMinRateOfAccel <- which.min(rateOfAccel)

          pred0 <- predict(loe, 281)
          pred10 <- predict(loe, 281 + 100)
          pred20 <- predict(loe, 281 + 200)
          pred30 <- predict(loe, 281 + 300)
          pred40 <- predict(loe, 281 + 400)
          pred60 <- predict(loe, 281 + 600)
          pred90 <- predict(loe, 281 + 900)
          pred100 <- predict(loe, 281 + 1000)
          pred120 <- predict(loe, 281 + 1200)
          pred180 <- predict(loe, 281 + 1800)

          if (is.na(pred120)) {
            pred120 <- 0
          }
          if (is.na(pred100)) {
            pred100 <- 0
          }
          if (is.na(pred10)) {
            pred120 <- 0
          }
          if (is.na(pred60)) {
            pred100 <- 0
          }
          if (is.na(pred60)) {
            pred90 <- 0
          }

          Drop <- round(((maxPoint - pred100)/(maxPoint- minPoint))*100)
          AbsoluteDrop <- round(maxPoint - pred100)

          RoC10 <- round((pred10 - pred0) / 10, 2)
          RoC20 <- round((pred20 - pred10) / 10, 2)
          RoC30 <- round((pred30 - pred20) / 10, 2)
          RoC40 <- round((pred40 - pred30) / 10, 2)
          RoC60 <- round((pred60 - pred40) / 20, 2)
          RoC90 <- round((pred90 - pred60) / 30, 2)
          RoC120 <- round((pred120 - pred90) / 30, 2)
          RoC180 <- round((pred180 - pred120) / 60, 2)

          ##Add a new label, early speed
          LabelRate <- "NA"
          LabelRateMetric <- round((pred30 - pred0) / 30, 2)
          if (is.na(LabelRateMetric)) {
            LabelRateMetric <- 0
          }
          if (LabelRateMetric > 80) {
               LabelRate <- "fast"
          }
          if (LabelRateMetric > 40 & LabelRateMetric <= 80) {
            LabelRate <- "medium"
          }
          if (LabelRateMetric <= 40) {
            LabelRate <- "slow"
          }
          if (LabelRateMetric <= 10) {
            LabelRate <- "no change"
          }

          Category <- "NA"
          CategoryMetric <- 0

          if (LabelRateMetric <= 40 &
              AbsoluteDrop < 1000) {
              Category <- "linear"
          } else
            {if (maxTime > 20 &
            maxTime < 100 &
            AbsoluteDrop > 2000) {
            Category <- "Peak"} else {
              Category <- "Sustained"
            }
            }

          #Now lets plot it all
          mainTitle = paste(wellDatas[k, "conc"], wellDatas[k, "units"], " ", wellDatas[k, "agonist"])

          if (measure == "FL1") {
            ylimit <- c(1, 800000)
          } else {
            ylimit <- c(1, 80000)
          }

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

          ylimit <- c(1, 20000)
          xlimit <- c(0, 140)

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
          abline(h = maxPoint) #These are new
          abline(v = maxTime, col = "blue")
          #text(120, 15000, LabelEarlyRateMetric)
          text(80, 12000, paste(LabelRate, "Rate of change"))
          #text(120, 10000, AbsoluteDrop)
          #text(120, 8000, Drop)
          #text(120, 6000, maxTime)
          text(80, 4000, paste("Category: ", Category))


          #Add a line to the summary file
          summaryFile[summaryFileCount, 1] <-
            donorDirList[i]
          summaryFile[summaryFileCount, 2] <- "calcium"
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
          summaryFile[summaryFileCount, 24] <- Drop
          summaryFile[summaryFileCount, 25] <-
            AbsoluteDrop
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
          summaryFile[summaryFileCount, 34] <- Category
          summaryFile[summaryFileCount, 35] <- CategoryMetric
          summaryFile[summaryFileCount, 36] <- LabelRate
          summaryFile[summaryFileCount, 37] <- LabelRateMetric
          summaryFileCount <- summaryFileCount + 1

        } else {
          print(paste(destfile, " does not exist"))
        }

      }
      dev.off()

    }

    #Save the summary
    summaryFileName <-
      paste(outputLoc, "summary", measure, "CalciumFile.csv", sep = "")

    write.csv(summaryFile, file = summaryFileName, row.names = FALSE)

    print("finished")

  }
