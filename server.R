library("shiny")

library("MALDIquant")

load("spectra.RData")

source("plotSlice.R")

s <- spectra

x <- unlist(lapply(s, function(x)metaData(x)$imaging$pos["x"]))
y <- unlist(lapply(s, function(x)metaData(x)$imaging$pos["y"]))

rx <- range(x)
ry <- range(y)

nr <- diff(ry)+1
nc <- diff(rx)+1

print(nr)
print(nc)

x <- x-(rx[1]-1)
y <- y-(ry[1]-1)

pixel <- matrix(FALSE, nrow=nr, ncol=nc, byrow=TRUE)
pixel[matrix(c(y, x), ncol=2)] <- TRUE

shinyServer(function(input, output, session) {

  output$selectROI <- renderUI({
    selectInput(inputId="selROI",
                label="Select Region Of Interest:",
                choices=roi(), selected=roi()[which.max(intensity(peaks()))],
                multiple=TRUE)
  })

  smoothedSpectra <- reactive({
    if (is.null(input$smHws)) {
      return(s)
    } else {
      return(suppressWarnings(smoothIntensity(s, method="SavitzkyGolay", halfWindowSize=input$smHws)))
    }
  })

  baselineCorrectedSpectra <- reactive({
    if (is.null(input$bcHws)) {
      return(smoothedSpectra())
    } else {
      return(removeBaseline(smoothedSpectra(), method="TopHat",
                            halfWindowSize=input$bcHws))
    }
  })

  meanSpectrum <- reactive({
    return(averageMassSpectra(baselineCorrectedSpectra()))
  })

  peaks <- reactive({
    return(detectPeaks(meanSpectrum(), method="MAD", SNR=input$SNR,
                       halfWindowSize=input$bcHws))
  })

  roi <- reactive({
    return(round(mass(peaks()), digits=3))
  })

  selRoi <- reactive({
    if (is.null(input$selROI)) {
      return(NULL)
    } else {
      i <- MALDIquant:::.which.closest(as.double(input$selROI), roi())
      return(roi()[i])
    }
  })

  roiRange <- reactive({
    if (is.null(input$roiTolerance)) {
      return(NULL)
    } else {
      return(lapply(selRoi(), function(x)c(x-input$roiTolerance,
                                           x+input$roiTolerance)))
    }
  })

  iMatrix <- reactive({
    return(.intensityMatrix(baselineCorrectedSpectra()))
  })

  ## taken from https://gist.github.com/wch/5436415
  output$plotSlices <- renderUI({
    plotOutputList <- lapply(seq_along(selRoi()), function(i) {
      plotname <- paste0("plot", i)
      plotOutput(plotname)
    })

    for (i in seq_along(selRoi())) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
      local({
        my_i <- i
        roiR <- roiRange()
        plotname <- paste0("plot", my_i)
        output[[plotname]] <- renderPlot({
          plotSlice(iMatrix(), roiR[[my_i]], pixel=pixel,
                    main=paste0("ROI: ", selRoi()[my_i], " +/- ",
                                input$roiTolerance, " m/z"))
        })
      })
    }

    return(plotOutputList)
  })

  output$plotMeanSpectrum <- renderPlot({
    plot(meanSpectrum(), main="Mean Spectrum")
    abline(v=selRoi(), col=2)
  })
})

