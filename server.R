library("shiny")

library("MALDIquant")

load("spectra.RData")

s <- spectra

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
        plotname <- paste0("plot", my_i)
        output[[plotname]] <- renderPlot({
          par(mar=c(0L, 0L, 2L, 0L))
          plotMsiSlice(baselineCorrectedSpectra(), center=selRoi()[my_i],
                       tolerance=input$roiTolerance,
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

