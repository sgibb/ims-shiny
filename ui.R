library("shiny")

shinyUI(pageWithSidebar(

  headerPanel("MALDIquant - MSI example"),

  sidebarPanel(
    sliderInput(inputId="smHws", label="Savitzky-Golay Filter - halfWindowSize:",
                min=1, max=125, value=8),
    sliderInput(inputId="bcHws", label="TopHat Baseline - halfWindowSize:",
                min=1, max=250, value=16),
    sliderInput(inputId="SNR", label="Peak Detection - SNR:",
                min=1, max=100, value=15),
    sliderInput(inputId="roiTolerance",
                label="Region Of Interest - Range (in +/- m/z):",
                min=0, max=10, step=0.25, value=0.5),
    uiOutput("selectROI")
  ),

  mainPanel(
    plotOutput("plotMeanSpectrum"),
    uiOutput("plotSlices")
  )
))
