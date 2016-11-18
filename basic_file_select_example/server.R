
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyFiles)

shinyServer(function(input, output, session) {

  output$distPlot <- renderPlot({

    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')

  })
  volumes <- c("C:"="C:\\")
  shinyDirChoose(input, 'directory', roots=volumes, session=session, restrictions=system.file(package='base'))

    # make dynamic file list
  filelist = reactive({ list.files(parseDirPath(volumes, input$directory), pattern = "*.cnv") })

    # update select input when filelist changes
  observe({ updateSelectInput(session, "files", choices = filelist() )})
  output$rawfilelist <- renderPrint(filelist())
})
