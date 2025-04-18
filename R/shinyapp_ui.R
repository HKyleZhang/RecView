#' Function building dashboard UI, used in Shiny app
#'
#' @noRd
recview_ui <- function(request) {
  fluidPage(
    shinyjs::useShinyjs(),
    h6("version 1.1.0"),
    imageOutput("title", height = "70px", width = "383px"),
    h5("Maintainer: Hongkai Zhang"),
    sidebarPanel(width = 3,
      shinyFiles::shinyFilesButton(id = "genofile", label = "Select genotype file", title = "Please select a file", multiple = FALSE, viewtype = "detail", buttonType = "primary btn-block"),
      conditionalPanel(
        condition = "output.genoexist == 1",
        uiOutput(outputId = "genofile_attached"),
      ),
      p(),
      shinyFiles::shinyFilesButton(id = "scfile", label = "Select scaffold file", title = "Please select a file", multiple = FALSE, viewtype = "detail", buttonType = "primary btn-block"),
      conditionalPanel(
        condition = "output.scexist == 1",
        uiOutput(outputId = "scfile_attached"),
      ),
      conditionalPanel(
        condition = "output.genoexist == 1",
        uiOutput(outputId = "off"),
      ),
      conditionalPanel(
        condition = "output.scexist == 1",
        uiOutput(outputId = "chr"),
      ),
      hr(style="border-color: black;"),
      checkboxGroupInput(inputId = "resolution", label = "Figure resolution:", choices =  "Low", selected = "Low"),
      radioButtons(inputId = "locate", label = "Locate recombination positions?", choices =  c("Yes", "No"), selected = "No"),
      conditionalPanel(
        condition = "input.locate != 'Yes'",
        checkboxGroupInput(inputId = "save_opt", label = "Saving options:", choices = c("GoO Inferences", "Plots"))
      ),
      conditionalPanel(
        condition = "input.locate == 'Yes'",
        radioButtons(inputId = "algorithm", label = "Algorithm:", choices = c("PD", "CCS"), selected = "CCS"),
        conditionalPanel(
          condition = "input.algorithm == 'PD'",
          textInput(inputId = "radius", label = "Radius", width = "100%", placeholder = "default: 550")
        ),
        conditionalPanel(
          condition = "input.algorithm == 'PD'",
          textInput(inputId = "step", label = "Step", width = "100%", placeholder = "default: 17")
        ),
        conditionalPanel(
          condition = "input.algorithm == 'PD'",
          textInput(inputId = "finer_step", label = "Finer step", width = "100%", placeholder = "default: 1")
        ),
        conditionalPanel(
          condition = "input.algorithm == 'PD'",
          textInput(inputId = "pd_threshold", label = "Threshold to initiate finer step", width = "100%", placeholder = "default: 0.9")
        ),
        conditionalPanel(
          condition = "input.algorithm == 'CCS'",
          textInput(inputId = "threshold", label = "Threshold value", width = "100%", placeholder = "default: 50")
        ),
        checkboxGroupInput(inputId = "save_opt2", label = "Saving options", choices = c("GoO Inferences", "Plots", "Locations"))
      ),
      actionButton(inputId = "click", "Run analysis")
    ),
    h6(),
    mainPanel(width = 9,
      conditionalPanel(
        condition = "input.click > 0",
        uiOutput(outputId = "show_chr")),
      
      conditionalPanel(
        condition = "input.click > 0",
        uiOutput(outputId = "ui_out") %>% shinycssloaders::withSpinner(color="#008AA4", type = 6, size = 1.8)
      )
    )
  )
}
