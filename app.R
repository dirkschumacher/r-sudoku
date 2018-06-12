# Authored by Dirk Schumacher
# License: GPLv3
library(ROI)
library(ROI.plugin.glpk)
library(ompr) # devtools::install_github("dirkschumacher/ompr")
library(ompr.roi) # devtools::install_github("dirkschumacher/ompr.roi")
library(dplyr)
library(purrr)
library(purrrlyr)
library(lazyeval)

# the base model
n <- 9
if (!file.exists("base_model.rds")) {
  base_model <- MIPModel() %>%

    # The number k stored in position i,j
    add_variable(
      x[i, j, k],
      i = 1:n,
      j = 1:n,
      k = 1:9,
      type = "binary"
    ) %>%

    # no objective
    set_objective(0) %>%

    # only one number can be assigned per cell
    add_constraint(sum_expr(x[i, j, k], k = 1:9) == 1, i = 1:n, j = 1:n) %>%

    # each number is exactly once in a row
    add_constraint(sum_expr(x[i, j, k], j = 1:n) == 1, i = 1:n, k = 1:9) %>%

    # each number is exactly once in a column
    add_constraint(sum_expr(x[i, j, k], i = 1:n) == 1, j = 1:n, k = 1:9) %>%

    # each 3x3 square must have all numbers
    add_constraint(
      sum_expr(x[i, j, k], i = 1:3 + sx, j = 1:3 + sy) == 1,
      sx = seq(0, n - 3, 3),
      sy = seq(0, n - 3, 3),
      k = 1:9
    )

  saveRDS(base_model, file = "base_model.rds")
} else {
  base_model <- readRDS("base_model.rds")
}

# helper function to create dynamic sudoku tables
create_sudoku_table <- function(code) {
  tags$table(map(seq_len(n), function(i) {
    tags$tr(map(seq_len(n), function(j) {
      tags$td(code(i, j))
    }))
  }))
}

server <- function(input, output) {

  # this constructs the current model
  # and adds constraints dynamically
  current_model <- eventReactive(input$solve, {
    combinations <- expand.grid(i = seq_len(n),
                                j = seq_len(n))
    constraints <-
      map2(combinations$i, combinations$j, function(i, j) {
        input_name <- paste0("n", i, j, collapse = "_")
        value <- as.numeric(isolate(input[[input_name]]))
        if (is.numeric(value) && !is.na(value)) {
          interp(lazy(x[i, j, v] == 1),
                 i = i,
                 j = j,
                 v = value)
        } else {
          NULL
        }
      }) %>% discard(is.null)
    reduce(constraints, function(acc, constr) {
      add_constraint_(acc, constr)
    }, .init = base_model)
  })

  # the current solution
  current_solution <- reactive({
    if (is.null(current_model())) {
      model <- base_model
    } else {
      model <- current_model()
    }
    withProgress(solve_model(current_model(), with_ROI("glpk")),
                 message = "Finding solution...",
                 value = 1)
  })

  # register output grid
  in_memory <- reactiveValues()
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      o_key <- paste0("o", i, j, collapse = "_")
      output[[o_key]] <-
        f_eval(interp(~ renderText({
          in_memory[[o_key]]
        })), list(o_key = o_key))
    }
  }

  solution_status <- reactive({
    solver_status(current_solution())
  })

  # update the right table whenever
  # GLPK found a valid solution
  observe({
    solution <- current_solution()
    status <- solution_status()
    if (status == "optimal") {
      solution %>%
        get_solution(x[i, j, k]) %>%
        filter(value > 0.9) %>%
        select(i, j, k) %>%
        by_row(function(row) {
          i <- row$i
          j <- row$j
          k <- row$k
          o_key <- paste0("o", i, j, collapse = "_")
          in_memory[[o_key]] <- k
        })
    }
  })

  output$solver_status <- renderText({
    if (solution_status() == "optimal") {
      "Found solution"
    } else {
      "Could not find a solution. Displaying the old result."
    }
  })

  output$sudoku_constraints <- renderUI({
    create_sudoku_table(function(i, j) {
      numericInput(
        paste0("n", i, j, collapse = "_"),
        "",
        if (i == j)
          i
        else
          NULL,
        min = 1,
        max = 9,
        step = 1,
        width = 50
      )
    })
  })

  output$sudoku_solution <- renderUI({
    create_sudoku_table(function(i, j) {
      textOutput(paste0("o", i, j, collapse = "_"), inline = TRUE)
    })
  })
}

ui <- fluidPage(
  includeHTML("github.html"),
  tags$style(
    # table css by GCyrillus (https://stackoverflow.com/users/2442099/gcyrillus)
    # http://stackoverflow.com/a/23497700
    "table {
    margin:1em auto;
    }
    td {
    height:30px;
    width:40px;
    border:1px solid;
    text-align:center;
    }
    td:first-child {
    border-left:solid;
    }
    td:nth-child(3n) {
    border-right:solid ;
    }
    tr:first-child {
    border-top:solid;
    }
    tr:nth-child(3n) td {
    border-bottom:solid ;
    }
    #sudoku_constraints .form-group {margin: 0}
    #sudoku_constraints label {display: none;}"
  ),
  h1("Solve Sudokus with GLPK and R"),
  shiny::fixedRow(
    column(
      5,
      h2("Enter your constraints"),
      uiOutput("sudoku_constraints")
    ),
    column(
      2,
      actionButton("solve", "Solve =>", class = "btn-success"),
      textOutput("solver_status")
    ),
    column(4,
           h2("Solution"),
           uiOutput("sudoku_solution")),
    height = "100px"
  ),
  hr(),
  p(
    "Solved with the R package ",
    a("ompr", href = "https://github.com/dirkschumacher/ompr")
  ),
  p(
    "Code is on ",
    a("Github", href = "https://github.com/dirkschumacher/r-sudoku")
  )
)

shinyApp(ui = ui, server = server)
