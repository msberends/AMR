#' AMR Extensions for Tidymodels
#'
#' This family of functions allows using AMR-specific data types such as `<sir>` and `<mic>` inside `tidymodels` pipelines.
#' @inheritParams recipes::step_center
#' @details
#' You can read more in our online [AMR with tidymodels introduction](https://amr-for-r.org/articles/AMR_with_tidymodels.html).
#'
#' Tidyselect helpers include:
#' - [all_sir()] and [all_sir_predictors()] to select [`<sir>`][as.sir()] columns
#' - [all_mic()] and [all_mic_predictors()] to select [`<mic>`][as.mic()] columns
#' - [all_disk()] and [all_disk_predictors()] to select [`<disk>`][as.disk()] columns
#'
#' Pre-processing pipeline steps include:
#' - [step_sir_numeric()] to convert SIR columns to numeric (via `as.numeric()`), to be used with [all_sir_predictors()]: `"S"` = 1, `"I"`/`"SDD"` = 2, `"R"` = 3. All other values are rendered `NA`. Keep this in mind for further processing, especially if the model does not allow for `NA` values.
#' - [step_mic_log2()] to convert MIC columns to numeric (via `as.numeric()`) and apply a log2 transform, to be used with [all_mic_predictors()]
#'
#' These steps integrate with `recipes::recipe()` and work like standard preprocessing steps. They are useful for preparing data for modelling, especially with classification models.
#' @seealso [recipes::recipe()], [as.sir()], [as.mic()], [as.disk()]
#' @name amr-tidymodels
#' @keywords internal
#' @export
#' @examples
#' if (require("tidymodels")) {
#'
#'   # The below approach formed the basis for this paper: DOI 10.3389/fmicb.2025.1582703
#'   # Presence of ESBL genes was predicted based on raw MIC values.
#'
#'
#'   # example data set in the AMR package
#'   esbl_isolates
#'
#'   # Prepare a binary outcome and convert to ordered factor
#'   data <- esbl_isolates %>%
#'     mutate(esbl = factor(esbl, levels = c(FALSE, TRUE), ordered = TRUE))
#'
#'   # Split into training and testing sets
#'   split <- initial_split(data)
#'   training_data <- training(split)
#'   testing_data <- testing(split)
#'
#'   # Create and prep a recipe with MIC log2 transformation
#'   mic_recipe <- recipe(esbl ~ ., data = training_data) %>%
#'
#'     # Optionally remove non-predictive variables
#'     remove_role(genus, old_role = "predictor") %>%
#'
#'     # Apply the log2 transformation to all MIC predictors
#'     step_mic_log2(all_mic_predictors()) %>%
#'
#'     # And apply the preparation steps
#'     prep()
#'
#'   # View prepped recipe
#'   mic_recipe
#'
#'   # Apply the recipe to training and testing data
#'   out_training <- bake(mic_recipe, new_data = NULL)
#'   out_testing <- bake(mic_recipe, new_data = testing_data)
#'
#'   # Fit a logistic regression model
#'   fitted <- logistic_reg(mode = "classification") %>%
#'     set_engine("glm") %>%
#'     fit(esbl ~ ., data = out_training)
#'
#'   # Generate predictions on the test set
#'   predictions <- predict(fitted, out_testing) %>%
#'     bind_cols(out_testing)
#'
#'   # Evaluate predictions using standard classification metrics
#'   our_metrics <- metric_set(accuracy,
#'                             recall,
#'                             precision,
#'                             sensitivity,
#'                             specificity,
#'                             ppv,
#'                             npv)
#'   metrics <- our_metrics(predictions, truth = esbl, estimate = .pred_class)
#'
#'   # Show performance
#'   metrics
#' }
all_sir <- function() {
  x <- tidymodels_amr_select(class = "sir")
  names(x)
}

#' @rdname amr-tidymodels
#' @export
all_sir_predictors <- function() {
  x <- tidymodels_amr_select(class = "sir")
  intersect(x, recipes::has_role("predictor"))
}

#' @rdname amr-tidymodels
#' @export
all_mic <- function() {
  x <- tidymodels_amr_select(class = "mic")
  names(x)
}

#' @rdname amr-tidymodels
#' @export
all_mic_predictors <- function() {
  x <- tidymodels_amr_select(class = "mic")
  intersect(x, recipes::has_role("predictor"))
}

#' @rdname amr-tidymodels
#' @export
all_disk <- function() {
  x <- tidymodels_amr_select(class = "disk")
  names(x)
}

#' @rdname amr-tidymodels
#' @export
all_disk_predictors <- function() {
  x <- tidymodels_amr_select(class = "disk")
  intersect(x, recipes::has_role("predictor"))
}

#' @rdname amr-tidymodels
#' @export
step_mic_log2 <- function(
    recipe,
    ...,
    role = NA,
    trained = FALSE,
    columns = NULL,
    skip = FALSE,
    id = recipes::rand_id("mic_log2")) {
  recipes::add_step(
    recipe,
    step_mic_log2_new(
      terms = rlang::enquos(...),
      role = role,
      trained = trained,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

step_mic_log2_new <- function(terms, role, trained, columns, skip, id) {
  recipes::step(
    subclass = "mic_log2",
    terms = terms,
    role = role,
    trained = trained,
    columns = columns,
    skip = skip,
    id = id
  )
}

#' @rawNamespace if(getRversion() >= "3.0.0") S3method(recipes::prep, step_mic_log2)
prep.step_mic_log2 <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], types = "ordered")
  step_mic_log2_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}

#' @rawNamespace if(getRversion() >= "3.0.0") S3method(recipes::bake, step_mic_log2)
bake.step_mic_log2 <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)
  for (col in object$columns) {
    new_data[[col]] <- log2(as.numeric(as.mic(new_data[[col]])))
  }
  new_data
}

#' @export
print.step_mic_log2 <- function(x, width = max(20, options()$width - 35), ...) {
  title <- "Log2 transformation of MIC columns"
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
}

#' @rawNamespace if(getRversion() >= "3.0.0") S3method(recipes::tidy, step_mic_log2)
tidy.step_mic_log2 <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble::tibble(terms = x$columns)
  } else {
    res <- tibble::tibble(terms = recipes::sel2char(x$terms))
  }
  res$id <- x$id
  res
}

#' @rdname amr-tidymodels
#' @export
step_sir_numeric <- function(
    recipe,
    ...,
    role = NA,
    trained = FALSE,
    columns = NULL,
    skip = FALSE,
    id = recipes::rand_id("sir_numeric")) {
  recipes::add_step(
    recipe,
    step_sir_numeric_new(
      terms = rlang::enquos(...),
      role = role,
      trained = trained,
      columns = columns,
      skip = skip,
      id = id
    )
  )
}

step_sir_numeric_new <- function(terms, role, trained, columns, skip, id) {
  recipes::step(
    subclass = "sir_numeric",
    terms = terms,
    role = role,
    trained = trained,
    columns = columns,
    skip = skip,
    id = id
  )
}

#' @rawNamespace if(getRversion() >= "3.0.0") S3method(recipes::prep, step_sir_numeric)
prep.step_sir_numeric <- function(x, training, info = NULL, ...) {
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names], types = "ordered")
  step_sir_numeric_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}

#' @rawNamespace if(getRversion() >= "3.0.0") S3method(recipes::bake, step_sir_numeric)
bake.step_sir_numeric <- function(object, new_data, ...) {
  recipes::check_new_data(object$columns, object, new_data)
  for (col in object$columns) {
    new_data[[col]] <- as.numeric(as.sir(new_data[[col]]))
  }
  new_data
}

#' @export
print.step_sir_numeric <- function(x, width = max(20, options()$width - 35), ...) {
  title <- "Numeric transformation of SIR columns"
  recipes::print_step(x$columns, x$terms, x$trained, title, width)
}

#' @rawNamespace if(getRversion() >= "3.0.0") S3method(recipes::tidy, step_sir_numeric)
tidy.step_sir_numeric <- function(x, ...) {
  if (recipes::is_trained(x)) {
    res <- tibble::tibble(terms = x$columns)
  } else {
    res <- tibble::tibble(terms = recipes::sel2char(x$terms))
  }
  res$id <- x$id
  res
}

tidymodels_amr_select <- function(class) {
  df <- get_current_data()
  ind <- which(
    vapply(
      FUN.VALUE = logical(1),
      df,
      function(x) inherits(x, class),
      USE.NAMES = TRUE
    ),
    useNames = TRUE
  )
  ind
}
