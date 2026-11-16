# AMR Extensions for Tidymodels

This family of functions allows using AMR-specific data types such as
`<sir>` and `<mic>` inside `tidymodels` pipelines.

## Usage

``` r
all_sir()

all_sir_predictors()

all_mic()

all_mic_predictors()

all_disk()

all_disk_predictors()

step_mic_log2(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  skip = FALSE,
  id = recipes::rand_id("mic_log2")
)

step_sir_numeric(
  recipe,
  ...,
  role = NA,
  trained = FALSE,
  columns = NULL,
  skip = FALSE,
  id = recipes::rand_id("sir_numeric")
)
```

## Arguments

- recipe:

  A recipe object. The step will be added to the sequence of operations
  for this recipe.

- ...:

  One or more selector functions to choose variables for this step. See
  [`selections()`](https://recipes.tidymodels.org/reference/selections.html)
  for more details.

- role:

  Not used by this step since no new variables are created.

- trained:

  A logical to indicate if the quantities for preprocessing have been
  estimated.

- skip:

  A logical. Should the step be skipped when the recipe is baked by
  [`bake()`](https://recipes.tidymodels.org/reference/bake.html)? While
  all operations are baked when
  [`prep()`](https://recipes.tidymodels.org/reference/prep.html) is run,
  some operations may not be able to be conducted on new data (e.g.
  processing the outcome variable(s)). Care should be taken when using
  `skip = TRUE` as it may affect the computations for subsequent
  operations.

- id:

  A character string that is unique to this step to identify it.

## Details

You can read more in our online [AMR with tidymodels
introduction](https://amr-for-r.org/articles/AMR_with_tidymodels.html).

Tidyselect helpers include:

- `all_sir()` and `all_sir_predictors()` to select
  [`<sir>`](https://amr-for-r.org/release/3.1.0/reference/as.sir.md)
  columns

- `all_mic()` and `all_mic_predictors()` to select
  [`<mic>`](https://amr-for-r.org/release/3.1.0/reference/as.mic.md)
  columns

- `all_disk()` and `all_disk_predictors()` to select
  [`<disk>`](https://amr-for-r.org/release/3.1.0/reference/as.disk.md)
  columns

Pre-processing pipeline steps include:

- `step_sir_numeric()` to convert SIR columns to numeric (via
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html)), to be used
  with `all_sir_predictors()`: `"S"` = 1, `"I"`/`"SDD"` = 2, `"R"` = 3.
  All other values are rendered `NA`. Keep this in mind for further
  processing, especially if the model does not allow for `NA` values.

- `step_mic_log2()` to convert MIC columns to numeric (via
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html)) and apply a
  log2 transform, to be used with `all_mic_predictors()`

These steps integrate with
[`recipes::recipe()`](https://recipes.tidymodels.org/reference/recipe.html)
and work like standard preprocessing steps. They are useful for
preparing data for modelling, especially with classification models.

## See also

[`recipes::recipe()`](https://recipes.tidymodels.org/reference/recipe.html),
[`as.sir()`](https://amr-for-r.org/release/3.1.0/reference/as.sir.md),
[`as.mic()`](https://amr-for-r.org/release/3.1.0/reference/as.mic.md),
[`as.disk()`](https://amr-for-r.org/release/3.1.0/reference/as.disk.md)

## Examples

``` r
if (require("tidymodels")) {
  # The below approach formed the basis for this paper: DOI 10.3389/fmicb.2025.1582703
  # Presence of ESBL genes was predicted based on raw MIC values.


  # example data set in the AMR package
  esbl_isolates

  # Prepare a binary outcome and convert to ordered factor
  data <- esbl_isolates %>%
    mutate(esbl = factor(esbl, levels = c(FALSE, TRUE), ordered = TRUE))

  # Split into training and testing sets
  split <- initial_split(data)
  training_data <- training(split)
  testing_data <- testing(split)

  # Create and prep a recipe with MIC log2 transformation
  mic_recipe <- recipe(esbl ~ ., data = training_data) %>%
    # Optionally remove non-predictive variables
    remove_role(genus, old_role = "predictor") %>%
    # Apply the log2 transformation to all MIC predictors
    step_mic_log2(all_mic_predictors()) %>%
    # And apply the preparation steps
    prep()

  # View prepped recipe
  mic_recipe

  # Apply the recipe to training and testing data
  out_training <- bake(mic_recipe, new_data = NULL)
  out_testing <- bake(mic_recipe, new_data = testing_data)

  # Fit a logistic regression model
  fitted <- logistic_reg(mode = "classification") %>%
    set_engine("glm") %>%
    fit(esbl ~ ., data = out_training)

  # Generate predictions on the test set
  predictions <- predict(fitted, out_testing) %>%
    bind_cols(out_testing)

  # Evaluate predictions using standard classification metrics
  our_metrics <- metric_set(
    accuracy,
    recall,
    precision,
    sensitivity,
    specificity,
    ppv,
    npv
  )
  metrics <- our_metrics(predictions, truth = esbl, estimate = .pred_class)

  # Show performance
  metrics
}
```
