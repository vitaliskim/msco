#' A generalised B-spline model evaluation using Cross-validation approaches
#'
#' This function implements four different cross-validation techniques to evaluate the
#'  predictive ability of the generalised B-spline model (*sensu* Lagat *et al.,* 2021a).
#'   The four different techniques implemented are:
#'  * validation set approach;
#'  * k -fold;
#'  * Leave-one-out-cross-validation (LOOCV), and
#'  * Repeated k-fold.
#'
#' The k-fold cross-validation approach is highly recommended due to its computational
#'  efficiency and an acceptable bias-variance trade-off, subject to the value of `k`
#'  chosen to be either 5 or 10 (Lagat *et al.,* 2021a). For more details on the other
#'   cross-validation approaches, see Lagat *et al.* (2021b).
#'
#' @param gbsm_obj An object of `class` `"gbsm"` (i.e., assigned to \link[msco]{gbsm} function).
#' @param type The type of the cross-validation approach used. It must be
#'  \eqn{\in \{}"validation.set", "k-fold", "LOOCV", "repeated.k-fold" \eqn{\}}
#' @param p The percentage of data used in training the model. The value is used if the
#' cross-validation approach implemented is "validation.set".
#' @param k The value of `k` used in both "k-fold" and "repeated.k-fold" types of
#'  cross-validation. This value represents the number of subsets or groups that a
#'   given sample of data is to be split into. A value of 5 or 10 is used in practice,
#'    as it leads to an ideal bias-variance trade-off (Lagat *et al.,* 2021a).
#' @param k_fold.repeats The number of replicates used in "repeated.k-fold" type of
#'  cross-validation.
#' @return Depending on the type of cross-validation approach implemented, the `cross_valid`
#'  function returns:
#' * a `data.frame` with the following test errors (for "validation.set"):
#'   + `RMSE`: &nbsp; A root mean squared error;
#'   + `R_squared`: &nbsp; the Pearson's \eqn{r^2}, and
#'   + `MAE`: &nbsp; the mean absolute error.
#' * an `array` with test errors as above including the type of the
#'  regression model used, size of the samples, number of predictors,
#'   type of cross-validation performed, and summary of sample sizes
#'    (for "k-fold", "LOOCV", and "repeated.k-fold").
#' @references
#' \enumerate{
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021a). *Dissecting the effects of
#'   neutral encounter versus functional traits on multi-order species interactions
#'    and co-occurrence with generalised B-spline modelling*. Submitted.}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). *`msco`: an R software package
#'   for null model testing of multi-species interactions and interference with
#'    covariates*. Submitted.}
#'  }
#'
#' @examples
#' \dontrun{
#'
#' my.path <- system.file("extdata/gsmdat", package = "msco")
#' setwd(my.path)
#' s.data <- get(load("s.data.csv")) ## Species-by-site matrix
#' t.data <- get(load("t.data.csv")) ## Species-by-trait matrix
#' p.d.mat <- get(load("p.d.mat.csv")) ## Species-by-species phylogenetic distance matrix
#'
#' gbsm_obj <- msco::gbsm(s.data, t.data, p.d.mat, metric= "Simpson_eqn", d.f=4,
#'  order.jo=3, degree=3, n=1000, b.plots=FALSE, bsplines="single", scat.plot=FALSE,
#'   response.curves=FALSE, leg=1, start=seq(-0.1, 0, length.out=(ncol(t.data)+2)*4+1))
#'
#' val.set <- msco::cross_valid(gbsm_obj, type="validation.set", p=0.8)
#' val.set
#'
#' kfold <- msco::cross_valid(gbsm_obj, type="k-fold", k=5)
#' kfold
#'
#' loocv <- msco::cross_valid(gbsm_obj, type="LOOCV")
#' loocv
#'
#' repeated.kfold <- msco::cross_valid(gbsm_obj, type="repeated.k-fold", k=5, k_fold.repeats=100)
#' repeated.kfold
#'
#' }
#' @export
#' @md

cross_valid <- function(gbsm_obj, type="k-fold", p, k, k_fold.repeats){

  data <- gbsm_obj$bs_pred
  newdat <- `names<-`(data.frame(gbsm_obj$j.occs, data), c("j.occs", names(data)))

  cvad <- data.frame()
  if(type=="validation.set"){

    #split the dataset into a training set ((p*100)%) and test set ((1-(p*100))%).
    `%>%` <- dplyr::`%>%`
    training_obs <- newdat$j.occs %>% caret::createDataPartition(p = p, list = FALSE)

    train <- newdat[training_obs, ]
    test <- newdat[-training_obs, ]

    # Build the generalized linear regression model on the training set
    model <- suppressWarnings(glm2::glm2(j.occs ~ ., family=stats::binomial(link="log"), data = train, start = gbsm_obj$start))

    # Use the model to make predictions on the test set
    predictions <- suppressWarnings(stats::predict.glm(model, newdata = test, type = "response"))

    #Examine R-squared, RMSE, and MAE of predictions
    cvad <- data.frame(RMSE = caret::RMSE(predictions, test$j.occs),
               R_squared = caret::R2(predictions, test$j.occs),
               MAE = caret::MAE(predictions, test$j.occs))

  }else if(type=="k-fold"){

    #define the number of subsets (or "folds") to use
    train_control <- caret::trainControl(method = "cv", number = k)

    #train the model
    model <- suppressWarnings(caret::train(j.occs ~ ., data = newdat, method = "glm", trControl = train_control))

    #Summarize the results
    cvad <- model

  }else if(type=="LOOCV"){

    #specify that we want to use LOOCV
    train_control <- caret::trainControl(method = "LOOCV")

    #train the model
    model <- suppressWarnings(caret::train(j.occs ~ ., data = newdat, method = "glm", trControl = train_control))

    #summarize the results
    cvad <- model

  }else if(type=="repeated.k-fold"){

    #define the number of subsets to use and number of times to repeat k-fold CV
    train_control <- caret::trainControl(method = "repeatedcv", number = k, repeats = k_fold.repeats)

    #train the model
    model <- suppressWarnings(caret::train(j.occs ~ ., data = newdat, method = "glm", trControl = train_control))

    #summarize the results
    cvad <- model
  }else{
    stop("Wrong type of cross-validation! The only available types are 'validation.set', 'k-fold', 'LOOCV', or 'repeated.k-fold'.")
  }
  return(cvad)
}

