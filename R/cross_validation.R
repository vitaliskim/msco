#' Cross validation of the generalised B-spline model
#'
#' This function implements four different cross-validation techniques to evaluate the
#'  predictive ability of the generalised B-spline model (*sensu* Lagat *et al.,* 2021b).
#'   The four different techniques implemented are:
#'  * validation set approach;
#'  * k -fold;
#'  * Leave-one-out-cross-validation (LOOCV), and
#'  * Repeated k-fold.
#'
#' The k-fold cross-validation approach is highly recommended due to its computational
#'  efficiency and an acceptable bias-variance trade-off, subject to the value of `k`
#'  chosen to be either 5 or 10 (Lagat *et al.,* 2021b). For more details on the other
#'   cross-validation approaches, see Lagat *et al.* (2021c).
#'
#' @param gbsm_obj An object of `class` `"gbsm"` (i.e., assigned to \link[msco]{gbsm} function).
#' @param type The type of the cross-validation approach used. It must be
#'  \eqn{\in \{}"validation.set", "k-fold", "LOOCV", "repeated.k-fold" \eqn{\}}
#' @param p The percentage (in decimal form) of data used in training the model. The value is used if the
#' cross-validation approach implemented is "validation.set".
#' @param k The value of `k` used in both "k-fold" and "repeated.k-fold" types of
#'  cross-validation. This value represents the number of subsets or groups that a
#'   given sample of data is to be split into. A value of 5 or 10 is used in practice,
#'    as it leads to an ideal bias-variance trade-off (Lagat *et al.,* 2021b).
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
#'  \item{Fushiki, T. (2011). Estimation of prediction error by using K-fold cross-validation.
#'   *Stat. Comput.* **21**, 137-146. <https://doi.org/10.1007/s11222-009-9153-8>}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021b). *Dissecting the effects of random
#'   encounter versus functional trait mismatching on multi-species co-occurrence and
#'    interference with generalised B-spline modelling*. DOI: `<To be added>`.}
#'
#'  \item{Lagat, V. K., Latombe, G. and Hui, C. (2021c). *`msco`: an R software package
#'   for null model testing of multi-species interactions and interference with
#'    covariates*. DOI: `<To be added>`.}
#'
#'  \item{Pearson, K. (1895) VII. Note on regression and inheritance in the
#'  case of two parents. *proceedings of the royal society of London,* **58**:240-242.
#'   <https://doi.org/10.1098/rspl.1895.0041>}
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
#'   response.curves=FALSE, leg=1, max.vif, max.vif2, start.range=c(-0.1,0))
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

  if(sum(newdat$j.occs)==0){
    stop("The model cannot be trained since all `j.occs` are zero")
  }else if(sum(newdat$j.occs)>0){
    cvad <- data.frame()
    if(type=="validation.set"){

      #split the dataset into a training set ((p*100)%) and test set ((1-(p*100))%).
      `%>%` <- dplyr::`%>%`
      training_obs <- newdat$j.occs %>% caret::createDataPartition(p = p, list = FALSE)

      train <- newdat[training_obs, ]
      test <- newdat[-training_obs, ]

      # Build the generalized linear regression model on the training set
      model <- suppressWarnings(glm2::glm2(j.occs ~ ., family=stats::binomial(link="log"), data = train,
                                           start = seq(gbsm_obj$start.range[1], gbsm_obj$start.range[2], length.out=(ncol(data))+1)))

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
}

