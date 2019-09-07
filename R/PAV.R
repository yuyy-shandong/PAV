
#' Eliminating the Unmeasured Confounders and Estimating Causal Effect.
#'
#' @param formula an object of class "formula" : a symbolic description of the model to be fitted
#' @param data an optional data frame containing the variables in the model.
#' @param sdmethod Choosing a method for estimating variance
#' @param x1.name the name of exposure
#' @param x2.name the name of post-exposure
#' @param boots.no the number of bootstrap
#'
#'
#' @return  \code{formula} the model
#' @return  \code{bootstrap} the casual effect coefficients when the variance is calculated by bootstrap
#' @return  \code{normal} the casual effect coefficients when the variance is calculated in theory
#' @return  \code{Quant} the confidence interval calculated by bootstrap
#' @examples
#' PAV <-  function(y~x1+x2, data = data, sdmethod ="normal",x1.name = "x1",x2.name = "x2",boots.no = NULL)
#'
#'
#'
#'

PAV <-  function(formula,
                 data = data,
                 sdmethod = c("normal", "bootstrap", "all"),
                 x1.name = "x1",
                 x2.name = "x2",
                 boots.no = NULL) {

  data <- na.omit(data)
  model_PAV <- lm(formula, data = data)
  coef1_PAV <- summary(model_PAV)$coefficients[2, 1]
  coef2_PAV <- summary(model_PAV)$coefficients[3, 1]
  result1_PAV <-  coef1_PAV - coef2_PAV


  if (sdmethod == "bootstrap") {
    result_PAV.b <- NULL
    for (k in 1:boots.no) {
      sample_boot <- data[sample(1:dim(data)[1],
                                 dim(data), replace = T), ]
      model_PAV_b <- lm (formula = formula,
                         data = sample_boot)
      PAV_b_1 <- summary (model_PAV_b)$coefficients[2, 1]
      PAV_b_2 <- summary (model_PAV_b)$coefficients[3, 1]

      result_PAV.b[k] <- PAV_b_1 - PAV_b_2

    }

    sd_PAV.p <- sd(result_PAV.b)
    p1.b <-  pnorm(result1_PAV,
                   mean = 0,
                   sd = sd_PAV.p,
                   lower.tail = TRUE, log.p = FALSE)
    p2.b <- 1 - pnorm(result1_PAV, mean = 0,
                      sd = sd_PAV.p,
                      lower.tail = TRUE,
                      log.p = FALSE)
    p.b <- 2 * min(p1.b, p2.b)
    conup.b  <- qnorm( 0.975,
                       mean = result1_PAV,
                       sd = sd_PAV.p,
                       lower.tail = TRUE,
                       log.p = FALSE)
    conlower.b <- qnorm( 0.025,
                         mean = result1_PAV,
                         sd = sd_PAV.p,
                         lower.tail = TRUE,
                         log.p = FALSE )
    conup.b_q <- quantile(result_PAV.b, c(0.025, 0.975))[2]
    conlower.b_q <- quantile(result_PAV.b, c(0.025, 0.975))[1]
    result_boot <- list()
    result_boot$formula <- formula
    boos_coef <- cbind(result1_PAV, sd_PAV.p, p.b, conlower.b, conup.b)
    colnames(boos_coef) <-  c("Estimation", "sd",
                              "P_value", "2.5%", "97.5%")
    result_boot$bootstrap <- boos_coef
    quant <- cbind(conlower.b_q, conup.b_q)
    colnames(quant) <- c("2.5%", "97.5%")
    result_boot$Quant <- quant



    return(result_boot)

  } else{
    if (sdmethod == "normal") {
      covqian <- cbind(data[, x1.name], data[, x2.name])
      cov <- cbind(1, covqian)
      r <- qr(cov)$rank
      n <- nrow(cov)
      resivar <- sum(summary(model_PAV)$residuals ^ 2) / (n - r)
      cov2 <- solve(t(cov) %*% cov)
      c <- cbind(0, 1, -1)
      varco <- c %*% cov2 %*% t(c)
      sd_PAV <- (resivar * varco) ^ (1 / 2)
      z_PAV <- result1_PAV / sd_PAV

      p1_PAV <- pnorm(
        z_PAV,
        mean = 0,
        sd = 1,
        lower.tail = TRUE,
        log.p = FALSE)
      p2_PAV <- 1 - pnorm(
        z_PAV,
        mean = 0,
        sd = 1,
        lower.tail = TRUE,
        log.p = FALSE)

      p_PAV <- 2 * min(p1_PAV, p2_PAV)

      conup_PAV <- qnorm(0.975, result1_PAV, sd_PAV)
      conlower_PAV <- qnorm(0.025, result1_PAV, sd_PAV)


      result_normal <- list()
      result_normal$formula <- formula
      coef <- cbind(result1_PAV, sd_PAV, z_PAV,
                    p_PAV, conlower_PAV, conup_PAV)
      colnames(coef) <- c("Estimation", "sd", "Z",
                          "P_value", "2.5%", "97.5%")
      result_normal$normal <- coef

      return(result_normal)

    } else{
      if (sdmethod == "all"){
        result_PAV.b <- NULL
        for (k in 1:boots.no) {
          sample_boot <- data[sample(1:dim(data)[1],
                                     dim(data), replace = T), ]
          model_PAV_b <- lm (formula,
                             data = sample_boot)
          PAV_b_1 <- summary (model_PAV_b)$coefficients[2, 1]
          PAV_b_2 <- summary (model_PAV_b)$coefficients[3, 1]

          result_PAV.b[k] <- PAV_b_1 - PAV_b_2

        }

        sd_PAV.p <- sd(result_PAV.b)
        p1.b <-  pnorm(result1_PAV,
                       mean = 0,
                       sd = sd_PAV.p,
                       lower.tail = TRUE, log.p = FALSE)
        p2.b <- 1 - pnorm(result1_PAV, mean = 0,
                          sd = sd_PAV.p,
                          lower.tail = TRUE,
                          log.p = FALSE)
        p.b <- 2 * min(p1.b, p2.b)
        conup.b  <- qnorm( 0.975,
                           mean = result1_PAV,
                           sd = sd_PAV.p,
                           lower.tail = TRUE,
                           log.p = FALSE)
        conlower.b <- qnorm( 0.025,
                             mean = result1_PAV,
                             sd = sd_PAV.p,
                             lower.tail = TRUE,
                             log.p = FALSE )
        conup.b_q <- quantile(result_PAV.b, c(0.025, 0.975))[2]
        conlower.b_q <- quantile(result_PAV.b, c(0.025, 0.975))[1]
        result_boot <- list()
        result_boot$formula <- formula
        boos_coef <- cbind(result1_PAV, sd_PAV.p, p.b, conlower.b, conup.b)
        colnames(boos_coef) <-  c("Estimation", "sd",
                                  "P_value", "2.5%", "97.5%")
        result_boot$Coefficients <- boos_coef
        quant <- cbind(conlower.b_q, conup.b_q)
        colnames(quant) <- c("2.5%", "97.5%")
        result_boot$Quant <- quant


        covqian <- cbind(data[, x1.name], data[, x2.name])
        cov <- cbind(1, covqian)
        r <- qr(cov)$rank
        n <- nrow(cov)
        resivar <- sum(summary(model_PAV)$residuals ^ 2) / (n - r)
        cov2 <- solve(t(cov) %*% cov)
        c <- cbind(0, 1, -1)
        varco <- c %*% cov2 %*% t(c)
        sd_PAV <- (resivar * varco) ^ (1 / 2)
        z_PAV <- result1_PAV / sd_PAV

        p1_PAV <- pnorm(
          z_PAV,
          mean = 0,
          sd = 1,
          lower.tail = TRUE,
          log.p = FALSE)
        p2_PAV <- 1 - pnorm(
          z_PAV,
          mean = 0,
          sd = 1,
          lower.tail = TRUE,
          log.p = FALSE)

        p_PAV <- 2 * min(p1_PAV, p2_PAV)

        conup_PAV <- qnorm(0.975, result1_PAV, sd_PAV)
        conlower_PAV <- qnorm(0.025, result1_PAV, sd_PAV)


        result_normal <- list()
        result_normal$formula <- formula
        coef <- cbind(result1_PAV, sd_PAV, z_PAV, p_PAV,
                      conlower_PAV, conup_PAV)
        colnames(coef) <- c("Estimation", "sd", "Z",
                            "P_value", "2.5%", "97.5%")
        result_normal$Coefficients <- coef

        result_all <- list()
        result_all$formula <- formula
        result_all$bootstrap <- boos_coef
        result_all$normal <- coef
        result_all$Quant <- quant


        return(result_all)
      }

    }

  }








}
