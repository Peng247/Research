
### Beta smooth only function for quantitative outcome
#' Generalized linear function regression model
#' @description
#' To fit generalized linear function regression model
#'
#' @param Y a numeric vector of response, which can be continuous or binary
#' @param X a numeric matrix of mixture exposure
#' @param Z a data frame of covariates, default value is NULL for no covariates model
#' @param q a number to decide quantile setting
#' @param sort.mx logical value to decide if smoothing function for exposure is needed
#' @param clustering logical value to decide if hierarchical clustering is needed for exposure matrix to sort variables in the way that correlated ones are closer
#' @param n.order.beta order of expansion model for mixture effect function
#' @param n.basis.beta number of basis functions for mixture effect function
#' @param base.system.beta basis function system for mixture effect function
#' @param n.order.mxex order of expansion model for mixture exposure function
#' @param n.basis.mxex number of basis functions for mixture exposure function
#' @param base.system.mxex basis function system for mixture exposure function
#' @param family specify link function, options are "gaussian" and "binomial"
#' @param screening logical value to specify outputs, if TRUE, mixture effect function is estimated
#'
#' @return Not Applicable
#' @export
#'
#' @examples Y=y;X=x;Z=z;q=4;sort.mx="ASC";mxex.smooth = FALSE;n.order.beta=4;n.basis.beta=7;base.system.beta="B-spline";n.order.mxex=4;n.basis.mxex=7;base.system.mxex="B-spline";family="binomial"
glfm <- function(Y, X, Z = NULL, q = 4, mxex.smooth = FALSE, sort.mx = "CR-direct",
                 n.order.beta = 4, n.basis.beta = 7, base.system.beta = "B-spline",
                 n.order.mxex = 4, n.basis.mxex = 7, base.system.mxex = "B-spline",
                 family = "gaussian", screening = TRUE){




  #Y=y;X=x;Z=z;q=4;sort.mx="ASC";mxex.smooth = FALSE;n.order.beta=4;n.basis.beta=7;base.system.beta="B-spline";n.order.mxex=4;n.basis.mxex=7;base.system.mxex="B-spline";family="binomial"
  # preliminary check
  if (sum(is.na(Y))>0) stop("'Y' must be complete vector")
  if (sum(is.na(X)) > 0) stop("'X' must be complete data frame with no missing values.")
  if (!is.null(Z)){
    if (!(length(Y) == nrow(X) && nrow(X) == nrow(Z))) stop("Please check if Y, X, and Z contain the same sample.")
  }

  # Options to perform hierarchical clustering to sort mixture exposure substance
  if (sort.mx == "CR-direct") {
    cor_matrix <- cor(X)
    hc <- hclust(as.dist(1 - cor_matrix))
    sort.list <- colnames(X)[hc$order]
    X.sort <- X[,sort.list]
  } else if (sort.mx == "CR-abs") {
    cor_matrix <- abs(cor(X))
    hc <- hclust(as.dist(1 - cor_matrix))
    sort.list <- colnames(X)[hc$order]
    X.sort <- X[,sort.list]
  } else if (sort.mx == "ASC") {
    cor.coe <- sapply(X, function(x) cor(Y, x)) %>%
      sort(decreasing = TRUE)
    X.sort <- X[, names(cor.coe)]
  } else {
    X.sort <- X
  }

  # convert X into quantiles
  X.sort.qtl <- X.sort %>%
    dplyr::mutate(across(where(is.numeric), ~{
      q <- unique(quantile(., probs = seq(0, 1, 1/q), na.rm = TRUE))
      cut(., breaks = q, labels = FALSE, include.lowest = TRUE)
    })) %>%
    '-'(1) %>%
    as.matrix()

  #QR decomposition
  dqr <- qr(X.sort.qtl)
  cpn <- dqr$pivot[1:dqr$rank]
  X3 <- X.sort.qtl[, cpn]
  pos <- c(1:ncol(X3))[cpn]
  nsample <- nrow(X3)
  nexp   <- ncol(X3)

  # normalize the position information to [0, 1]
  pos2 <- (pos - min(pos)) / (max(pos) - min(pos))

  if (mxex.smooth){
    # create spline function to fit mixture effect function
    if (base.system.beta == "B-spline"){
      betabasis <- create.bspline.basis(norder = n.order.beta, nbasis = n.basis.beta)
    } else if (base.system.beta == "Fourier"){
      betabasis  <- create.fourier.basis(c(0,1), nbasis = n.basis.beta)
    }

    # create spline function to fit mixture exposure function
    if (base.system.mxex == "B-spline"){
      mxexbasis <- create.bspline.basis(norder = n.order.mxex, nbasis = n.basis.mxex)
    } else if (base.system.mxex == "Fourier"){
      mxexbasis  <- create.fourier.basis(c(0,1), nbasis = n.basis.mxex)
    }

    # Functional manipulation to mixture exposure substance
    B      <- eval.basis(pos, mxexbasis)
    to_mul <- ginv(t(B) %*% B) %*% t(B)
    U      <- X3 %*% t( to_mul )
    J      <- inprod(mxexbasis, betabasis)   ### added by Fan ###
    UJ     <- matrix( U %*% J, ncol = ncol(J) )
  } else {
    # create spline function to fit mixture effect function
    if (base.system.beta == "B-spline"){
      betabasis <- create.bspline.basis(norder = n.order.beta, nbasis = n.basis.beta)
    } else if (base.system.beta == "Fourier"){
      betabasis  <- create.fourier.basis(c(0,1), nbasis = n.basis.beta)
    }

    # Calculate W (please refer to the paper)
    B <- eval.basis(pos2, betabasis)
    UJ <- X3 %*% B
  }

  # collect output results
  output.result <- list()

  # build formula
  if (is.null(Z)) {
    reduce.list <- c("Y ~ 1")
    full.list <- c("Y ~ 1 + UJ")
  } else {
    reduce.list <- c("Y ~ Z[,1]")
    for (i in 2: ncol(Z)){
      reduce.list <- paste(reduce.list, "+ Z[,", i, "]")
    }
    full.list <- paste(reduce.list, "+ UJ")
  }

  # Calculate p values for LRT, F and global test
  lm.fit.full <- glm(as.formula(full.list), family = family)
  if (family == "gaussian") {
    gt.fit <- gt(as.formula(reduce.list), as.formula(full.list), model = "linear")
    output.result$p.value <- c(anova(lm.fit.full, test = "LRT")[(2+ncol(Z)),5], anova(lm.fit.full, test = "F")[(2+ncol(Z)),6], p.value(gt.fit))
    names(output.result$p.value) <- c("Likelihood ratio test", "F distribution test", "Global test")
  } else if (family == "binomial") {
    gt.fit <- gt(as.formula(reduce.list), as.formula(full.list), model = "logistic")
    output.result$p.value <- c(anova(lm.fit.full, test = "LRT")[(2+ncol(Z)),5], p.value(gt.fit))
    names(output.result$p.value) <- c("Likelihood ratio test", "Global test")
  }

  cef.list <- summary(lm.fit.full)$coefficient %>%
    '['(c((nrow(.)-n.basis.beta+1):nrow(.)), 1)
  Beta.t <- B %*% cef.list
  output.result$est.effect <- c(sum(Beta.t), sum(abs(Beta.t)), sum(Beta.t[Beta.t>0]))
  names(output.result$est.effect) <- c("Direct sum", "Absolute sum", "Positive sum")

  # Selection for output
  if (screening) {
    return(output.result)
  } else {
    # Calculate effect size
    lm.fit.part <- glm(as.formula(reduce.list), family = family)
    SSuj <- sum((fitted(lm.fit.full) - fitted(lm.fit.part))^2)
    SSEfull <- sum(residuals(lm.fit.full)^2)
    output.result$effect.size <- eta_square <- SSuj/(SSuj + SSEfull)

    # Plot the effect function beta(t)
    plot.dt <- data.frame(x = 1:ncol(X.sort.qtl),
                          y = Beta.t)

    output.result$beta.t.plot <-ggplot(plot.dt, aes(x = x, y = y)) +
      geom_point() +
      geom_smooth(method = "loess", se = FALSE) +
      geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
      xlab("Exposure substance number") +
      ylab("Exposure effect") +
      ggtitle("Beta(t)")

    # Ordered effect estimate
    ord.eff <- data.frame(exp.name = colnames(X.sort.qtl),
                          exp.coef = Beta.t,
                          exp.effc = abs(Beta.t))
    ord.eff$exp.name <- factor(ord.eff$exp.name, levels = ord.eff$exp.name)
    output.result$effect.plot <- ggplot(ord.eff, aes(x = exp.name, y = exp.coef, fill = exp.coef > 0)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("#7c9d97","#e9b383"), labels = c("Negative", "Positive"), name = "Coefficient Sign") +
      coord_flip() +
      labs(x = "Exposure Name", y = "Exposure Effect", title = "Mixture Exposure Effect Estimated by Beta(x)") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 7),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 9))
    return(output.result)
  }
}

