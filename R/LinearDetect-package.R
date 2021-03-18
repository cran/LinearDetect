## usethis namespace: start
#' @useDynLib LinearDetect, .registration = TRUE
## usethis namespace: end
NULL
#' Plot the cross-validation score
#'
#' @param pred.error prediction error
#' @param lambda indice of tuning parameter lambda
#' @return No return value, called for plot
#' @import graphics
mspe.plot <- function(pred.error, lambda){
    plot( lambda, pred.error, type = 'o', col = "blue")
}

#' BIC  and HBIC function
#'
#' @param residual residual matrix
#' @param phi estimated coefficient matrix of the model
#' @param gamma.val hyperparameter for HBIC, if HBIC == TRUE.
#' @param method method name for the model: MLR: Multiple Linear Regression; VAR: Vector autoregression;
#' @return A list object, which contains the followings
#' \describe{
#'   \item{BIC}{BIC value}
#'   \item{HBIC}{HBIC value}
#' }
BIC <- function(residual, phi, gamma.val = 1, method = 'MLR'){
  p.y <- length(phi[, 1]);
  p.x <- length(phi[1, ]);
  T.new <- length(residual[1, ]);
  # print("p.x"); print(p.x);
  # print("p.y"); print(p.y);
  # print("T.new"); print(T.new);
  # count : non-zero coefficient
  count <- 0;
  for (i in 1:p.y){
    for (j in 1:p.x){
      if(phi[i,j] != 0){
        count <- count + 1;
      }
    }
  }
  # print("nonzero count"); print(count)

  sigma.hat <- 0*diag(p.y);
  for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[, i]%*%t(residual[, i]);  }
  sigma.hat <- (1/(T.new))*sigma.hat;
  ee.temp <- min(eigen(sigma.hat)$values);
  if(ee.temp <= 10^(-8)){
    # print("nonpositive eigen values!")
    sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp) + 10^(-3))*diag(p.y);
  }

  log.det <- log(det(sigma.hat));
  # print(log.det)
  # print(log(T.new)*count/T.new)
  # print(2*gamma.val*log(p.x*p.y)*count/T.new)
  if(method == 'VAR'){
    # print("this is only for VAR model!!!!! count/p")
    count <- count/p.y
  }
  return(list(BIC = log.det + log(T.new)*count/T.new , HBIC = log.det + 2*gamma.val*log(p.x*p.y)*count/T.new))
}

#' Generate the linear regression model data with break points
#'
#' @param nobs number of time points
#' @param px the number of features
#' @param cnst the constant
#' @param phi parameter coefficient matrix of the linear model
#' @param sigma covariance matrix of the white noise
#' @param sigma_x variance of the predictor variable x
#' @param brk vector of break points
#' @return A list object, which contains the followings
#' \describe{
#'   \item{series_y}{matrix of response data}
#'   \item{series_x}{matrix of predictor data}
#'   \item{noises}{matrix of white noise error}
#' }
#' @import mvtnorm
#' @export
lm.sim.break <- function (nobs, px, cnst = NULL, phi = NULL, sigma, sigma_x =1, brk = nobs+1) {
    if (!is.matrix(sigma))
        sigma = as.matrix(sigma)
    #k = nrow(sigma)
    p.y <-  nrow(sigma)
    p.x <- px
    m <- length(brk)
    nT <- nobs

    # error term
    data_e = rmvnorm(nT, rep(0, p.y), sigma)
    # data x
    data_x = rmvnorm(nT, rep(0, p.x), sigma_x*diag(p.x))
    # data y
    data_y = matrix(0, nT, p.y)
    if (length(cnst) == 0)
        cnst = rep(0, p.y)
    if (m == 1){
        for (it in 1:nT) {
            tmp = matrix(data_e[it, ], 1, p.y)
            tmp_x = matrix(data_x[it, ], 1, p.x)
            phj = phi[, 1:p.x]
            if (p.y == 1){
                tmp = tmp + tmp_x %*% phj
            }else{
                tmp = tmp + tmp_x %*% t(phj)
            }
            data_y[it, ] = cnst + tmp
        }
    }

    if (m > 1){
        for (it in 1:(brk[1]-1)) {
            tmp = matrix(data_e[it, ], 1, p.y)
            tmp_x = matrix(data_x[it, ], 1, p.x)
            #idx = (i - 1) * p.x
            phj = phi[, 1:p.x]
            if (p.y == 1){
                tmp = tmp + tmp_x %*% phj
            }else{
                tmp = tmp + tmp_x %*% t(phj)
            }
            #tmp = tmp + tmp_x %*% t(phj)
            data_y[it, ] = cnst + tmp
        }
        for ( mm in 1:(m-1)){
            for (it in (brk[mm]):(brk[mm+1]-1) ) {
                tmp = matrix(data_e[it, ], 1, p.y)
                tmp_x = matrix(data_x[it, ], 1, p.x)
                phj = phi[, (mm*p.x + 1):(mm*p.x+ p.x)]
                if (p.y == 1){
                    tmp = tmp + tmp_x %*% phj
                }else{
                    tmp = tmp + tmp_x %*% t(phj)
                }
                #tmp = tmp + tmp_x %*% t(phj)
                data_y[it, ] = cnst + tmp
            }
        }
    }

    data_y = data_y[1:nT, ]
    data_x = data_x[1:nT, ]
    data_e = data_e[1:nT, ]
    lmsim <- list(series_y = data_y, series_x = data_x, noises = data_e)
}


#' Generate the constant model data with break points
#'
#' @param nobs number of time points
#' @param cnst the constant
#' @param sigma covariance matrix of the white noise
#' @param brk vector of break points
#' @return A list object, which contains the followings
#' \describe{
#'   \item{series_y}{matrix of response data}
#'   \item{noises}{matrix of white noise error}
#' }
#' @import mvtnorm
#' @export
constant.sim.break <- function (nobs, cnst, sigma, brk = nobs+1) {
  if (!is.matrix(sigma))
    sigma = as.matrix(sigma)
  p.y <-  nrow(sigma)
  m <- length(brk)
  nT <- nobs

  # error term
  data_e = rmvnorm(nT, rep(0, p.y), sigma)
  # data y
  data_y = matrix(0, nT, p.y)
  if (length(cnst) == 0)
    cnst = rep(0, p.y)
  if (m == 1){
    for (it in 1:nT) {
      tmp = matrix(data_e[it, ], 1, p.y)
      cnst.temp = cnst
      data_y[it, ] = cnst.temp + tmp
    }
  }

  if (m > 1){
    for (it in 1:(brk[1]-1)) {
      tmp = matrix(data_e[it, ], 1, p.y)
      cnst.temp = cnst[, 1]
      data_y[it, ] = tmp + cnst.temp
    }
    for ( mm in 1:(m-1)){
      for (it in (brk[mm]):(brk[mm+1]-1) ) {
        tmp = matrix(data_e[it, ], 1, p.y)
        cnst.temp = cnst[, mm + 1]
        data_y[it, ] = tmp + cnst.temp
      }
    }
  }

  data_y = data_y[1:nT, ]
  data_e = data_e[1:nT, ]
  cnstsim <- list(series_y = data_y, noises = data_e)
}


#' Generate the gaussian graphical model data with break points
#'
#' @param nobs number of time points
#' @param px the number of features
#' @param sigma covariance matrix of the X matrix
#' @param brk vector of break points
#' @return A list object, which contains the followings
#' \describe{
#'   \item{series_x}{matrix of data}
#' }
#' @import mvtnorm
#' @export
ggm.sim.break <- function (nobs, px, sigma, brk = nobs+1) {
  if (!is.matrix(sigma))
    sigma = as.matrix(sigma)
  p.x <- px
  m <- length(brk)
  nT <- nobs

  sigma.temp = sigma[, 1:p.x]
  data_x = rmvnorm(nT, rep(0, p.x), sigma.temp)
  if (m == 1){
    sigma.temp = sigma[, 1:p.x]
    data_x = rmvnorm(nT, rep(0, p.x), sigma.temp)
  }

  if (m > 1){
    for (it in 1:(brk[1]-1)) {
      sigma.temp = sigma[, 1:p.x]
      data_x[it, ] = rmvnorm(1, rep(0, p.x), sigma.temp)
    }
    for ( mm in 1:(m-1)){
      for (it in (brk[mm]):(brk[mm+1]-1) ) {
        sigma.temp = sigma[, (mm*p.x + 1):(mm*p.x+ p.x)]
        data_x[it, ] = rmvnorm(1, rep(0, p.x), sigma.temp)
      }
    }
  }

  data_x = data_x[1:nT, ]
  ggmsim <- list(series_x = data_x)
}



################################################################
#' Threshold block fused lasso (TBFL) algorithm for change point detection
#'
#' @description Perform the threshold block fused lasso (TBFL) algorithm to detect the structural breaks
#' in large scale high-dimensional non-stationary linear regression models.
#'
#' @param method method name for the model: Constant: Mean-shift Model; MvLR: Multivariate Linear Regression; MLR: Multiple Linear Regression; VAR: Vector autoregression; GGM: Gaussian graphical model
#' @param data_y input data matrix (response), with each column representing the time series component
#' @param data_x input data matrix (predictor), with each column representing the time series component
#' @param lambda.1.cv tuning parmaeter lambda_1 for fused lasso
#' @param lambda.2.cv tuning parmaeter lambda_2 for fused lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param block.size the block size
#' @param blocks the blocks
#' @param refit logical; if TRUE, refit the model, if FALSE, use BIC to find a thresholding value and then output the parameter estimates without refitting. Default is FALSE.
#' @param fixed_index index for linear regression model with only partial compoenents change.
#' @param HBIC logical; if TRUE, use high-dimensional BIC, if FALSE, use orginal BIC. Default is FALSE.
#' @param gamma.val hyperparameter for HBIC, if HBIC == TRUE.
#' @param optimal.block logical; if TRUE, grid search to find optimal block size, if FALSE, directly use the default block size. Default is TRUE.
#' @param optimal.gamma.val hyperparameter for optimal block size, if optimal.blocks == TRUE. Default is 1.5.
#' @param block.range the search domain for optimal block size.
#' @return A list object, which contains the followings
#' \describe{
#'   \item{cp.first}{a set of selected break point after the first block fused lasso step}
#'   \item{cp.final}{a set of selected break point after the final exhaustive search step}
#'   \item{beta.hat.list}{a list of estimated parameter coefficient matrices for each stationary segementation}
#'   \item{beta.est}{a list of estimated parameter coefficient matrices for each block}
#'   \item{beta.final}{a list of estimated parameter coefficient matrices for each stationary segementation, using BIC thresholding or refitting the model.}
#'   \item{beta.full.final}{For GGM only. A list of \eqn{p \times p} matrices for each stationary segementation. The off-diagonal entries are same as the beta.final.}
#'   \item{jumps}{The change (jump) of the values in estimated parameter coefficient matrix.}
#'   \item{bn.optimal}{The optimal block size.}
#'   \item{bn.range}{The values of block size in grid search.}
#'   \item{HBIC.full}{The HBIC values.}
#'   \item{pts.full}{The selected change points for each block size.}
#' }
#' @author Yue Bai, \email{baiyue@@ufl.edu}
#' @importFrom Rcpp sourceCpp
#' @importFrom glmnet cv.glmnet
#' @importFrom sparsevar fitVAR
#' @export
#' @examples
#' #### constant model
#' TT <- 10^3; # number of observations/samples
#' p.y <- 50; # dimension of observed Y
#' brk <- c(floor(TT/3),floor(2*TT/3), TT+1)
#' m <- length(brk)
#' d <- 5 #number of non-zero coefficient
#' ### generate coefficient
#' constant.full <- matrix(0, p.y, m)
#' set.seed(1)
#' constant.full[sample(1:p.y, d, replace = FALSE), 1] <- runif(d, -1, -0.5);
#' constant.full[sample(1:p.y, d, replace = FALSE), 2] <- runif(d, 0.5, 1);
#' constant.full[sample(1:p.y, d, replace = FALSE), 3] <- runif(d, -1, -0.5);
#' e.sigma <- as.matrix(1*diag(p.y))
#' try <- constant.sim.break(nobs = TT, cnst = constant.full, sigma = e.sigma, brk = brk)
#' data_y <- try$series_y; data_y <- as.matrix(data_y, ncol = p.y)
#' ### Fit the model
#' method <- c("Constant")
#' temp <- tbfl(method, data_y, block.size = 40, optimal.block = FALSE) #use a single block size
#' temp$cp.final
#' temp$beta.final
#' \donttest{temp <- tbfl(method, data_y) # using optimal block size}
#'
#'
#'
#'
#' #### multiple linear regression
#' TT <- 2*10^3; # number of observations/samples
#' p.y <- 1; # dimension of observed Y
#' p.x <- 20
#' brk <- c(floor(TT/4), floor(2*TT/4), floor(3*TT/4), TT+1)
#' m <- length(brk)
#' d <- 15 #number of non-zero coefficient
#' ###generate coefficient beta
#' beta.full <- matrix(0, p.y, p.x*m)
#' set.seed(1)
#' aa <- c(-3, 5, -3, 3)
#' for(i in 1:m){beta.full[1, (i-1)*p.x+sample(1:p.x, d, replace = FALSE)] <- aa[i] + runif(d, -1, 1);}
#' e.sigma <- as.matrix(1*diag(p.y))
#' try <- lm.sim.break(nobs = TT, px = p.x, phi = beta.full, sigma = e.sigma, sigma_x = 1, brk = brk)
#' data_y <- try$series_y; data_y <- as.matrix(data_y, ncol = p.y)
#' data_x <- try$series_x; data_x <- as.matrix(data_x)
#' ### Fit the model
#' method <- c("MLR")
#' \donttest{temp <- tbfl(method, data_y, data_x)}
#' \donttest{temp$cp.final} #change points
#' \donttest{temp$beta.final} #final estimated parameters (after BIC threshold)
#' \donttest{temp_refit <- tbfl(method, data_y, data_x, refit = TRUE)}
#' \donttest{temp_refit$beta.final} #final estimated parameters (refitting the model)
#'
#'
#'
#'
#' #### Gaussian Graphical model
#' TT <- 3*10^3; # number of observations/samples
#' p.x <- 20 # dimension of obsrved X
#' # TRUE BREAK POINTS WITH T+1 AS THE LAST ELEMENT
#' brk <- c(floor(TT/3), floor(2*TT/3), TT+1)
#' m <- length(brk)
#' ###generate precision matrix and covariance matrix
#' eta = 0.1
#' d <- ceiling(p.x*eta)
#' sigma.full <- matrix(0, p.x, p.x*m)
#' omega.full <- matrix(0, p.x, p.x*m)
#' aa <- 1/d
#' for(i in 1:m){
#' if(i%%2==1){
#' ajmatrix <- matrix(0, p.x, p.x)
#' for(j in 1:(floor(p.x/5)) ){
#' ajmatrix[ ((j-1)*5+1): (5*j), ((j-1)*5+1): (5*j)] <- 1
#' }
#' }
#' if(i%%2==0){
#' ajmatrix <- matrix(0, p.x, p.x)
#' for(j in 1:(floor(p.x/10)) ){
#' ajmatrix[ seq(((j-1)*10+1), (10*j), 2), seq(((j-1)*10+1), (10*j), 2)] <- 1
#' ajmatrix[ seq(((j-1)*10+2), (10*j), 2), seq(((j-1)*10+2), (10*j), 2)] <- 1
#' }
#' }
#' theta <- aa* ajmatrix
#' # force it to be positive definite
#' if(min(eigen(theta)$values) <= 0){
#' print('add noise')
#' theta = theta - (min(eigen(theta)$values)-0.05) * diag(p.x)
#' }
#' sigma.full[, ((i-1)*p.x+1):(i*p.x)]  <- as.matrix(solve(theta))
#' omega.full[, ((i-1)*p.x+1):(i*p.x)]  <- as.matrix(theta)
#' }
#' # simulate data
#' try <- ggm.sim.break(nobs = TT, px = p.x, sigma = sigma.full, brk = brk)
#' data_y <- try$series_x; data_y <- as.matrix(data_y)
#' ### Fit the model
#' method <- c("GGM")
#' #use a single block size
#' \donttest{temp <- tbfl(method,data_y = data_y,block.size = 80,optimal.block = FALSE)}
#' \donttest{temp$cp.final #change points}
#' \donttest{temp$beta.final}
#'
tbfl <- function(method, data_y, data_x = NULL, lambda.1.cv = NULL, lambda.2.cv = NULL, q = 1,
                 max.iteration = 100, tol = 10^(-2), block.size = NULL, blocks = NULL, refit = FALSE,
                 fixed_index = NULL, HBIC = FALSE, gamma.val = NULL, optimal.block = TRUE, optimal.gamma.val = 1.5, block.range = NULL){
  method.full <- c("Constant","MvLR", "MLR", "VAR", "GGM", "MvLR_fixed", "MLR_fixed");
  if ( !(method %in% method.full) ){
    stop("Incorrect method name!")
  }
  TT  <- length(data_y[, 1]);
  ############# block size and blocks ###########
  if(method == 'Constant'){
    p.y <- length(data_y[1, ]);
    data_x <- matrix(1, nrow = TT)
    p.x <- 1
  }
  if(method == 'GGM'){
    p.y <- length(data_y[1, ]);
    data_x <- data_y
    p.x <- p.y
  }
  if (method == 'MvLR' |  method =="MLR"){
    if (is.null(data_x)){
      stop("Empty predictor data!")
    }
    p.y <- length(data_y[1, ]); p.x <- length(data_x[1, ]);
  }
  if(method == 'VAR'){
    p <- length(data_y[1,]);
    p.y <- p
    p.x <- p
  }
  if( (!is.null(block.size) || !is.null(blocks)) && optimal.block == TRUE){
    stop("Customed blocksize/blocks setting found. Set optimal.block to be FALSE!")
  }
  if(optimal.block == TRUE){
    if(method == 'VAR'){
      stop("For VAR model, set optimal.block = FALSE!")
    }
    if(!is.null(block.range)){
        b_n.range <- block.range[block.range > 1]
        print("block size:")
        print(b_n.range)

    }else{
        if(method == 'GGM'){
            if(sqrt(TT) > p.x){
                b_n.max <- ceiling(min(sqrt(TT), TT/20))
            }else{
                b_n.max <- ceiling(min(sqrt(TT)*log(p.x), TT/20))
            }
            b_n.min <- floor(min(log(TT)*log(p.x), TT/20 ))
            print("maximum block size:")
            print(b_n.max)
            print("minimum block size:")
            print(b_n.min)
            b_n.range <- round(seq(b_n.min, b_n.max, length.out = 5))
            # b_n.range <- round(seq(b_n.min, b_n.max, length.out = 5)/5)*5
            b_n.range <- unique(b_n.range)
            b_n.range <- b_n.range[b_n.range > 1]
            print("block size:")
            print(b_n.range)

        }else{
            if(sqrt(TT) > p.x*p.y){
                b_n.max <- ceiling(min(sqrt(TT), TT/20))
            }else{
                b_n.max <- ceiling(min(sqrt(TT)*log(p.x*p.y), TT/20))
            }
            b_n.min <- floor(min(log(TT)*log(p.x*p.y), TT/20 ))
            print("maximum block size:")
            print(b_n.max)
            print("minimum block size:")
            print(b_n.min)
            b_n.range <- round(seq(b_n.min, b_n.max, length.out = 5))
            # b_n.range <- round(seq(b_n.min, b_n.max, length.out = 5)/5)*5
            b_n.range <- unique(b_n.range)
            b_n.range <- b_n.range[b_n.range > 1]
            print("block size:")
            print(b_n.range)
        }

    }


    n.method <- length(b_n.range) #number fo methods
    temp.full <- vector("list", n.method);
    pts.full <- vector("list", n.method);
    for(j.2 in 1:n.method){
      block.size = b_n.range[j.2]
      b_n_bound = 2*block.size  #block size for boundary
      blocks <- c(seq(1, b_n_bound*2+1 , b_n_bound),
                  seq(b_n_bound*2+block.size+1, TT+1-2*b_n_bound, block.size),
                  seq(TT+1-b_n_bound, TT+1,  b_n_bound));
      # blocks <- seq(1, TT + 1, block.size);
      if(blocks[length(blocks)] < TT+1){
        blocks <- c(blocks[-length(blocks)], TT + 1)
      }

      n.new <- length(blocks) - 1;
      blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj + 1] - blocks[jjj]  );
      print("block sizes")
      print(blocks.size)

      #sample the cv index for cross-validation
      #bbb <- floor(n.new/5);
      #aaa <- sample(1:5, 1);
      bbb <- floor(n.new/4);
      aaa <- 4;
      cv.index <- seq(aaa,n.new,floor(n.new/bbb));
      cv.l <- length(cv.index);
      print("cv.index:"); print(cv.index)

      if(method == 'Constant'){
        p.y <- length(data_y[1, ]);
        data_x <- matrix(1, nrow = TT)
        p.x <- 1
        ############# Tuning parameter ################
        if(is.null(lambda.1.cv)){
          lambda.1.max <- lambda_warm_up_lm(data_y, data_x, blocks, cv.index)$lambda_1_max
          if(blocks[2] <= (p.x + p.y) ){
            epsilon <-  10^(-3)
          }
          if(blocks[2] >= (p.x + p.y) ){
            epsilon <-  10^(-4)
          }
          nlam <- 10
          lambda.1.min <-  lambda.1.max*epsilon
          delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
          lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
        }

        if(is.null(lambda.2.cv)){
          lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
        }

        ####################################################
        ########## first step       #######################
        ####################################################
        temp.first <- lm.first.step.blocks(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration, tol = tol, blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
        cp.first <- temp.first$pts.list;
        cl.number <- length(cp.first);
        beta.est <- temp.first$beta.full
        temp.first.all<- temp.first
        jumps.l2 <- temp.first$jumps.l2

        ####################################################
        ########## second step       #######################
        ####################################################
        if(length(cp.first) > 0){
          temp.second<- lm.second.step.search(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
          temp.second.all<- temp.second
          cp.final<- temp.second$cp.final;
          beta.hat.list <- temp.second$beta.hat.list
        }else{
          cp.final <- c()
          print('no change points!')
          beta.hat.list <- beta.est[[floor(n.new/2)]]
        }

        if(refit){
          print('refit the model!')
          cp.full <- c(1, cp.final, TT+1)
          m <- length(cp.final) + 1
          temp_beta <-  list(m)
          for(i in 1:m){
            # data_x_temp <- matrix(1, nrow = TT)
            data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
            temp_beta[[i]] <- apply(data_y_temp, 2, mean)
          }

          beta.final.temp <- matrix(c(do.call("cbind", temp_beta)), nrow = 1)
          lambda.val.best <- BIC.threshold(method, beta.final.temp, p.y,
                                           length(c(cp.final, TT+1)), c(cp.final, TT+1), data_y , b_n = block.size)

          for(j in 1:length(c(cp.final, TT+1))){
            temp_beta[[j]][abs(temp_beta[[j]]) < lambda.val.best[j]] <- 0
          }
          beta.final <- temp_beta

        }else{
          print('no refit!')
          # BIC threholding step!!!
          beta.final.temp <- matrix(c(do.call("cbind", beta.hat.list)), nrow= 1)
          lambda.val.best <- BIC.threshold(method, beta.final.temp, p.y, length(cp.final) + 1, c(cp.final, TT+1),
                                           data_y,  b_n = block.size, nlam = 20)
          # print(lambda.val.best)
          temp <- beta.hat.list
          for(j in 1:(length(cp.final) + 1) ){
            temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
          }
          beta.final <- temp

        }

        temp.full[[j.2]] <- beta.final
        pts.full[[j.2]] <- cp.final
        # return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
        #             beta.est = beta.est, beta.final = beta.final, jumps = jumps.l2))
        # return(list(cp.first = cp.first, cp.final = cp.final, beta.hat.list = beta.hat.list, beta.est = beta.est ))


      }
      if(method == 'GGM'){
        p.y <- length(data_y[1, ]);
        data_x <- data_y
        p.x <- p.y
        ############# Tuning parameter ################
        if(is.null(lambda.1.cv)){
          lambda.1.max <- lambda_warm_up_lm(data_y, data_x, blocks, cv.index)$lambda_1_max
          if(blocks[2] <= (p.x + p.y) ){
            epsilon <-  10^(-3)
          }
          if(blocks[2] >= (p.x + p.y) ){
            epsilon <-  10^(-4)
          }
          nlam <- 10
          lambda.1.min <-  lambda.1.max*epsilon
          delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
          lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
        }

        if(is.null(lambda.2.cv)){
          lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
        }

        ####################################################
        ########## first step       #######################
        ####################################################
        temp.first <- ggm.first.step.blocks(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration, tol = tol, blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
        cp.first <- temp.first$pts.list;
        cl.number <- length(cp.first);
        beta.est <- temp.first$beta.full
        temp.first.all<- temp.first
        jumps.l2 <- temp.first$jumps.l2


        ####################################################
        ########## second step       #######################
        ####################################################
        if(length(cp.first) > 0){
          temp.second<- ggm.second.step.search(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
          temp.second.all<- temp.second
          cp.final<- temp.second$cp.final;
          beta.hat.list <- temp.second$beta.hat.list
        }else{
          cp.final <- c()
          print('no change points!')
          beta.hat.list <- beta.est[[floor(n.new/2)]]

        }

        if(refit){
          print('refit the model!')
          temp <- beta.hat.list
          cp.full <- c(1, cp.final, TT+1)
          for(i in 1:(length(cp.final) + 1) ){
            data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
            tmp_coef <- c()
            for(j.1 in 1:p.y){
              data_x_temp <- data_y_temp[, -j.1]
              cvfit = cv.glmnet(data_x_temp, data_y_temp[, j.1], intercept = FALSE)
              opt_lambda  <- cvfit$lambda.1se  # Optimal Lambda
              # fit = glmnet(data_x_temp, data_y_temp, intercept = FALSE, lambda = opt_lambda)
              tmp_coef <- rbind(tmp_coef, coef(cvfit)[-1])

            }
            rownames(tmp_coef) <- NULL
            temp[[i]] <- tmp_coef


          }
          beta.final <- temp

        }else{
          print('no refit!')
          # BIC threholding step!!!
          beta.final.temp <- as.matrix(do.call("cbind", beta.hat.list))
          lambda.val.best <- BIC.threshold.ggm(beta.final.temp, p.x-1, length(cp.final) + 1, c(cp.final, TT+1),
                                               data_y, b_n = block.size, nlam = 50 )$lambda.val.best
          temp <- beta.hat.list
          for(j in 1:(length(cp.final) + 1)){
            temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
          }
          beta.final <- temp

        }

        m <- length(cp.final) + 1
        beta.full.final <- vector("list", m);
        for(i in 1:m){
          beta.full.final[[i]] <- matrix(0, p.y, p.y)
          beta.full.final[[i]][1, 2:p.y] <- beta.final[[i]][1, 1:(p.y-1) ]
          beta.full.final[[i]][p.y,  1:(p.y-1)] <- beta.final[[i]][p.y, 1:(p.y-1)]
          for(j in 2:(p.y-1)){
            beta.full.final[[i]][j, (j+1):p.y] <- beta.final[[i]][j, (j):(p.y-1) ]
            beta.full.final[[i]][j, 1:(j-1)] <- beta.final[[i]][j, 1:(j-1) ]
          }
        }

        temp.full[[j.2]] <- beta.final
        pts.full[[j.2]] <- cp.final
        # return(list(cp.first = cp.first, cp.final = cp.final, beta.hat.list = beta.hat.list ))
        # return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
        #             beta.est = beta.est, beta.final = beta.final, beta.full.final= beta.full.final, jumps = jumps.l2))
      }
      if (method == 'MvLR' |  method =="MLR"){
        if (is.null(data_x)){
          stop("Empty predictor data!")
        }
        p.y <- length(data_y[1, ]); p.x <- length(data_x[1, ]);

        ############# Tuning parameter ################
        if(is.null(lambda.1.cv)){
          lambda.1.max <- lambda_warm_up_lm(data_y, data_x, blocks, cv.index)$lambda_1_max
          if(blocks[2] <= (p.x + p.y) ){
            epsilon <-  10^(-3)
          }
          if(blocks[2] >= (p.x + p.y) ){
            epsilon <-  10^(-4)
          }
          nlam <- 10
          lambda.1.min <-  lambda.1.max*epsilon
          delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
          lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
        }

        if(is.null(lambda.2.cv)){
          lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
        }

        ####################################################
        ########## first step       #######################
        ####################################################
        temp.first <- lm.first.step.blocks(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration, tol = tol, blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
        cp.first <- temp.first$pts.list;
        cl.number <- length(cp.first);
        beta.est <- temp.first$beta.full
        temp.first.all<- temp.first
        jumps.l2 <- temp.first$jumps.l2
        ####################################################
        ########## second step       #######################
        ####################################################
        if(length(cp.first) > 0){
          temp.second<- lm.second.step.search(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
          temp.second.all<- temp.second
          cp.final<- temp.second$cp.final;
          beta.hat.list <- temp.second$beta.hat.list
        }else{
          cp.final <- c()
          print('no change points!')
          beta.hat.list <- beta.est[[floor(n.new/2)]]
        }

        #final estimates ( BIC thresholding or refitting the model)
        if(refit){
          print('refit the model!')
          temp <- beta.hat.list
          cp.full <- c(1, cp.final, TT+1)
          for(i in 1:(length(cp.final) + 1) ){
            data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
            data_x_temp <- as.matrix(data_x[cp.full[i]: (cp.full[i+1]-1), ])
            if(method == 'MvLR'){
              cvfit = cv.glmnet(data_x_temp, data_y_temp, intercept = FALSE,
                                family = "mgaussian")
              opt_lambda  <- cvfit$lambda.1se  # Optimal Lambda
              tmp_coeffs <- coef(cvfit, s = "lambda.1se")
              tmp_coeffs <- unlist(tmp_coeffs, use.names = FALSE)
              tmp_coef <- tmp_coeffs[[1]][-1]
              if(length(tmp_coeffs) >= 2){
                for(j in 2: length(tmp_coeffs)){
                  tmp_coef <- rbind(tmp_coef, tmp_coeffs[[j]][-1])
                }
              }

              rownames(tmp_coef) <- NULL
              temp[[i]] <- tmp_coef
            }else{
              cvfit = cv.glmnet(data_x_temp, data_y_temp, intercept = FALSE)
              opt_lambda  <- cvfit$lambda.1se  # Optimal Lambda
              # fit = glmnet(data_x_temp, data_y_temp, intercept = FALSE, lambda = opt_lambda)
              temp[[i]] <- as.matrix(coef(cvfit)[-1])
            }


          }
          beta.final <- temp

        }else{
          print('no refit!')
          # BIC threholding step!!!
          if(method == 'MLR'){
            beta.final.temp <- matrix(c(do.call("cbind", beta.hat.list)), nrow= 1)
            print(block.size)
            lambda.val.best <- BIC.threshold(method, beta.final.temp, p.x, length(cp.final) + 1, c(cp.final, TT+1),
                                             data_y, data_x,  b_n = block.size, nlam = 20)

          }else{
            print(beta.hat.list)
            beta.final.temp <- as.matrix(do.call("cbind", beta.hat.list))
            print(dim(beta.final.temp))
            lambda.val.best <- BIC.threshold(method, beta.final.temp, p.x, length(cp.final) + 1, c(cp.final, TT+1),
                                             data_y, data_x,  b_n = block.size, nlam = 20)

          }

          # print(lambda.val.best)
          temp <- beta.hat.list
          for(j in 1:(length(cp.final) + 1) ){
            temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
          }
          beta.final <- temp
        }


        temp.full[[j.2]] <- beta.final
        pts.full[[j.2]] <- cp.final
        # return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
        #             beta.est = beta.est, beta.final = beta.final, jumps = jumps.l2))
      }


    }

    HBIC.full <- c()
    BIC.full <- c()
    for(j.2 in 1:n.method){
      brk.temp <- pts.full[[j.2]]
      brk.full <- c(1, brk.temp, TT+1)
      m.hat <- length(brk.full) - 1

      if(method == 'Constant'){
        k <- p.y
        beta.final <- matrix(c(do.call("cbind", temp.full[[j.2]])), nrow = 1)
        HBIC.res <- c()
        # BIC.res <- c()
        residual.full <- c()
        for(i in 1:m.hat){
          beta.temp <- beta.final[((i-1)*k+1):(i*k)]
          data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
          data_x.temp <- matrix(1, nrow = brk.full[i+1] - brk.full[i])
          data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, nrow = 1)
          residual.temp <- data_y.temp - data_y.est
          residual.full <- rbind(residual.full, residual.temp)
        }

        phi = do.call("cbind", temp.full[[j.2]])
        residual <- t(residual.full)
        p.y <- length(phi[, 1]);
        p.x <- length(phi[1, ]);
        T.new <- length(residual[1, ]);
        count <- 0;
        for (i in 1:p.y){
          for (j in 1:p.x){
            if(phi[i,j] != 0){
              count <- count + 1;
            }
          }
        }

        sigma.hat <- 0*diag(p.y);
        for(i in 1:T.new){sigma.hat <- sigma.hat +  residual[, i]%*%t(residual[, i]);  }
        sigma.hat <- (1/(T.new))*sigma.hat;
        ee.temp <- min(eigen(sigma.hat)$values);
        if(ee.temp <= 10^(-8)){
          print("nonpositive eigen values!")
          sigma.hat <- sigma.hat + (2.0)*(abs(ee.temp) + 10^(-3))*diag(p.y);
        }

        log.det <- log(det(sigma.hat));
        HBIC.res <- log.det*TT + 2*optimal.gamma.val*log(p.y)*count
        # BIC.res <- log.det*TT + (count+3*(p.x-1))*log(TT)

      }
      if(method == 'MLR'){
        k <- p.x
        beta.final <- matrix(c(do.call("cbind", temp.full[[j.2]])), nrow = 1)
        HBIC.res <- c()
        residual.full <- c()
        for(i in 1:m.hat){
          beta.temp <- beta.final[((i-1)*k+1):(i*k)]
          data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
          data_x.temp <- data_x[brk.full[i]:(brk.full[i+1]-1), ]
          data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, ncol = 1)
          residual.temp <- data_y.temp - data_y.est
          residual.full <- rbind(residual.full, residual.temp)
        }

        phi =  t(do.call("cbind", temp.full[[j.2]]))
        residual <- t(residual.full)
        p.y <- length(phi[, 1]);
        p.x <- length(phi[1, ]);
        T.new <- length(residual[1, ]);
        print("p.x"); print(p.x);
        print("p.y"); print(p.y);
        print("T.new"); print(T.new);
        # count : non-zero coefficient
        count <- 0;
        for (i in 1:p.y){
          for (j in 1:p.x){
            if(phi[i,j] != 0){
              count <- count + 1;
            }
          }
        }

        sigma.hat <- sum(residual^2)/T.new

        log.det <- log(sigma.hat);
        p.y <- 1
        HBIC.res <- log.det*TT + 2*optimal.gamma.val*log(p.x)*count

      }
      if(method == 'GGM'){
        beta.final.temp <- temp.full[[j.2]]
        k <- p.x - 1
        HBIC.res <- c()
        for(jjj in 1:p.x){
          beta.final <- matrix(c(do.call("cbind", beta.final.temp)[jjj, ]), nrow = 1)

          residual.full <- c()
          data_x <- data_y[, -jjj]
          for(i in 1:m.hat){
            beta.temp <- beta.final[((i-1)*k+1):(i*k)]
            data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), jjj ] )
            data_x.temp <- data_x[brk.full[i]:(brk.full[i+1]-1), ]
            data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, ncol = 1)
            residual.temp <- data_y.temp - data_y.est
            residual.full <- rbind(residual.full, residual.temp)
          }
          phi =  t(do.call("cbind", beta.final.temp))
          residual <- t(residual.full)
          p_y <- length(phi[, 1])
          p_x <- length(phi[1, ]);
          T.new <- length(residual[1, ]);
          # count : non-zero coefficient
          count <- 0;
          for (i in 1:p_y){
            if(phi[i, jjj] != 0){
              count <- count + 1;
            }
          }

          sigma.hat <- sum(residual^2)/T.new

          log.det <- log(sigma.hat);
          HBIC.res <- c(HBIC.res, log.det*TT + 2*optimal.gamma.val*log(p.x-1)*count)
        }
      }
      HBIC.full <- c(HBIC.full, sum(HBIC.res))

    }

    # res.HBIC <- which.min(HBIC.full)
    pts.res.HBIC<- pts.full[[which.min(HBIC.full)]]
    bn.res.HBIC <- b_n.range[which.min(HBIC.full)]
    print("HBIC.full:")
    print(HBIC.full)

    return(list(cp.final = pts.res.HBIC, beta.final = temp.full[[which.min(HBIC.full)]], bn.optimal = bn.res.HBIC,
                bn.range = b_n.range, HBIC.full = HBIC.full, pts.full = pts.full))



  }else{
    if(is.null(block.size) && is.null(blocks) ){
      # block.size = floor(sqrt(TT));
      block.size = floor(log(TT)*log(p.y*p.x));
      # blocks <- seq(1, TT + 1, block.size);
      b_n_bound = 2*block.size  #block size for boundary
      blocks <- c(seq(1, b_n_bound*2+1 , b_n_bound),
                  seq(b_n_bound*2+block.size+1, TT+1-2*b_n_bound, block.size),
                  seq(TT+1-b_n_bound, TT+1,  b_n_bound));
    }else if( !is.null(block.size) && is.null(blocks)){
      # blocks <- seq(1, TT + 1, block.size);
      b_n_bound = 2*block.size  #block size for boundary
      blocks <- c(seq(1, b_n_bound*2+1 , b_n_bound),
                  seq(b_n_bound*2+block.size+1, TT+1-2*b_n_bound, block.size),
                  seq(TT+1-b_n_bound, TT+1,  b_n_bound));

    }else if(!is.null(block.size) && !is.null(blocks)){
      #check if the block.size and blocks match
      n.new <- length(blocks) - 1;
      blocks.size.check <- sapply(c(1:n.new), function(jjj) blocks[jjj + 1] - blocks[jjj]  );
      if( sum(blocks.size.check[1: (length(blocks.size.check) - 1 )] != block.size ) > 0 ){
        stop("Error: The block.size and blocks can't match!")
      }
    }

    if(blocks[length(blocks)] < TT+1){
      blocks <- c(blocks[-length(blocks)], TT + 1)
    }

    n.new <- length(blocks) - 1;
    blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj + 1] - blocks[jjj]  );
    print("block sizes")
    print(blocks.size)
    if(is.null(block.size)){
        block.size <- median(blocks.size)
    }
    print("block size")
    print(block.size)
    #sample the cv index for cross-validation
    #bbb <- floor(n.new/5);
    #aaa <- sample(1:5, 1);
    bbb <- floor(n.new/4);
    aaa <- 4;
    cv.index <- seq(aaa,n.new,floor(n.new/bbb));
    cv.l <- length(cv.index);
    print("cv.index:"); print(cv.index)

    if(method == 'Constant'){
      p.y <- length(data_y[1, ]);
      data_x <- matrix(1, nrow = TT)
      p.x <- 1
      ############# Tuning parameter ################
      if(is.null(lambda.1.cv)){
        lambda.1.max <- lambda_warm_up_lm(data_y, data_x, blocks, cv.index)$lambda_1_max
        if(blocks[2] <= (p.x + p.y) ){
          epsilon <-  10^(-3)
        }
        if(blocks[2] >= (p.x + p.y) ){
          epsilon <-  10^(-4)
        }
        nlam <- 10
        lambda.1.min <-  lambda.1.max*epsilon
        delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
        lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
      }

      if(is.null(lambda.2.cv)){
        lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
      }

      ####################################################
      ########## first step       #######################
      ####################################################
      temp.first <- lm.first.step.blocks(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration, tol = tol, blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
      cp.first <- temp.first$pts.list;
      cl.number <- length(cp.first);
      beta.est <- temp.first$beta.full
      temp.first.all<- temp.first
      jumps.l2 <- temp.first$jumps.l2

      ####################################################
      ########## second step       #######################
      ####################################################
      if(length(cp.first) > 0){
        temp.second<- lm.second.step.search(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
        temp.second.all<- temp.second
        cp.final<- temp.second$cp.final;
        beta.hat.list <- temp.second$beta.hat.list
      }else{
        cp.final <- c()
        print('no change points!')
        beta.hat.list <- beta.est[[floor(n.new/2)]]
      }

      if(refit){
        print('refit the model!')
        cp.full <- c(1, cp.final, TT+1)
        m <- length(cp.final) + 1
        temp_beta <-  list(m)
        for(i in 1:m){
          # data_x_temp <- matrix(1, nrow = TT)
          data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
          temp_beta[[i]] <- apply(data_y_temp, 2, mean)
        }

        beta.final.temp <- matrix(c(do.call("cbind", temp_beta)), nrow = 1)
        lambda.val.best <- BIC.threshold(method, beta.final.temp, p.y,
                                         length(c(cp.final, TT+1)), c(cp.final, TT+1), data_y , b_n = block.size)

        for(j in 1:length(c(cp.final, TT+1))){
          temp_beta[[j]][abs(temp_beta[[j]]) < lambda.val.best[j]] <- 0
        }
        beta.final <- temp_beta

      }else{
        print('no refit!')
        # BIC threholding step!!!
        beta.final.temp <- matrix(c(do.call("cbind", beta.hat.list)), nrow= 1)
        lambda.val.best <- BIC.threshold(method, beta.final.temp, p.y, length(cp.final) + 1, c(cp.final, TT+1),
                                         data_y,  b_n = block.size, nlam = 20)
        # print(lambda.val.best)
        temp <- beta.hat.list
        for(j in 1:(length(cp.final) + 1) ){
          temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
        }
        beta.final <- temp

      }



      return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
                  beta.est = beta.est, beta.final = beta.final, jumps = jumps.l2))
      # return(list(cp.first = cp.first, cp.final = cp.final, beta.hat.list = beta.hat.list, beta.est = beta.est ))



    }
    if(method == 'GGM'){
      p.y <- length(data_y[1, ]);
      data_x <- data_y
      p.x <- p.y
      ############# Tuning parameter ################
      if(is.null(lambda.1.cv)){
        lambda.1.max <- lambda_warm_up_lm(data_y, data_x, blocks, cv.index)$lambda_1_max
        if(blocks[2] <= (p.x + p.y) ){
          epsilon <-  10^(-3)
        }
        if(blocks[2] >= (p.x + p.y) ){
          epsilon <-  10^(-4)
        }
        nlam <- 10
        lambda.1.min <-  lambda.1.max*epsilon
        delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
        lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
      }

      if(is.null(lambda.2.cv)){
        lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
      }

      ####################################################
      ########## first step       #######################
      ####################################################
      temp.first <- ggm.first.step.blocks(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration, tol = tol, blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
      cp.first <- temp.first$pts.list;
      cl.number <- length(cp.first);
      beta.est <- temp.first$beta.full
      temp.first.all<- temp.first
      jumps.l2 <- temp.first$jumps.l2


      ####################################################
      ########## second step       #######################
      ####################################################
      if(length(cp.first) > 0){
        temp.second<- ggm.second.step.search(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
        temp.second.all<- temp.second
        cp.final<- temp.second$cp.final;
        beta.hat.list <- temp.second$beta.hat.list
      }else{
        cp.final <- c()
        print('no change points!')
        beta.hat.list <- beta.est[[floor(n.new/2)]]

      }

      if(refit){
        print('refit the model!')
        temp <- beta.hat.list
        cp.full <- c(1, cp.final, TT+1)
        for(i in 1:(length(cp.final) + 1) ){
          data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
          tmp_coef <- c()
          for(j.1 in 1:p.y){
            data_x_temp <- data_y_temp[, -j.1]
            cvfit = cv.glmnet(data_x_temp, data_y_temp[, j.1], intercept = FALSE)
            opt_lambda  <- cvfit$lambda.1se  # Optimal Lambda
            # fit = glmnet(data_x_temp, data_y_temp, intercept = FALSE, lambda = opt_lambda)
            tmp_coef <- rbind(tmp_coef, coef(cvfit)[-1])

          }
          rownames(tmp_coef) <- NULL
          temp[[i]] <- tmp_coef


        }
        beta.final <- temp

      }else{
        print('no refit!')
        # BIC threholding step!!!
        beta.final.temp <- as.matrix(do.call("cbind", beta.hat.list))
        lambda.val.best <- BIC.threshold.ggm(beta.final.temp, p.x-1, length(cp.final) + 1, c(cp.final, TT+1),
                                             data_y, b_n = block.size, nlam = 50 )$lambda.val.best
        temp <- beta.hat.list
        for(j in 1:(length(cp.final) + 1)){
          temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
        }
        beta.final <- temp

      }

      m <- length(cp.final) + 1
      beta.full.final <- vector("list", m);
      for(i in 1:m){
        beta.full.final[[i]] <- matrix(0, p.y, p.y)
        beta.full.final[[i]][1, 2:p.y] <- beta.final[[i]][1, 1:(p.y-1) ]
        beta.full.final[[i]][p.y,  1:(p.y-1)] <- beta.final[[i]][p.y, 1:(p.y-1)]
        for(j in 2:(p.y-1)){
          beta.full.final[[i]][j, (j+1):p.y] <- beta.final[[i]][j, (j):(p.y-1) ]
          beta.full.final[[i]][j, 1:(j-1)] <- beta.final[[i]][j, 1:(j-1) ]
        }
      }

      # return(list(cp.first = cp.first, cp.final = cp.final, beta.hat.list = beta.hat.list ))
      return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
                  beta.est = beta.est, beta.final = beta.final, beta.full.final= beta.full.final, jumps = jumps.l2))
    }
    if (method == 'MvLR' |  method =="MLR"){
      if (is.null(data_x)){
        stop("Empty predictor data!")
      }
      p.y <- length(data_y[1, ]); p.x <- length(data_x[1, ]);

      ############# Tuning parameter ################
      if(is.null(lambda.1.cv)){
        lambda.1.max <- lambda_warm_up_lm(data_y, data_x, blocks, cv.index)$lambda_1_max
        if(blocks[2] <= (p.x + p.y) ){
          epsilon <-  10^(-3)
        }
        if(blocks[2] >= (p.x + p.y) ){
          epsilon <-  10^(-4)
        }
        nlam <- 10
        lambda.1.min <-  lambda.1.max*epsilon
        delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
        lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
      }

      if(is.null(lambda.2.cv)){
        lambda.2.cv <-  c(10*sqrt( (log(p.x) + log(p.y)  )/TT), 1*sqrt((log(p.x) + log(p.y)  )/TT), 0.10*sqrt((log(p.x) + log(p.y)  )/TT))
      }

      ####################################################
      ########## first step       #######################
      ####################################################
      temp.first <- lm.first.step.blocks(data_y, data_x, lambda.1.cv, lambda.2.cv, max.iteration = max.iteration, tol = tol, blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
      cp.first <- temp.first$pts.list;
      cl.number <- length(cp.first);
      beta.est <- temp.first$beta.full
      temp.first.all<- temp.first
      jumps.l2 <- temp.first$jumps.l2
      ####################################################
      ########## second step       #######################
      ####################################################
      if(length(cp.first) > 0){
        temp.second<- lm.second.step.search(data_y, data_x, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
        temp.second.all<- temp.second
        cp.final<- temp.second$cp.final;
        beta.hat.list <- temp.second$beta.hat.list
      }else{
        cp.final <- c()
        print('no change points!')
        beta.hat.list <- beta.est[[floor(n.new/2)]]
      }

      #final estimates ( BIC thresholding or refitting the model)
      if(refit){
        print('refit the model!')
        temp <- beta.hat.list
        cp.full <- c(1, cp.final, TT+1)
        for(i in 1:(length(cp.final) + 1) ){
          data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
          data_x_temp <- as.matrix(data_x[cp.full[i]: (cp.full[i+1]-1), ])
          if(method == 'MvLR'){
            cvfit = cv.glmnet(data_x_temp, data_y_temp, intercept = FALSE,
                              family = "mgaussian")
            opt_lambda  <- cvfit$lambda.1se  # Optimal Lambda
            tmp_coeffs <- coef(cvfit, s = "lambda.1se")
            tmp_coeffs <- unlist(tmp_coeffs, use.names = FALSE)
            tmp_coef <- tmp_coeffs[[1]][-1]
            if(length(tmp_coeffs) >= 2){
              for(j in 2: length(tmp_coeffs)){
                tmp_coef <- rbind(tmp_coef, tmp_coeffs[[j]][-1])
              }
            }

            rownames(tmp_coef) <- NULL
            temp[[i]] <- tmp_coef
          }else{
            cvfit = cv.glmnet(data_x_temp, data_y_temp, intercept = FALSE)
            opt_lambda  <- cvfit$lambda.1se  # Optimal Lambda
            # fit = glmnet(data_x_temp, data_y_temp, intercept = FALSE, lambda = opt_lambda)
            temp[[i]] <- as.matrix(coef(cvfit)[-1])
          }


        }
        beta.final <- temp

      }else{
        print('no refit!')
        # BIC threholding step!!!
        if(method == 'MLR'){
          beta.final.temp <- matrix(c(do.call("cbind", beta.hat.list)), nrow= 1)
          lambda.val.best <- BIC.threshold(method, beta.final.temp, p.x, length(cp.final) + 1, c(cp.final, TT+1),
                                           data_y, data_x,  b_n = block.size, nlam = 20)

        }else{
          print(beta.hat.list)
          beta.final.temp <- as.matrix(do.call("cbind", beta.hat.list))
          print(dim(beta.final.temp))
          lambda.val.best <- BIC.threshold(method, beta.final.temp, p.x, length(cp.final) + 1, c(cp.final, TT+1),
                                           data_y, data_x,  b_n = block.size, nlam = 20)

        }

        # print(lambda.val.best)
        temp <- beta.hat.list
        for(j in 1:(length(cp.final) + 1) ){
          temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
        }
        beta.final <- temp
      }


      HBIC.full <- c()
      brk.temp <- cp.final
      brk.full <- c(1, brk.temp, TT+1)
      m.hat <- length(brk.full) - 1


      if(method == 'MLR'){
          k <- p.x
          beta.final.mat <- matrix(c(do.call("cbind", beta.final)), nrow = 1)
          HBIC.res <- c()
          residual.full <- c()
          for(i in 1:m.hat){
              beta.temp <- beta.final.mat[((i-1)*k+1):(i*k)]
              data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
              data_x.temp <- data_x[brk.full[i]:(brk.full[i+1]-1), ]
              data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, ncol = 1)
              residual.temp <- data_y.temp - data_y.est
              residual.full <- rbind(residual.full, residual.temp)
          }

          phi =  t(do.call("cbind", beta.final))
          residual <- t(residual.full)
          p.y <- length(phi[, 1]);
          p.x <- length(phi[1, ]);
          T.new <- length(residual[1, ]);
          print("p.x"); print(p.x);
          print("p.y"); print(p.y);
          print("T.new"); print(T.new);
          # count : non-zero coefficient
          count <- 0;
          for (i in 1:p.y){
              for (j in 1:p.x){
                  if(phi[i,j] != 0){
                      count <- count + 1;
                  }
              }
          }

          sigma.hat <- sum(residual^2)/T.new

          log.det <- log(sigma.hat);
          p.y <- 1
          HBIC.res <- log.det*TT + 2*optimal.gamma.val*log(p.x)*count

      }

      HBIC.full <- sum(HBIC.res)


      print("HBIC.full:")
      print(HBIC.full)



      return(list(cp.first = cp.first,  cp.final = cp.final,  beta.hat.list = beta.hat.list,
                  beta.est = beta.est, beta.final = beta.final, jumps = jumps.l2, HBIC.full = HBIC.full))
    }
    if(method == 'VAR'){
      p <- length(data_y[1,]);

      ############# Tuning parameter ################
      if(is.null(lambda.1.cv)){
        lambda.1.max <- lambda_warm_up_var(data_y, q, blocks, cv.index)$lambda_1_max
        if(blocks[2] <= 2*p ){
          epsilon <-  10^(-3)
        }
        if(blocks[2] >= 2*p ){
          epsilon <-  10^(-4)
        }
        nlam <- 10
        lambda.1.min <-  lambda.1.max*epsilon
        delata.lam <- (log(lambda.1.max)-log(lambda.1.min))/(nlam -1)
        lambda.1.cv <-  sapply(1:(nlam), function(jjj) lambda.1.min*exp(delata.lam*(nlam-jjj)))
      }

      if(is.null(lambda.2.cv)){
        lambda.2.cv <-  c(10*sqrt(log(p)/TT),1*sqrt(log(p)/TT),0.10*sqrt(log(p)/TT))
      }

      ####################################################
      ########## first step       #######################
      ####################################################
      temp.first <- var.first.step.blocks(data_y, lambda.1.cv, lambda.2.cv, q = q,  max.iteration = max.iteration, tol = tol,
                                          blocks , cv.index, HBIC = HBIC, gamma.val = gamma.val)
      cp.first <- temp.first$pts.list;
      cl.number <- length(cp.first);
      beta.est <- temp.first$phi.full
      temp.first.all<- temp.first
      jumps.l2 <- temp.first$jumps.l2
      ####################################################
      ########## second step       #######################
      ####################################################
      if(length(cp.first) > 0){
        temp.second<- var.second.step.search(data_y, q = q, max.iteration = max.iteration, tol = tol, cp.first, beta.est, blocks)
        temp.second.all<- temp.second
        cp.final<- temp.second$cp.final;
        beta.hat.list <- temp.second$phi.hat.list
      }else{
        cp.final <- c()
        print('no change points!')
        beta.hat.list <- beta.est[[floor(n.new/2)]]
      }


      #final estimates ( BIC thresholding or refitting the model)
      if(refit){
        print('refit the model!')
        temp <- beta.hat.list
        cp.full <- c(1, cp.final, TT+1)
        for(i in 1:(length(cp.final) + 1) ){
          data_y_temp <- as.matrix(data_y[cp.full[i]: (cp.full[i+1]-1), ])
          # data_x_temp <- as.matrix(data_x[cp.full[i]: (cp.full[i+1]-1), ])
          fit <- fitVAR(data_y_temp, p = 2)
          temp[[i]] <- do.call(cbind, fit$A)
        }
        beta.final <- temp

      }else{
        print('no refit!')
        # BIC threholding step!!!
        beta.final.temp <- as.matrix(do.call("cbind", beta.hat.list))
        print(dim(beta.final.temp))
        lambda.val.best <- BIC.threshold(method, beta.final.temp, p.x*q, length(cp.final) + 1, c(cp.final, TT+1),
                                         data_y, data_x,  b_n = block.size, nlam = 20)

        # print(lambda.val.best)
        temp <- beta.hat.list
        for(j in 1:(length(cp.final) + 1) ){
          temp[[j]][abs(temp[[j]]) < lambda.val.best[j]] <- 0
        }
        beta.final <- temp
      }



      return(list(cp.first = cp.first, beta.est = beta.est, cp.final = cp.final,
                  beta.hat.list = beta.hat.list, jumps.l2 = jumps.l2, beta.final = beta.final ))
    }


  }



}


#' Threshold block fused lasso step for linear regression model.
#'
#' @description Perform the block fused lasso with thresholding to detect candidate break points.
#'
#' @param data_y input data matrix Y, with each column representing the time series component
#' @param data_x input data matrix X, with each column representing the time series component
#' @param lambda1 tuning parmaeter lambda_1 for fused lasso
#' @param lambda2 tuning parmaeter lambda_2 for fused lasso
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cv.index the index of time points for cross-validation
#' @param blocks the blocks
#' @param fixed_index index for linear regression model with only partial compoenents change.
#' @param nonfixed_index index for linear regression model with only partial compoenents change.
#' @param HBIC logical; if TRUE, use high-dimensional BIC, if FALSE, use orginal BIC. Default is FALSE.
#' @param gamma.val hyperparameter for HBIC, if HBIC == TRUE.
#' @return A list object, which contains the followings
#' \describe{
#'   \item{jump.l2}{estimated jump size in L2 norm}
#'   \item{jump.l1}{estimated jump size in L1 norm}
#'   \item{pts.list}{estimated change points in the first step}
#'   \item{beta.full}{estimated parameters in the first step}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
#' @import factoextra
#'
lm.first.step.blocks <- function(data_y, data_x, lambda1, lambda2, max.iteration = max.iteration, tol = tol,
                                 blocks, cv.index, fixed_index = NULL, nonfixed_index = NULL,
                                 HBIC = FALSE, gamma.val = NULL){

  #create the tuning parmaeter combination of lambda1 and lambda2
  lambda.full <- expand.grid(lambda1, lambda2)
  kk <- length(lambda.full[, 1]);

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  cv.l <- length(cv.index);
  cv <- rep(0, kk);
  phi.final <- vector("list", kk);
  phi.final.2 <- vector("list", kk);

  data.y.temp <- data_y; data.x.temp <- data_x
  TT <- length(data.y.temp[,1]);
  p.y <- length(data.y.temp[1,]); p.x.all <- length(data.x.temp[1, ]);
  p.x <- p.x.all - length(fixed_index);

  flag.full <- rep(0, kk);

  for (i in 1:kk) {
    print(i)
    print("##################lambda1:")
    print(lambda.full[i,1])
    print("##################lambda2:")
    print(lambda.full[i,2])
    if(!is.null(fixed_index)){
      print("Some coefficients are constant!")
      if ( i == 1){
        test <- lm_partial_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1],lambda.full[i,2], max_iteration = max.iteration, tol = tol,
        initial_phi =  0.0+matrix(0.0,p.y,p.x*n.new), initial_phi_2 =  0.0+matrix(0.0,p.y,(p.x.all-p.x)), blocks = blocks, cv.index, fixed_index, nonfixed_index)
        flag.full[i] <- test$flag;
      }
      if ( i > 1 ){
        initial.phi <- phi.final[[i-1]]
        initial.phi.2 <- phi.final.2[[i-1]]
        if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
        test <- lm_partial_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1], lambda.full[i,2], max_iteration = max.iteration, tol = tol,
                                           initial_phi =  initial.phi, initial_phi_2 =  initial.phi.2, blocks = blocks, cv.index, fixed_index, nonfixed_index)
        flag.full[i] <- test$flag;
      }
      #print(test$phi.hat)
      #print(test$phi.hat.2)
      phi.hat.full <- test$phi.hat;
      phi.final[[i]] <- phi.hat.full;
      phi.hat.full.2 <- test$phi.hat.2;
      phi.final.2[[i]] <- phi.hat.full.2;

      phi.full.all <- vector("list",n.new);
      forecast <- matrix(0,p.y,TT);
      forecast.new <- matrix(0,p.y,cv.l);

      phi.hat <- phi.hat.full;
      phi.full.all[[1]] <- matrix(phi.hat[,(1):(p.x)],ncol = p.x);
      for(i.1 in 2:n.new){
        phi.full.all[[i.1]] <- matrix(phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        forecast[,(blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x[,nonfixed_index]), phi.full.all[[i.1-1]], blocks[i.1], p.x, p.y, blocks[i.1+1]-blocks[i.1]);
      }

      forecast.new <- matrix(0,p.y,cv.l);
      for(j in (1):cv.l){
        forecast.new[,j] <- pred(t(data_x[,nonfixed_index]),phi.full.all[[(cv.index[j])]],blocks[cv.index[j]+1]-1,p.x, p.y)
        forecast.new[,j] <- forecast.new[,j] + phi.hat.full.2%*%as.matrix(data_x[blocks[cv.index[j]+1]-1,fixed_index])
      }

      temp.index <- rep(0,cv.l);
      for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1]-1;}
      cv[i] <- (1/(p.y*cv.l))*sum( (forecast.new - t(data_y[temp.index,])  )^2 );

      print("============cv-result=======================")
      print(cv[i])
      print("====================================")


    }
    if(is.null(fixed_index)){
      print("All coefficients are piecewise constant!")
      if ( i == 1){
        test <- lm_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1], lambda.full[i,2], max_iteration = max.iteration, tol = tol, initial_phi =  0.0+matrix(0.0,p.y,p.x*n.new), blocks = blocks, cv.index)
        flag.full[i] <- test$flag;
      }
      if ( i > 1 ){
        initial.phi <- phi.final[[i-1]]
        if(max(abs(phi.final[[(i-1)]])) > 10^3 | flag.full[i-1] == 1  ){
          print("bad initial!")
          initial.phi <- 0*phi.final[[(i-1)]];
        }
        test <- lm_break_fit_block(data.y.temp,data.x.temp, lambda.full[i,1], lambda.full[i,2], max_iteration = max.iteration, tol = tol, initial_phi =  initial.phi, blocks = blocks, cv.index)
        flag.full[i] <- test$flag;
      }
      #print(test$phi.hat)
      phi.hat.full <- test$phi.hat;
      phi.final[[i]] <- phi.hat.full;


      #forecast the time series based on the estimated matrix Phi (beta for the linear regression model)
      #and compute the forecast error
      phi.full.all <- vector("list", n.new);
      forecast <- matrix(0, p.y, TT);
      forecast.new <- matrix(0, p.y, cv.l);

      phi.hat <- phi.hat.full;
      phi.full.all[[1]] <- matrix(phi.hat[,(1):(p.x)], ncol = p.x);
      forecast[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data_x), phi.full.all[[1]], blocks[1], p.x, p.y, blocks[2]-blocks[1]);
      for(i.1 in 2:n.new){
        phi.full.all[[i.1]] <- matrix(phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        forecast[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x), phi.full.all[[i.1]], blocks[i.1], p.x, p.y, blocks[i.1+1]-blocks[i.1]);
      }
      forecast.new <- matrix(0, p.y, cv.l);
      for(j in (1):cv.l){
        forecast.new[, j] <- pred(t(data_x), phi.full.all[[(cv.index[j])]], blocks[cv.index[j]+1]-1, p.x, p.y)
      }
      temp.index <- rep(0, cv.l);
      for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1]-1;}
      # print(phi.full.all)
      # print(forecast.new)
      # print(data_y[temp.index,])
      cv[i] <- (1/(p.y*cv.l))*sum( (forecast.new - t(data_y[temp.index, ])  )^2 );
      # cv[i] <- (1/(p.y*cv.l))*sum( ( (forecast.new - t(data_y[temp.index,]))/t(data_y[temp.index,])  )^2);

      print("============cv-result=======================")
      print(cv[i])
      print("====================================")
    }

  }



  lll <- min(which(cv==min(cv)));

  mspe.plot(cv, c(1:kk))
  abline(v = seq(length(lambda1), length(lambda1)*(length(lambda2)-1), length.out =length(lambda2)-1)+0.5)
  abline(v= lll, col="red")

  print("All cv:")
  print(cv)

  phi.hat.full <- phi.final[[lll]];
  beta.fixed.full <- phi.final.2[[lll]];

  #jumps.sq is the L2 norm square
  #jumps.l1 is the L1 norm
  # again, note that here the phi.hat.full is the estimated theta in the paper
  jumps.l2 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l2[i] <- (sum((phi.hat.full[,((i-1)*p.x+1):(i*p.x)] )^2 ));
  }
  jumps.l1 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l1[i] <- (sum(abs(phi.hat.full[,((i-1)*p.x+1):(i*p.x)] ) ));
  }

  # print(jumps.l2)
  # print(jumps.l1)
  # plot(jumps.l2,type = 'o', main= 'l2 norm')
  # plot(jumps.l1,type = 'o', main= 'l1 norm')

  #ignore the large jump at the boundary!!!!!!!!!!!!!!!!!!!!
  # print('use l2 norm!')
  jumps <- jumps.l2
  # print('use l1 norm!')
  # jumps <- jumps.l1
  ignore_num <- min(5 , round(200/mean(blocks.size)))
  jumps[1:ignore_num]  <- 0
  jumps[(length(jumps)-(ignore_num-1)):length(jumps)]  <- 0


  plot(jumps, type = 'o', main= 'l2 norm')
  ##################################################################
  # use BIC to determine the k-means
  BIC.diff <- 1
  BIC.old <- 10^8
  pts.sel <- c()
  loc.block.full <- c()
  while(BIC.diff > 0 & length(unique(jumps)) > 1 ){
    pts.sel.old <- pts.sel
    #use kmeans to hard threshold the jumps
    if( length(unique(jumps)) > 2 ){
      # print("consider 2 clusters for fit.2")
      clus.2 <- kmeans(jumps, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss;
      # print(fit.2);
      if(fit.2 < 0.20){
        print("no significant jumps!!")
        pts.sel <- c(pts.sel);
      }
      if( fit.2 >= 0.20 ){
        loc <- clus.2$cluster;
        if( clus.2$centers[1] > clus.2$centers[2]  ){
          loc.block <- which(loc==1);
        }
        if( clus.2$centers[1] < clus.2$centers[2]  ){
          loc.block <- which(loc==2);
        }
        pts.sel <- c(pts.sel, blocks[loc.block]);
        loc.block.full <- c(loc.block.full, loc.block)
      }
    }
    if( length(unique(jumps)) <= 2 ){
      if(length(unique(jumps)) == 2){
        loc.block <- which.max(jumps)
        pts.sel <- c(pts.sel, blocks[loc.block])
        loc.block.full <- c(loc.block.full, loc.block)
      }else{
        pts.sel <- c(pts.sel);
      }
    }

    # print("pts.sel:"); print(sort(pts.sel))

    phi.hat.full.new <- phi.hat.full
    for(i in 2:n.new){
      if(!(i %in% loc.block.full)){
        phi.hat.full.new[, ((i-1)*p.x+1):(i*p.x)] <- matrix(0, p.y, p.x)
      }
    }


    phi.full.all.new <- vector("list", n.new);
    phi.full.all.new.temp <- vector("list", n.new);
    forecast.all.new <- matrix(0, p.y, TT);

    phi.hat.new <- phi.hat.full.new;
    phi.full.all.new[[1]] <- matrix(phi.hat.new[, (1):(p.x)], ncol = p.x);
    phi.full.all.new.temp[[1]] <- matrix(phi.hat.new[, (1):(p.x)], ncol = p.x);

    if(!is.null(fixed_index)){
      forecast.all.new[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data_x[, nonfixed_index]), phi.full.all.new[[1]], blocks[1], p.x, p.y, blocks[2] - blocks[1]);
      forecast.all.new[, (blocks[1]):(blocks[2]-1)] <- forecast.all.new[, (blocks[1]):(blocks[2]-1)] + beta.fixed.full %*%t(as.matrix(data_x[(blocks[1]):(blocks[2]-1), fixed_index]))

      for(i.1 in 2:n.new){
        phi.full.all.new[[i.1]] <- matrix(phi.full.all.new[[i.1-1]] + phi.hat.new[,((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x[, nonfixed_index]), phi.full.all.new[[i.1]], blocks[i.1], p.x, p.y, blocks[i.1+1] - blocks[i.1]);
        forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] + beta.fixed.full %*%t(as.matrix(data_x[(blocks[i.1]):(blocks[i.1+1]-1), fixed_index]))
      }
    }
    if(is.null(fixed_index)){

      forecast.all.new[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data_x), phi.full.all.new[[1]],
                                                                  blocks[1], p.x, p.y, blocks[2] - blocks[1]);
      for(i.1 in 2:n.new){
        #phi.full.all.new.temp keeps adding
        phi.full.all.new.temp[[i.1]] <- matrix(phi.full.all.new.temp[[i.1-1]] + phi.hat.full[, ((i.1-1)*p.x+1):(i.1*p.x)], ncol = p.x);
        if((i.1 %in% loc.block.full)){
          phi.full.all.new[[i.1]] <- phi.full.all.new.temp[[i.1]]
        }
        if(!(i.1 %in% loc.block.full)){
          phi.full.all.new[[i.1]] <- phi.full.all.new[[i.1-1]]
        }
        forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data_x), phi.full.all.new[[i.1]],
                                                                          blocks[i.1], p.x, p.y, blocks[i.1+1] - blocks[i.1]);
      }
    }
    residual <- t(data.y.temp[(1:TT), ]) - forecast.all.new;
    # print(phi.full.all.new)
    # print(phi.full.all.new.temp)
    # print(residual)

    if(HBIC == TRUE){
      print("Use HBIC!")
      # print('gamma value:')
      # print(gamma.val)
      if(is.null(gamma.val)){
        BIC.new <- BIC(residual, phi = phi.hat.full.new)$HBIC
      }else{
        BIC.new <- BIC(residual, phi = phi.hat.full.new, gamma.val = gamma.val)$HBIC
      }

    }else{
      print("Use BIC!")
      BIC.new <- BIC(residual, phi = phi.hat.full.new )$BIC
    }
    # print("BIC.new:"); print(BIC.new)
    BIC.diff <- BIC.old - BIC.new
    # print("BIC.diff:");print(BIC.diff)
    BIC.old <- BIC.new
    if(BIC.diff <= 0){
      pts.sel <- sort(pts.sel.old)
      break
    }
    jumps[loc.block] <- 0
    plot(jumps, type = 'o', main= 'l2 norm')
  }

  # print(pts.sel)


  if ( length(pts.sel) == 0){cp.final <- c(); pts.list <-  vector("list", 0);}
  if ( length(pts.sel) > 0){
    cp.final <- pts.sel;
    # REMOVE BOUNDARY POINTS and other cleaning
    # cp.final <- cp.final[which(cp.final > 3*mean(blocks.size))];
    # cp.final <- cp.final[which(cp.final < (TT-3*mean(blocks.size)))];
    cp.final <- cp.final[which(cp.final > sum(blocks.size[1:3]))];
    cp.final <- cp.final[which(cp.final < (TT-sum(blocks.size[(length(blocks.size)-2):length(blocks.size)])))];
    cp.final <- sort(cp.final);
    # print("cp.final:"); print(cp.final)

    if(length(cp.final) >=2){
      gap.temp <- sapply(2:length(cp.final), function(jjj) cp.final[jjj]-cp.final[jjj-1])
      # print("gap.temp:");print(gap.temp)
    }


    #if there are multipler change points
    # use gap stat and the k-means clustering
    if(length(cp.final) > 5){
      if(length(unique(gap.temp)) > 1  ){

        print(fviz_nbclust(matrix(cp.final, length(cp.final), 1), kmeans, nstart = 25,  method = "gap_stat",
                           k.max = min(50, length(cp.final)-1), nboot = 100)+
                labs(subtitle = "Gap statistic method"))
        cl <- fviz_nbclust(matrix(cp.final, length(cp.final), 1), kmeans, nstart = 25,  method = "gap_stat",
                           k.max = min(50, length(cp.final)-1), nboot = 100)+
          labs(subtitle = "Gap statistic method")
        cl.data <- cl$data;
        gap <- cl.data$gap;
        # se <- cl.data$SE.sim;
        # i.cl <- 0;
        print("choose the maximum gap stat")
        cl.number <- which.max(gap)
        gap.order <- order(gap, decreasing = TRUE)
        print("check if the diameter of clusters is small enough")

        if(median(blocks.size) <= sqrt(TT)/4){
          cnst <- 2*4+1
        }else if(median(blocks.size) <= 2*sqrt(TT)/4){
          cnst <- 2*3+1
        }else{
          cnst <- 2*2+1
        }

        flag = TRUE
        idx <- 1
        while(flag & idx <= length(gap.order)){
          cl.number <- gap.order[idx]
          print("cl.number:")
          print(cl.number)
          if(cl.number > 1){
            cluster.pos <- sort(order(gap.temp, decreasing = TRUE)[1:(cl.number-1)])
            cluster.pos <- c(0, cluster.pos, length(cp.final))
          }else{
            cluster.pos <- c(0, length(cp.final))
          }

          wide <- 0
          for (i in c(1:cl.number)) {
            pts.i <- cp.final[(cluster.pos[i]+1): cluster.pos[i+1]]
            print(paste0("Cluster: ", i))
            print(pts.i)
            block_avg <- mean(blocks.size[match(pts.i, blocks)])
            # if ( max(pts.i) - min(pts.i) > 5*mean(blocks.size) ){
            if ( max(pts.i) - min(pts.i) > cnst*block_avg ){
              print("Too wide, need seperate!")
              wide <- 1
              break
            }
          }
          # print(idx)
          if (wide){
            idx <- idx + 1
          }else{
            idx <- idx
            flag <- FALSE
          }
        }
        if(flag){
          print("they are seperate change points! not use gap stat!")
          cl.number <- length(gap.temp) + 1
        }

        # print(cl.number)
        print("check if distance between two clusters are large enough")
        if(median(blocks.size) <= sqrt(TT)/4){
          cnst2 <- 7
        }else if(median(blocks.size) <= sqrt(TT)/2){
          cnst2 <- 5
        }else{
          cnst2 <- 3
        }

        if(cl.number > 1){
          cl.number <- sum(sort(gap.temp, decreasing = TRUE)[1:(cl.number - 1)] > cnst2*median(blocks.size)) + 1
        }
      }else if(unique(gap.temp) == median(blocks.size) ){
          print("one single cluster")
          cl.number <- 1
      }else{
        print("they are seperate change points! not use gap stat!")
        cl.number <- length(gap.temp) + 1
      }
      print('number of clusters:')
      print(cl.number)

      if(cl.number > 1){
        cluster.pos <- sort(order(gap.temp, decreasing = TRUE)[1:(cl.number-1)])
        cluster.pos <- c(0, cluster.pos, length(cp.final))
      }else{
        cluster.pos <- c(0, length(cp.final))
      }

      pts.list <-  vector("list", cl.number);
      for (i in c(1:cl.number)) {
        pts.i <- cp.final[(cluster.pos[i]+1): cluster.pos[i+1]]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }

    #   # three times in case the kmeans gives wrong clusters
    #   cl.final.1 <- kmeans(scale(cp.final), centers = cl.number);
    #   fit.1 <- cl.final.1$betweenss/cl.final.1$totss
    #   print(cl.final.1)
    #   cl.final.2 <- kmeans(scale(cp.final), centers = cl.number);
    #   fit.2 <- cl.final.2$betweenss/cl.final.2$totss
    #   print(cl.final.2)
    #   cl.final.3 <- kmeans(scale(cp.final), centers = cl.number);
    #   fit.3 <- cl.final.3$betweenss/cl.final.3$totss
    #   print(cl.final.3)
    #   if(fit.1 >= max(fit.2, fit.3) ){
    #     cl.final <- cl.final.1
    #   }else if(fit.2 >= fit.3){
    #     cl.final <- cl.final.2
    #   }else{
    #     cl.final <- cl.final.3
    #   }
    #   pts.list <-  vector("list",cl.number);
    #   loc.new <- cl.final$cluster;
    #   print(loc.new)
    #   cl.reorder = c(1:cl.number)[order(cl.final$centers)]
    #   for (i in c(1:cl.number)) {
    #     pts.i <- cp.final[which(loc.new==cl.reorder[i])]
    #     print(paste0("Cluster: ", i))
    #     print(pts.i)
    #     pts.list[[i]] <- pts.i
    #   }
    # }

    #!!!!!!!!!!!!!!!!need  ajustment!!!!!!
    if(length(cp.final) <= 5 & length(cp.final) > 1 ){
      print("small number of cp !!!!")
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))

      for (i in 2:length(cp.final)){
        if (cp.final[i]-cp.final[i-1]<= max(3*mean(blocks.size))){
          cl.number <-cl.number-1
          loc.new[i] <-loc.new[i-1]
        }else{
          loc.new[i] <- i
        }
      }

      pts.list <-  vector("list",cl.number);
      #for (i in unique(loc.new)) {
      loc.new.unique <- unique(loc.new)
      for (i in 1:length(loc.new.unique)) {
        pts.i <- cp.final[which(loc.new==loc.new.unique[i])]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }


    if(length(cp.final) == 1 ){
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))
      pts.list <-  vector("list",cl.number);
      for (i in unique(loc.new)) {
        pts.i <- cp.final[which(loc.new==i)]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }

    if(length(cp.final) == 0 ){
      pts.list <-  vector("list", 0);
    }

  }

  #compute the estimated beta
  phi.par.sum <- vector("list",n.new);
  phi.par.sum[[1]] <- phi.hat.full[,1:(p.x)];
  for(i in 2:n.new){
    phi.par.sum[[i]] <- phi.par.sum[[i-1]] + phi.hat.full[,((i-1)*p.x+1):(i*p.x)];
  }

  #plot(blocks[1:n.new],jumps,main = "JUMPS.FULL", type = "o")

  print("First step DONE!!!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  return(list(jumps.l2 = jumps.l2, jumps.l1 = jumps.l1, pts.list = pts.list, beta.full = phi.par.sum))

}



#' Exhaustive search step for linear regression model.
#'
#' @description Perform the exhaustive search to "thin out" redundant break points.
#'
#' @param data_y input data matrix, with each column representing the time series component
#' @param data_x input data matrix, with each column representing the time series component
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cp.first the selected break points after the first step
#' @param beta.est the estiamted parameters by block fused lasso
#' @param blocks the blocks
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{cp.final}{a set of selected break point after the exhaustive search step}
#'   \item{beta.hat.list}{the estimated coefficient matrix for each segmentation}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
lm.second.step.search <- function(data_y, data_x, max.iteration = max.iteration, tol = tol,
                                  cp.first, beta.est, blocks){

  TT <- length(data_y[,1]); p.y <- length(data_y[1,]); p.x <- length(data_x[1,]);

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );


  cp.search <- cp.first;
  cl.number <- length(cp.first)

  cp.list <- vector("list", cl.number + 2);
  cp.list[[1]] <- c(1);
  cp.list[[cl.number+2]] <- c(TT+1);

  cp.index.list <- vector("list", cl.number + 2);
  cp.index.list[[1]] <- c(1);
  cp.index.list[[cl.number+2]] <- c(n.new+1);

  for (i in 1:cl.number) {
    cp.list[[i+1]] <- cp.first[[i]]
    cp.index.list[[i+1]] <- match(cp.first[[i]], blocks)
  }


  cp.search <- rep(0,cl.number);
  cp.list.full <- cp.list

  beta.hat.list <- vector("list", cl.number+1)
  for(i in 1:(cl.number)){
    idx <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2);
    beta.hat.list[[i]] <- beta.est[[idx]]

    if(length(cp.list[[i+1]]) >1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] + 1)  :  (cp.list[[i+1]][length(cp.list[[i+1]])]-1  ) )
    }
    if(length(cp.list[[i+1]]) == 1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] -  (blocks.size[cp.index.list[[i+1]][1] ]) + 1) :  (cp.list[[i+1]][length(cp.list[[i+1]])] +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1   ) )
    }


    #compare the SSE of first num and last num
    num  <- cp.list.full[[i+1]][1]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    beta.hat <- beta.est[[idx.1]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.1 == 1){
      temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
    }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    beta.hat <- beta.est[[idx.2]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.2 == 1){
      temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
    }
    sse1 <- temp.1 + temp.2;


    num  <- cp.list.full[[i+1]][length(cp.list.full[[i+1]])]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    beta.hat <- beta.est[[idx.1]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.1 == 1){
      temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
    }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    beta.hat <- beta.est[[idx.2]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    if(len.2 == 1){
      temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
    }
    sse2 <- temp.1 + temp.2;



    if(sse1 <= sse2){
      sse.full <- 0;
      ii <- 0
      # print(cp.list.full[[i+1]] )
      for(num in cp.list.full[[i+1]]  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        beta.hat <- beta.est[[idx.1]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
        if(len.1 == 1){
          temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );

        }else{
          temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
        }


        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        beta.hat <- beta.est[[idx.2]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x),matrix(beta.hat,ncol = p.x),jjj, p.x, p.y)  )
        if(len.2 == 1){
          temp.2 <- sum( (data_y[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full,0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][min(which(sse.full == min(sse.full)))];

    }
    if(sse1 > sse2){
      sse.full <- 0;
      ii <- 0
      # print(rev(cp.list.full[[i+1]]) )
      for(num in rev(cp.list.full[[i+1]])  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        beta.hat <- beta.est[[idx.1]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
        if(len.1 == 1){
          temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );

        }else{
          temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
        }


        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        beta.hat <- beta.est[[idx.2]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x),matrix(beta.hat,ncol = p.x),jjj, p.x, p.y)  )
        if(len.2 == 1){
          temp.2 <- sum( (data_y[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full,0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][length(cp.list.full[[i+1]]) + 1 - min(which(sse.full == min(sse.full)))];

    }

  }

  print("cp.final:")
  print(cp.search)
  print("Second step DONE!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

  idx <- floor((min(cp.index.list[[cl.number+1+1]]) + max(cp.index.list[[cl.number+1]]))/2);
  beta.hat.list[[cl.number+1]] <- beta.est[[idx]]
  return(list(cp.final = cp.search, beta.hat.list = beta.hat.list))
}




#' Threshold block fused lasso step for gaussian graphical model.
#'
#' @description Perform the block fused lasso with thresholding to detect candidate break points.
#'
#' @param data_y input data matrix Y
#' @param data_x input data matrix X
#' @param lambda1 tuning parmaeter lambda_1 for fused lasso
#' @param lambda2 tuning parmaeter lambda_2 for fused lasso
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cv.index the index of time points for cross-validation
#' @param blocks the blocks
#' @param HBIC logical; if TRUE, use high-dimensional BIC, if FALSE, use orginal BIC. Default is FALSE.
#' @param gamma.val hyperparameter for HBIC, if HBIC == TRUE.
#' @return A list object, which contains the followings
#' \describe{
#'   \item{jump.l2}{estimated jump size in L2 norm}
#'   \item{jump.l1}{estimated jump size in L1 norm}
#'   \item{pts.list}{estimated change points in the first step}
#'   \item{beta.full}{estimated parameters in the first step}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
#' @import factoextra
ggm.first.step.blocks <- function(data_y, data_x, lambda1, lambda2, max.iteration = max.iteration, tol = tol,
                                  blocks, cv.index, HBIC = FALSE, gamma.val = NULL){

  #create the tuning parmaeter combination of lambda1 and lambda2
  lambda.full <- expand.grid(lambda1, lambda2)
  kk <- length(lambda.full[, 1]);

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );
  cv.l <- length(cv.index);
  cv.all <- rep(0, kk);
  phi.final <- vector("list", kk);
  phi.final.2 <- vector("list", kk);
  TT <- length(data_y[, 1]);
  p.y <- length(data_y[1, ]); p.x <- length(data_x[1, ]);
  p.x.temp <- p.x - 1;
  p.y.temp <- 1;

  flag.full <- rep(0, kk);

  for (i in 1:kk){
    print(i)
    print("##################lambda1:")
    print(lambda.full[i, 1])
    print("##################lambda2:")
    print(lambda.full[i, 2])
    if(i == 1){
      test <- ggm_break_fit_block(data_y, data_x, lambda.full[i, 1], lambda.full[i, 2], max_iteration = max.iteration, tol = tol, initial_phi = 0.0+matrix(0.0, p.y.temp, p.x.temp*n.new), blocks = blocks, cv.index)
      phi.hat.full.all <- test$phi.hat
      phi.final[[i]] <- phi.hat.full.all;
    }else{
      initial.phi <- phi.final[[i-1]]
      if(max(abs(phi.final[[(i-1)]])) > 10^3){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- ggm_break_fit_block(data_y, data_x, lambda.full[i, 1], lambda.full[i, 2], max_iteration = max.iteration, tol = tol, initial_phi = initial.phi, blocks = blocks, cv.index)
      phi.hat.full.all <- test$phi.hat
      phi.final[[i]] <- phi.hat.full.all;
    }

    for(j in 1:p.y){
      cv <- rep(0, kk);
      data.y.temp <- as.matrix(data_y[, j])
      data.x.temp <- data_x[, -j]
      #forecast the time series based on the estimated matrix Phi (beta for the linear regression model)
      #and compute the forecast error
      phi.full.all <- vector("list", n.new);
      forecast <- matrix(0, p.y.temp, TT);
      forecast.new <- matrix(0, p.y.temp, cv.l);

      phi.hat.full <- matrix(phi.hat.full.all[j, ], nrow = p.y.temp)
      phi.hat <- phi.hat.full;
      phi.full.all[[1]] <- matrix(phi.hat[, (1):(p.x.temp)], ncol = p.x.temp);
      forecast[, (blocks[1]):(blocks[2]-1)] <- pred.block(t(data.x.temp), phi.full.all[[1]], blocks[1], p.x.temp, p.y.temp, blocks[2]-blocks[1]);
      for(i.1 in 2:n.new){
        phi.full.all[[i.1]] <- matrix(phi.full.all[[i.1-1]] + phi.hat[, ((i.1-1)*p.x.temp + 1):(i.1*p.x.temp)], ncol = p.x.temp);
        forecast[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block(t(data.x.temp), phi.full.all[[i.1]], blocks[i.1], p.x.temp, p.y.temp, blocks[i.1+1]-blocks[i.1]);
      }
      forecast.new <- matrix(0, p.y.temp, cv.l);
      for(j.1 in (1):cv.l){
        forecast.new[, j.1] <- pred(t(data.x.temp), phi.full.all[[(cv.index[j.1])]], blocks[cv.index[j.1]+1]-1, p.x.temp, p.y.temp)
      }
      temp.index <- rep(0, cv.l);
      for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1]-1;}
      cv[i] <- (1/(p.y.temp*cv.l))*sum( (forecast.new - t(data.y.temp[temp.index, ])  )^2 );

      print("============cv-result=======================")
      print(cv[i])
      cv.all[i] <- cv.all[i] + cv[i]
      print("====================================")

    }

  }


  lll <- min(which(cv.all == min(cv.all)));

  mspe.plot(cv.all, c(1:kk))
  abline(v = seq(length(lambda1), length(lambda1)*(length(lambda2)-1), length.out = length(lambda2)-1)+0.5)
  abline(v = lll, col = "red")

  phi.hat.full <- phi.final[[lll]];

  #jumps.sq is the L2 norm square
  #jumps.l1 is the L1 norm
  # again, note that here the phi.hat.full is the estimated theta in the paper
  jumps.l2 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l2[i] <- (sum((phi.hat.full[, ((i-1)*p.x.temp+1):(i*p.x.temp)] )^2 ));
  }
  jumps.l1 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l1[i] <- (sum(abs(phi.hat.full[, ((i-1)*p.x.temp+1):(i*p.x.temp)] ) ));
  }

  # print(jumps.l2)
  # print(jumps.l1)
  # plot(jumps.l2,type = 'o', main= 'l2 norm')
  # plot(jumps.l1,type = 'o', main= 'l1 norm')

  #ignore the large jump at the boundary!!!!!!!!!!!!!!!!!!!!
  # print('use l2 norm!')
  jumps <- jumps.l2
  # print('use l1 norm!')
  # jumps <- jumps.l1
  # ignore_num <- min(3 , round(100/mean(blocks.size)))
  #for eeg data set
  # ignore_num <- max(6 , round(250/min(blocks.size)))
  # for small TT
  if(TT < 250){
    ignore_num <- min(3 , round(100/mean(blocks.size)))
  }else{
    #for eeg data set
    ignore_num <- max(6 , round(250/min(blocks.size)))
  }

  jumps[1:ignore_num]  <- 0
  jumps[(length(jumps)-(ignore_num-1)):length(jumps)]  <- 0


  plot(jumps, type = 'o', main = 'l2 norm')
  ##################################################################
  # use BIC to determine the k-means
  BIC.diff <- 1
  BIC.old <- 10^8
  pts.sel <- c()
  loc.block.full <- c()


  # updated : 09/09/2020: instead of zero, set the jump of selected block to be the max value of jump among those non-selected blocks.
  # while(BIC.diff > 0 & length(unique(jumps)) > 1 ){
  #   pts.sel.old <- pts.sel
  #   #use kmeans to hard threshold the jumps
  #   if( length(unique(jumps)) > 2 ){
  #     # print("consider 2 clusters for fit.2")
  #     clus.2 <- kmeans(jumps, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss;
  #     # print(fit.2);
  #     if(fit.2 < 0.20){
  #       print("no significant jumps!!")
  #       pts.sel <- c(pts.sel);
  #     }
  #     if( fit.2 >= 0.20 ){
  #       loc <- clus.2$cluster;
  #       if( clus.2$centers[1] > clus.2$centers[2]  ){
  #         loc.block <- which(loc==1);
  #       }
  #       if( clus.2$centers[1] < clus.2$centers[2]  ){
  #         loc.block <- which(loc==2);
  #       }
  #       pts.sel <- blocks[loc.block];
  #       loc.block.full <- loc.block
  #     }
  #   }
  #   if( length(unique(jumps)) <= 2 ){
  #     if(length(unique(jumps)) == 2){
  #       loc.block <- which.max(jumps)
  #       pts.sel <- blocks[loc.block]
  #       loc.block.full <- loc.block
  #     }else{
  #       pts.sel <- c(pts.sel);
  #     }
  #   }
  #
  #   # print("pts.sel:"); print(sort(pts.sel))
  #
  #   phi.hat.full.new <- phi.hat.full
  #   for(i in 2:n.new){
  #     if(!(i %in% loc.block.full)){
  #       phi.hat.full.new[, ((i-1)*p.x.temp+1):(i*p.x.temp)] <- matrix(0, p.y, p.x.temp)
  #     }
  #   }
  #
  #   phi.full.all.new <- vector("list", n.new);
  #   phi.full.all.new.temp <- vector("list", n.new);
  #   forecast.all.new <- matrix(0, p.y, TT);
  #
  #   # update: 09/09/2020: use the estimation in middle index (same as second step)
  #   phi.hat.new <- phi.hat.full.new;
  #   phi.full.all.new[[1]] <- matrix(phi.hat.new[, (1):(p.x.temp)], ncol = p.x.temp);
  #   phi.full.all.new.temp[[1]] <- matrix(phi.hat.new[, (1):(p.x.temp)], ncol = p.x.temp);
  #   for(i.1 in 2:n.new){
  #     #phi.full.all.new.temp keeps adding
  #     phi.full.all.new.temp[[i.1]] <- matrix(phi.full.all.new.temp[[i.1-1]] + phi.hat.full[, ((i.1-1)*p.x.temp+1):(i.1*p.x.temp)], ncol = p.x.temp);
  #   }
  #
  #   beta.hat.list.new <- vector("list", length(loc.block.full) + 1) ;
  #   # includes the boundary block, and sorted
  #   loc.block.new <-  c(1, sort(loc.block.full), n.new+1)
  #   idx.1 <- floor((loc.block.new[1] + loc.block.new[2])/2) ;
  #   beta.hat.list.new[[1]] <- phi.full.all.new.temp[[idx.1]]
  #   lb.1 <- blocks[loc.block.new[1]];
  #   ub.1 <- blocks[loc.block.new[2]] - 1;
  #   # print("lb:");print(lb.1)
  #   # print("ub:");print(ub.1)
  #   for(j.1 in 1:p.y){
  #     data.x.temp <- data_x[, -j.1]
  #     forecast.all.new[j.1, lb.1:ub.1] <- pred.block(t(data.x.temp), matrix(beta.hat.list.new[[1]][j.1, ], ncol = p.x.temp),
  #                                                    lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1);
  #   }
  #
  #   # print(sort(loc.block.full))
  #   for(i.1 in 1:length(loc.block.full)){
  #     idx.1 <- floor((loc.block.new[i.1 + 1] + loc.block.new[i.1 + 2])/2);
  #     # print(idx.1)
  #     beta.hat.list.new[[i.1 + 1]] <- phi.full.all.new.temp[[idx.1]]
  #     lb.1 <- blocks[loc.block.new[i.1 + 1]];
  #     ub.1 <- blocks[loc.block.new[i.1 + 2]] - 1;
  #     # print("lb:");print(lb.1)
  #     # print("ub:");print(ub.1)
  #     for(j.1 in 1:p.y){
  #       data.x.temp <- data_x[, -j.1]
  #       forecast.all.new[j.1, lb.1:ub.1] <- pred.block(t(data.x.temp), matrix(beta.hat.list.new[[i.1 + 1]][j.1, ], ncol = p.x.temp),
  #                                                      lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1);
  #     }
  #   }
  #   # print(forecast.all.new[1, ])
  #
  #   residual <- t(data_y[(1:TT), ]) - forecast.all.new;
  #
  #   if(HBIC == TRUE){
  #     print("Use HBIC!")
  #     # print('gamma value:')
  #     # print(gamma.val)
  #     if(is.null(gamma.val)){
  #       BIC.new <- BIC(matrix(residual[1, ], nrow = 1), phi = matrix(phi.hat.full.new[1, ], nrow = 1))$HBIC
  #     }else{
  #       BIC.new <- BIC(matrix(residual[1, ], nrow = 1), phi = matrix(phi.hat.full.new[1, ], nrow = 1), gamma.val = gamma.val)$HBIC
  #     }
  #
  #   }else{
  #     print("Use BIC!")
  #     BIC.new <- BIC(matrix(residual[1, ], nrow = 1), phi = matrix(phi.hat.full.new[1, ], nrow = 1) )$BIC
  #   }
  #   # print("BIC.new:"); print(BIC.new)
  #   BIC.diff <- BIC.old - BIC.new
  #   # print("BIC.diff:"); print(BIC.diff)
  #   BIC.old <- BIC.new
  #   if(BIC.diff <= 0){
  #     pts.sel <- sort(pts.sel.old)
  #     break
  #   }
  #   jumps[loc.block] <- max(jumps[-loc.block])
  #   plot(jumps, type = 'o', main= 'l2 norm')
  # }

  while(BIC.diff > 0 & length(unique(jumps)) > 1 ){
    pts.sel.old <- pts.sel
    #use kmeans to hard threshold the jumps
    if( length(unique(jumps)) > 2 ){
      print("consider 2 clusters for fit.2")
      clus.2 <- kmeans(jumps, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss; print(fit.2);
      if(fit.2 < 0.20){
        print("no significant jumps!!")
        pts.sel <- c(pts.sel);
      }
      if( fit.2 >= 0.20 ){
        loc <- clus.2$cluster;
        if( clus.2$centers[1] > clus.2$centers[2]  ){
          loc.block <- which(loc==1);
        }
        if( clus.2$centers[1] < clus.2$centers[2]  ){
          loc.block <- which(loc==2);
        }
        # pts.sel <- blocks[loc.block];
        # loc.block.full <- loc.block
        pts.sel <- unique(c(pts.sel, blocks[loc.block]));
        loc.block.full <- unique(c(loc.block.full, loc.block))
      }
    }
    if( length(unique(jumps)) <= 2 ){
      if(length(unique(jumps)) == 2){
        loc.block <- which.max(jumps)
        # pts.sel <- blocks[loc.block]
        # loc.block.full <- loc.block
        pts.sel <- unique(c(pts.sel, blocks[loc.block]));
        loc.block.full <- unique(c(loc.block.full, loc.block))
      }else{
        pts.sel <- c(pts.sel);
      }
    }

    print("pts.sel:"); print(sort(pts.sel))

    phi.hat.full.new <- phi.hat.full
    for(i in 2:n.new){
      if(!(i %in% loc.block.full)){
        phi.hat.full.new[, ((i-1)*p.x.temp+1):(i*p.x.temp)] <- matrix(0, p.y, p.x.temp)
      }
    }

    phi.full.all.new <- vector("list", n.new);
    phi.full.all.new.temp <- vector("list", n.new);
    forecast.all.new <- matrix(0, p.y, TT);

    # update: 09/09/2020: use the estimation in middle index (same as second step)
    phi.hat.new <- phi.hat.full.new;
    phi.full.all.new[[1]] <- matrix(phi.hat.new[, (1):(p.x.temp)], ncol = p.x.temp);
    phi.full.all.new.temp[[1]] <- matrix(phi.hat.new[, (1):(p.x.temp)], ncol = p.x.temp);
    for(i.1 in 2:n.new){
      #phi.full.all.new.temp keeps adding
      phi.full.all.new.temp[[i.1]] <- matrix(phi.full.all.new.temp[[i.1-1]] + phi.hat.full[, ((i.1-1)*p.x.temp+1):(i.1*p.x.temp)], ncol = p.x.temp);
    }

    beta.hat.list.new <- vector("list", length(loc.block.full) + 1) ;
    # includes the boundary block, and sorted
    loc.block.new <-  c(1, sort(loc.block.full), n.new+1)
    idx.1 <- floor((loc.block.new[1] + loc.block.new[2])/2) ;
    beta.hat.list.new[[1]] <- phi.full.all.new.temp[[idx.1]]
    lb.1 <- blocks[loc.block.new[1]];
    ub.1 <- blocks[loc.block.new[2]] - 1;
    # print("lb:");print(lb.1)
    # print("ub:");print(ub.1)
    for(j.1 in 1:p.y){
      data.x.temp <- data_x[, -j.1]
      forecast.all.new[j.1, lb.1:ub.1] <- pred.block(t(data.x.temp), matrix(beta.hat.list.new[[1]][j.1, ], ncol = p.x.temp),
                                                     lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1);
    }

    print(sort(loc.block.full))
    for(i.1 in 1:length(loc.block.full)){
      idx.1 <- floor((loc.block.new[i.1 + 1] + loc.block.new[i.1 + 2])/2);
      # print(idx.1)
      beta.hat.list.new[[i.1 + 1]] <- phi.full.all.new.temp[[idx.1]]
      lb.1 <- blocks[loc.block.new[i.1 + 1]];
      ub.1 <- blocks[loc.block.new[i.1 + 2]] - 1;
      # print("lb:");print(lb.1)
      # print("ub:");print(ub.1)
      for(j.1 in 1:p.y){
        data.x.temp <- data_x[, -j.1]
        forecast.all.new[j.1, lb.1:ub.1] <- pred.block(t(data.x.temp), matrix(beta.hat.list.new[[i.1 + 1]][j.1, ], ncol = p.x.temp),
                                                       lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1);
      }
    }
    # print(forecast.all.new[1, ])

    residual <- t(data_y[(1:TT), ]) - forecast.all.new;

    if(HBIC == TRUE){
      print("Use HBIC!")
      print('gamma value:')
      print(gamma.val)

      if(is.null(gamma.val)){
        BIC.new <- 0
        for(jjj in 1:p.y){
          BIC.new <- BIC.new +
            BIC(matrix(residual[jjj, ], nrow = 1), phi = matrix(phi.hat.full.new[jjj, ], nrow = 1))$HBIC
        }
        # BIC.new <- BIC(matrix(residual[1, ], nrow = 1), phi = matrix(phi.hat.full.new[1, ], nrow = 1))$HBIC
      }else{
        BIC.new <- 0
        for(jjj in 1:p.y){
          BIC.new <- BIC.new +
            BIC(matrix(residual[jjj, ], nrow = 1), phi = matrix(phi.hat.full.new[jjj, ], nrow = 1), gamma.val = gamma.val)$HBIC
        }
        # BIC.new <- BIC(matrix(residual[1, ], nrow = 1), phi = matrix(phi.hat.full.new[1, ], nrow = 1), gamma.val = gamma.val)$HBIC
      }

    }else{
      print("Use BIC!")
      BIC.new <- 0
      for(jjj in 1:p.y){
        BIC.new <- BIC.new +
          BIC(matrix(residual[jjj, ], nrow = 1), phi = matrix(phi.hat.full.new[jjj, ], nrow = 1) )$BIC
      }
      # BIC.new <- BIC(matrix(residual[1, ], nrow = 1), phi = matrix(phi.hat.full.new[1, ], nrow = 1) )$BIC
    }
    print("BIC.new:"); print(BIC.new)
    BIC.diff <- BIC.old - BIC.new
    print("BIC.diff:"); print(BIC.diff)
    BIC.old <- BIC.new
    if(BIC.diff <= 0){
      pts.sel <- sort(pts.sel.old)
      break
    }
    # jumps[loc.block] <- max(jumps[-loc.block])
    jumps[loc.block] <- 0
    plot(jumps, type = 'o', main= 'squared L2 norm',lwd = 2,cex=2 )
  }

  # print(pts.sel)

  if ( length(pts.sel) == 0){cp.final <- c(); pts.list <-  vector("list", 0);}
  if ( length(pts.sel) > 0){
    cp.final <- pts.sel;
    # REMOVE BOUNDARY POINTS and other cleaning
    cp.final <- cp.final[which(cp.final > sum(blocks.size[1:3]))];
    cp.final <- cp.final[which(cp.final < (TT-sum(blocks.size[(length(blocks.size)-2):length(blocks.size)])))];
    cp.final <- sort(cp.final);
    # print("cp.final:"); print(cp.final)

    if(length(cp.final) >=2){
      gap.temp <- sapply(2:length(cp.final), function(jjj) cp.final[jjj]-cp.final[jjj-1])
      # print("gap.temp:");print(gap.temp)
    }


    # if there are multipler change points
    # use the k-means clustering
    if(length(cp.final) > 5  ){
      if(length(unique(gap.temp)) > 1  ){
        print(fviz_nbclust(matrix(cp.final, length(cp.final), 1), kmeans, nstart = 25,  method = "gap_stat",
                           k.max = min(50, length(cp.final)-1), nboot = 100)+
                labs(subtitle = "Gap statistic method"))
        cl <- fviz_nbclust(matrix(cp.final, length(cp.final), 1), kmeans, nstart = 25,  method = "gap_stat",
                           k.max = min(50, length(cp.final)-1), nboot = 100)+
          labs(subtitle = "Gap statistic method")
        cl.data <- cl$data;
        gap <- cl.data$gap;
        # se <- cl.data$SE.sim;
        # i.cl <- 0;
        print("choose the maximum gap stat")
        cl.number <- which.max(gap)
        gap.order <- order(gap, decreasing = TRUE)
        print("check if the diameter of clusters is small enough")

        if(median(blocks.size) <= sqrt(TT)/4){
          cnst <- 2*4+1
        }else if(median(blocks.size) <= 2*sqrt(TT)/4){
          cnst <- 2*3+1
        }else if(median(blocks.size) <= 3*sqrt(TT)/4){
          cnst <- 2*2+1
        }else{
          cnst <- 2*2+1
        }

        flag = TRUE
        idx <- 1
        while(flag & idx <= length(gap.order)){
          cl.number <- gap.order[idx]
          print("cl.number:")
          print(cl.number)
          if(cl.number > 1){
            cluster.pos <- sort(order(gap.temp, decreasing = TRUE)[1:(cl.number-1)])
            cluster.pos <- c(0, cluster.pos, length(cp.final))
          }else{
            cluster.pos <- c(0, length(cp.final))
          }

          wide <- 0
          for (i in c(1:cl.number)) {
            pts.i <- cp.final[(cluster.pos[i]+1): cluster.pos[i+1]]
            print(paste0("Cluster: ", i))
            print(pts.i)
            block_avg <- mean(blocks.size[match(pts.i, blocks)])
            # if ( max(pts.i) - min(pts.i) > 5*mean(blocks.size) ){
            if ( max(pts.i) - min(pts.i) > cnst*block_avg ){
              print("Too wide, need seperate!")
              wide <- 1
              break
            }
          }
          # print(idx)
          if (wide){
            idx <- idx + 1
          }else{
            idx <- idx
            flag <- FALSE
          }
        }
        if(flag){
          print("they are seperate change points! not use gap stat!")
          cl.number <- length(gap.temp) + 1
        }

        # print(cl.number)
        print("check if distance between two clusters are large enough")
        if(median(blocks.size) <= sqrt(TT)/4){
          cnst2 <- 8
        }else if(median(blocks.size) <= sqrt(TT)/2){
          cnst2 <- 5
        }else{
          cnst2 <- 3
        }

        if(cl.number > 1){
          cl.number <- sum(sort(gap.temp, decreasing = TRUE)[1:(cl.number - 1)] > cnst2*median(blocks.size)) + 1
        }

      }else if(unique(gap.temp) == median(blocks.size) ){
          print("one single cluster")
          cl.number <- 1
      }else{
        print("they are seperate change points! not use gap stat!")
        cl.number <- length(gap.temp) + 1
      }
      print('number of clusters:')
      print(cl.number)






      # if(length(unique(gap.temp)) > 1 ){
      #   print(fviz_nbclust(matrix(cp.final,length(cp.final),1), kmeans, nstart = 25,  method = "gap_stat", k.max = min(12, length(cp.final)-1), nboot = 100)+
      #           labs(subtitle = "Gap statistic method"))
      #   cl <- fviz_nbclust(matrix(cp.final,length(cp.final),1), kmeans, nstart = 25,  method = "gap_stat", k.max = min(12, length(cp.final)-1), nboot = 100)+
      #     labs(subtitle = "Gap statistic method")
      #
      #   cl.data <- cl$data;
      #   gap <- cl.data$gap;
      #   se <- cl.data$SE.sim;
      #   i.cl <- 0;
      #   print("choose the maximum gap stat")
      #   cl.number <- which.max(gap)
      #
      # }else{
      #   print("they are seperate change points! not use gap stat!")
      #   cl.number <- length(gap.temp) + 1
      # }

      #cl.number

      print("use gap instead of kmeans!")
      if(cl.number > 1){
        cluster.pos <- sort(order(gap.temp, decreasing = TRUE)[1:(cl.number-1)])
        cluster.pos <- c(0, cluster.pos, length(cp.final))

      }else{
        cluster.pos <- c(0, length(cp.final))
      }

      pts.list <-  vector("list", cl.number);
      for (i in c(1:cl.number)) {
        pts.i <- cp.final[(cluster.pos[i]+1): cluster.pos[i+1]]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }

      # cl.final <- kmeans(cp.final, centers = cl.number);
      # pts.list <-  vector("list",cl.number);
      # loc.new <- cl.final$cluster;
      # cl.reorder = c(1:cl.number)[order(cl.final$centers)]
      # for (i in c(1:cl.number)) {
      #   pts.i <- cp.final[which(loc.new==cl.reorder[i])]
      #   print(paste0("Cluster: ", i))
      #   print(pts.i)
      #   pts.list[[i]] <- pts.i
      # }
    }

    #!!!!!!!!!!!!!!!!need  ajustment!!!!!!
    if(length(cp.final) <= 5 & length(cp.final) > 1 ){
      print("small number of cp !!!!")
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))

      for (i in 2:length(cp.final)){
        if (cp.final[i]-cp.final[i-1]<= max(3*mean(blocks.size))){
          cl.number <-cl.number-1
          loc.new[i] <-loc.new[i-1]
        }else{
          loc.new[i] <- i
        }
      }

      pts.list <-  vector("list",cl.number);
      #for (i in unique(loc.new)) {
      loc.new.unique <- unique(loc.new)
      for (i in 1:length(loc.new.unique)) {
        pts.i <- cp.final[which(loc.new==loc.new.unique[i])]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }


    if(length(cp.final) == 1 ){
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))
      pts.list <-  vector("list",cl.number);
      for (i in unique(loc.new)) {
        pts.i <- cp.final[which(loc.new==i)]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }

    if(length(cp.final) == 0 ){
      pts.list <-  vector("list", 0);
    }

  }

  #compute the estimated beta
  phi.par.sum <- vector("list",n.new);
  phi.par.sum[[1]] <- phi.hat.full[, 1:(p.x.temp)];
  for(i in 2:n.new){
    phi.par.sum[[i]] <- phi.par.sum[[i-1]] + phi.hat.full[,((i-1)*p.x.temp+1):(i*p.x.temp)];
  }


  print("First step DONE!!!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  return(list(jumps.l2 = jumps.l2, jumps.l1 = jumps.l1,
              pts.list = pts.list, beta.full = phi.par.sum ))

}


#' Exhaustive search step for gaussian graphical model.
#'
#' @description Perform the exhaustive search to "thin out" redundant break points.
#'
#' @param data_y input data matrix, with each column representing the time series component
#' @param data_x input data matrix, with each column representing the time series component
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cp.first the selected break points after the first step
#' @param beta.est the estiamted parameters by block fused lasso
#' @param blocks the blocks
#' @return A list oject, which contains the followings
#' \describe{
#'   \item{cp.final}{a set of selected break point after the exhaustive search step}
#'   \item{beta.hat.list}{the estimated coefficient matrix for each segmentation}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
ggm.second.step.search <- function(data_y, data_x, max.iteration = max.iteration, tol = tol,  cp.first, beta.est, blocks){

  TT <- length(data_y[,1]); p.y <- length(data_y[1,]); p.x <- length(data_x[1, ]);
  p.x.temp <- p.x - 1;
  p.y.temp <- 1;


  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );

  cp.search <- cp.first;
  cl.number <- length(cp.first)

  cp.list <- vector("list", cl.number+2);
  cp.list[[1]] <- c(1);
  cp.list[[cl.number+2]] <- c(TT+1);

  cp.index.list <- vector("list", cl.number+2);
  cp.index.list[[1]] <- c(1);
  cp.index.list[[cl.number+2]] <- c(n.new+1);

  for (i in 1:cl.number) {
    cp.list[[i+1]] <- cp.first[[i]]
    cp.index.list[[i+1]] <- match(cp.first[[i]], blocks)
  }


  cp.search <- rep(0, cl.number);
  cp.list.full <- cp.list

  beta.hat.list <- vector("list", cl.number+1)


  for(i in 1:(cl.number)){
    idx <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2);
    beta.hat.list[[i]] <- beta.est[[idx]]

    if(length(cp.list[[i+1]]) > 1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] + 1)  :  (cp.list[[i+1]][length(cp.list[[i+1]])]-1  ) )
    }
    if(length(cp.list[[i+1]]) == 1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] -  (blocks.size[cp.index.list[[i+1]][1] ]) + 1) :  (cp.list[[i+1]][length(cp.list[[i+1]])] +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1   ) )
    }


    #compare the SSE of first num and last num
    # Start from the left
    num  <- cp.list.full[[i+1]][1]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    beta.hat <- beta.est[[idx.1]]

    forecast <- matrix(0, ncol = p.y, nrow = len.1)
    for(j.1 in 1:p.y){
      data.x.temp <- data_x[, -j.1]
      forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1)
    }
    temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    # forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x), jjj, p.x, p.y)  )
    # if(len.1 == 1){
    #   temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    # }else{
    #   temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
    # }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    beta.hat <- beta.est[[idx.2]]
    forecast <- matrix(0, ncol = p.y, nrow = len.2)
    for(j.1 in 1:p.y){
      data.x.temp <- data_x[, -j.1]
      forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.2, p.x.temp, p.y.temp, ub.2 - lb.2 + 1)
    }
    temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    # forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    # if(len.2 == 1){
    #   temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    # }else{
    #   temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
    # }
    sse1 <- temp.1 + temp.2;

    # Start from the right
    num  <- cp.list.full[[i+1]][length(cp.list.full[[i+1]])]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    beta.hat <- beta.est[[idx.1]]
    forecast <- matrix(0, ncol = p.y, nrow = len.1)
    for(j.1 in 1:p.y){
      data.x.temp <- data_x[, -j.1]
      forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1)
    }
    temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    # forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    # if(len.1 == 1){
    #   temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    # }else{
    #   temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
    # }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    beta.hat <- beta.est[[idx.2]]
    forecast <- matrix(0, ncol = p.y, nrow = len.2)
    for(j.1 in 1:p.y){
      data.x.temp <- data_x[, -j.1]
      forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.2, p.x.temp, p.y.temp, ub.2 - lb.2 + 1)
    }
    temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    # forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
    # if(len.2 == 1){
    #   temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    # }else{
    #   temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
    # }
    sse2 <- temp.1 + temp.2;



    if(sse1 <= sse2){
      sse.full <- 0;
      ii <- 0
      # print(cp.list.full[[i+1]] )
      for(num in cp.list.full[[i+1]]  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        beta.hat <- beta.est[[idx.1]]
        forecast <- matrix(0, ncol = p.y, nrow = len.1)
        for(j.1 in 1:p.y){
          data.x.temp <- data_x[, -j.1]
          forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1)
        }
        temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
        # forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
        # if(len.1 == 1){
        #   temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
        #
        # }else{
        #   temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
        # }


        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        beta.hat <- beta.est[[idx.2]]
        forecast <- matrix(0, ncol = p.y, nrow = len.2)
        for(j.1 in 1:p.y){
          data.x.temp <- data_x[, -j.1]
          forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.2, p.x.temp, p.y.temp, ub.2 - lb.2 + 1)
        }
        temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );

        # forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x),matrix(beta.hat,ncol = p.x),jjj, p.x, p.y)  )
        # if(len.2 == 1){
        #   temp.2 <- sum( (data_y[lb.2:ub.2,]-forecast)^2 );
        # }else{
        #   temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
        # }
        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full,0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][min(which(sse.full == min(sse.full)))];

    }
    if(sse1 > sse2){
      sse.full <- 0;
      ii <- 0
      # print(rev(cp.list.full[[i+1]]) )
      for(num in rev(cp.list.full[[i+1]])  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        beta.hat <- beta.est[[idx.1]]
        forecast <- matrix(0, ncol = p.y, nrow = len.1)
        for(j.1 in 1:p.y){
          data.x.temp <- data_x[, -j.1]
          forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.1, p.x.temp, p.y.temp, ub.1 - lb.1 + 1)
        }
        temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
        # forecast <- sapply(c(lb.1:ub.1), function(jjj) pred(t(data_x), matrix(beta.hat, ncol = p.x),jjj, p.x, p.y)  )
        # if(len.1 == 1){
        #   temp.1 <- sum( (data_y[lb.1:ub.1, ]-forecast)^2 );
        #
        # }else{
        #   temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
        # }


        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        beta.hat <- beta.est[[idx.2]]
        forecast <- matrix(0, ncol = p.y, nrow = len.2)
        for(j.1 in 1:p.y){
          data.x.temp <- data_x[, -j.1]
          forecast[, j.1] <- pred.block(t(data.x.temp), matrix(beta.hat[j.1, ], ncol = p.x.temp), lb.2, p.x.temp, p.y.temp, ub.2 - lb.2 + 1)
        }
        temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
        # forecast <- sapply(c(lb.2:ub.2), function(jjj) pred(t(data_x),matrix(beta.hat,ncol = p.x),jjj, p.x, p.y)  )
        # if(len.2 == 1){
        #   temp.2 <- sum( (data_y[lb.2:ub.2,]-forecast)^2 );
        # }else{
        #   temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
        # }
        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full, 0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][length(cp.list.full[[i+1]]) + 1 - min(which(sse.full == min(sse.full)))];

    }

  }

  print("cp.final:")
  print(cp.search)
  print("Second step DONE!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  idx <- floor((min(cp.index.list[[cl.number+1+1]]) + max(cp.index.list[[cl.number+1]]))/2);
  beta.hat.list[[cl.number+1]] <- beta.est[[idx]]
  return(list(cp.final = cp.search, beta.hat.list = beta.hat.list))
}


################################################
#' Generating non-stationary ARMA data.
#'
#' @param nobs number of time points
#' @param arlags the true AR order
#' @param malags the true MA order
#' @param cnst the constant
#' @param phi parameter matrix of the AR model
#' @param theta parameter matrix of the MA model
#' @param skip the number of time points to skip at the begining (for stable data)
#' @param sigma covariance matrix of the white noise
#' @param brk vector of break points
#' @return Matrice of time series data and white noise data
#' @import mvtnorm
#' @export
var.sim.break <- function (nobs, arlags = NULL, malags = NULL, cnst = NULL, phi = NULL, theta = NULL,
                           skip = 200, sigma, brk = nobs+1) {
  if (!is.matrix(sigma))
    sigma = as.matrix(sigma)
  k <- nrow(sigma); m <- length(brk); nT <- nobs + skip

  #generate multivariate normal distributed data as the white noise data
  at <- rmvnorm(nT, rep(0, k), sigma)

  #generate the ARMA time series data
  nar <- length(arlags); p <- 0
  if (nar > 0) {
    arlags <- sort(arlags)
    p <- arlags[nar]
  }

  nma <- length(malags); q <- 0
  if (nma > 0) {
    malags <- sort(malags)
    q <- malags[nma]
  }

  ist = max(p, q) + 1
  zt = matrix(0, nT, k)
  if (length(cnst) == 0)
    cnst = rep(0, k)
  if (m == 1){
    for (it in ist:nT) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
  }

  #if there are some break points
  if (m > 1){
    for (it in ist:(skip+brk[1]-1)) {
      tmp = matrix(at[it, ], 1, k)
      if (nma > 0) {
        for (j in 1:nma) {
          jdx = (j - 1) * k
          thej = theta[, (jdx + 1):(jdx + k)]
          atm = matrix(at[it - malags[j], ], 1, k)
          tmp = tmp - atm %*% t(thej)
        }
      }
      if (nar > 0) {
        for (i in 1:nar) {
          idx = (i - 1) * k
          phj = phi[, (idx + 1):(idx + k)]
          ztm = matrix(zt[it - arlags[i], ], 1, k)
          tmp = tmp + ztm %*% t(phj)
        }
      }
      zt[it, ] = cnst + tmp
    }
    for ( mm in 1:(m-1)){
      for (it in (skip+brk[mm]):(skip+brk[mm+1]-1) ) {
        tmp = matrix(at[it, ], 1, k)
        if (nma > 0) {
          for (j in 1:nma) {
            jdx = (j - 1) * k
            thej = theta[, (jdx + 1):(jdx + k)]
            atm = matrix(at[it - malags[j], ], 1, k)
            tmp = tmp - atm %*% t(thej)
          }
        }
        if (nar > 0) {
          for (i in 1:nar) {
            idx = (i - 1) * k
            phj = phi[, ((mm)*p*k+idx + 1):((mm)*p*k+idx + k)]
            ztm = matrix(zt[it - arlags[i], ], 1, k)
            tmp = tmp + ztm %*% t(phj)
          }
        }
        zt[it, ] = cnst + tmp
      }
    }
  }

  zt = zt[(1 + skip):nT, ]
  at = at[(1 + skip):nT, ]
  VARMAsim <- list(series = zt, noises = at)
}



################################################
#' Threshold block fused lasso step for linear regression model.
#'
#' @description Perform the block fused lasso with thresholding to detect candidate break points.
#'
#' @param data_y input data matrix Y, with each column representing the time series component
#' @param lambda1 tuning parmaeter lambda_1 for fused lasso
#' @param lambda2 tuning parmaeter lambda_2 for fused lasso
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cv.index the index of time points for cross-validation
#' @param blocks the blocks
#' @param HBIC logical; if TRUE, use high-dimensional BIC, if FALSE, use orginal BIC. Default is FALSE.
#' @param gamma.val hyperparameter for HBIC, if HBIC == TRUE.
#' @return A list object, which contains the followings
#' \describe{
#'   \item{jump.l2}{estimated jump size in L2 norm}
#'   \item{jump.l1}{estimated jump size in L1 norm}
#'   \item{pts.list}{estimated change points in the first step}
#'   \item{phi.full}{estimated parameters in the first step}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
#' @import factoextra
var.first.step.blocks <- function(data_y, lambda1, lambda2, q, max.iteration, tol,
                                  blocks, cv.index, HBIC = FALSE, gamma.val = NULL){

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );

  #create the tuning parmaeter combination of lambda1 and lambda2
  lambda.full <- expand.grid(lambda1,lambda2)
  kk <- length(lambda.full[,1]);

  cv.l <- length(cv.index);
  cv <- rep(0,kk);
  phi.final <- vector("list",kk);

  data.y.temp <- data_y;
  TT <- length(data.y.temp[,1]);
  p <- length(data.y.temp[1,]);

  flag.full <- rep(0,kk);

  for (i in 1:kk) {
    print(i)
    print("##################lambda1:")
    print(lambda.full[i,1])
    print("##################lambda2:")
    print(lambda.full[i,2])
    if ( i == 1){
      test <- var_break_fit_block(data.y.temp, lambda.full[i,1],lambda.full[i,2], q, max_iteration = max.iteration, tol = tol, initial_phi =  0.0+matrix(0.0,p,p*q*n.new), blocks = blocks, cv.index)
      flag.full[i] <- test$flag;
    }
    if ( i > 1 ){
      initial.phi <- phi.final[[i-1]]
      if(max(abs(phi.final[[(i-1)]])) > 10^3  ){initial.phi <- 0*phi.final[[(i-1)]];}
      test <- var_break_fit_block(data.y.temp, lambda.full[i,1], lambda.full[i,2], q, max_iteration = max.iteration, tol = tol, initial_phi =  initial.phi, blocks = blocks, cv.index)
      flag.full[i] <- test$flag;
    }

    phi.hat.full <- test$phi.hat;
    phi.final[[i]] <- phi.hat.full;


    #forecast the time series based on the estimated matrix Phi (beta for the linear regression model)
    #and compute the forecast error
    phi.full.all <- vector("list",n.new);
    forecast <- matrix(0,p,TT);
    phi.hat <- phi.hat.full;
    phi.full.all[[1]] <- phi.hat[,(1):(p*q)];
    for(i.1 in 2:n.new){
      phi.full.all[[i.1]] <- phi.full.all[[i.1-1]] + phi.hat[,((i.1-1)*p*q+1):(i.1*p*q)];
      #forecast[,(blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block.var(t(data_y), phi.full.all[[i.1-1]], q, blocks[i.1], p, blocks[i.1+1]-blocks[i.1]);
    }
    forecast.new <- matrix(0, p, cv.l);
    for(j in (1):cv.l){
      forecast.new[,j] <- pred.var(t(data_y), phi.full.all[[(cv.index[j])]], q, blocks[cv.index[j]+1]-1, p, 1)
    }
    temp.index <- rep(0,cv.l);
    for(ff in 1:cv.l){temp.index[ff] <- blocks[cv.index[ff]+1]-1;}
    cv[i] <- (1/(p*cv.l))*sum( (forecast.new - t(data_y[temp.index,])  )^2 );


    print("============cv-result=======================")
    print(cv[i])
    print("====================================")
  }



  lll <- min(which(cv==min(cv)));
  ind.new <- 0;
  lll <- lll - ind.new;

  mspe.plot(cv, c(1:kk))
  abline(v = seq(length(lambda1), length(lambda1)*(length(lambda2)-1), length.out =length(lambda2)-1)+0.5)
  abline(v= lll, col="red")


  phi.hat.full <- phi.final[[lll]];

  #jumps.sq is the L2 norm square
  #jumps.l1 is the L1 norm
  # again, note that here the phi.hat.full is the estimated theta in the paper
  jumps.l2 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l2[i] <- (sum((phi.hat.full[,((i-1)*p*q+1):(i*p*q)] )^2 ));
  }
  jumps.l1 <- rep(0,n.new);
  for(i in c(2:n.new)){
    jumps.l1[i] <- (sum(abs(phi.hat.full[,((i-1)*p*q+1):(i*p*q)] ) ));
  }

  # print(jumps.l2)
  #print(jumps.l1)
  # plot(jumps.l2,type = 'o', main= 'l2 norm')
  #plot(jumps.l1,type = 'o', main= 'l1 norm')

  #ignore the large jump at the boundary!!!!!!!!!!!!!!!!!!!!
  # print('use l2 norm!')
  jumps <- jumps.l2
  #print('use l1 norm!')
  #jumps <- jumps.l1
  ignore_num <- max(6, round(250/min(blocks.size)))
  jumps[1:ignore_num]  <- 0
  jumps[(length(jumps)-(ignore_num-1)):length(jumps)]  <- 0
  plot(jumps,type = 'o', main= 'l2 norm')
  ##################################################################
  #use kmeans to hard threshold the jumps
  ###### use BIC to determine the k-means
  BIC.diff <- 1
  BIC.old <- 10^8
  pts.sel <- c()
  loc.block.full <- c()
  while(BIC.diff > 0 & length(unique(jumps)) > 1 ){
    pts.sel.old <- pts.sel
    #use kmeans to hard threshold the jumps
    if( length(unique(jumps)) > 2 ){
      # print("consider 2 clusters for fit.2")
      clus.2 <- kmeans(jumps, centers = 2); fit.2 <- clus.2$betweenss/clus.2$totss;
      # print(fit.2);
      if(fit.2 < 0.20){
        print("no significant jumps!!")
        pts.sel <- c(pts.sel);
      }
      if( fit.2 >= 0.20 ){
        loc <- clus.2$cluster;
        if( clus.2$centers[1] > clus.2$centers[2]  ){
          loc.block <- which(loc==1);
        }
        if( clus.2$centers[1] < clus.2$centers[2]  ){
          loc.block <- which(loc==2);
        }
        pts.sel <- c(pts.sel, blocks[loc.block]);
        loc.block.full <- c(loc.block.full, loc.block)
      }
    }
    if( length(unique(jumps)) <= 2 ){
      if(length(unique(jumps)) == 2){
        loc.block <- which.max(jumps)
        pts.sel <- c(pts.sel, blocks[loc.block])
        loc.block.full <- c(loc.block.full, loc.block)
      }else{
        pts.sel <- c(pts.sel);
      }
    }

    # print("pts.sel:"); print(sort(pts.sel))

    phi.hat.full.new <- phi.hat.full
    for(i in 2:n.new){
      if(!(i %in% loc.block.full)){
        phi.hat.full.new[, ((i-1)*p*q+1):(i*p*q)] <- matrix(0, p, p*q)
      }
    }

    phi.full.all.new <- vector("list", n.new);
    phi.full.all.new.temp <- vector("list", n.new);
    forecast.all.new <- matrix(0, p, TT);

    phi.hat.new <- phi.hat.full.new;
    phi.full.all.new[[1]] <- matrix(phi.hat.new[, (1):(p*q)], ncol = p*q);
    phi.full.all.new.temp[[1]] <- matrix(phi.hat.new[, (1):(p*q)], ncol = p*q);



    forecast.all.new[, (blocks[1]+q):(blocks[2]-1)] <- pred.block.var(t(data_y), phi.full.all.new[[1]], q,
                                                                blocks[1]+q, p, blocks[2] - (blocks[1]+q) );
    for(i.1 in 2:n.new){
      #phi.full.all.new.temp keeps adding
      phi.full.all.new.temp[[i.1]] <- matrix(phi.full.all.new.temp[[i.1-1]] + phi.hat.full[, ((i.1-1)*p*q+1):(i.1*p*q)], ncol = p*q);
      if((i.1 %in% loc.block.full)){
        phi.full.all.new[[i.1]] <- phi.full.all.new.temp[[i.1]]
      }
      if(!(i.1 %in% loc.block.full)){
        phi.full.all.new[[i.1]] <- phi.full.all.new[[i.1-1]]
      }
      forecast.all.new[, (blocks[i.1]):(blocks[i.1+1]-1)] <- pred.block.var(t(data_y), phi.full.all.new[[i.1]], q,
                                                                        blocks[i.1], p, blocks[i.1+1] - blocks[i.1]);
    }

    residual <- t(data.y.temp[( (1+q) :TT), ]) - forecast.all.new[, (1+q) :TT];
    # print(phi.full.all.new)
    # print(phi.full.all.new.temp)
    # print(residual)

    if(HBIC == TRUE){
      print("Use HBIC!")
      # print('gamma value:')
      # print(gamma.val)
      if(is.null(gamma.val)){
        BIC.new <- BIC(residual, phi = phi.hat.full.new, method = 'VAR')$HBIC
      }else{
        BIC.new <- BIC(residual, phi = phi.hat.full.new, gamma.val = gamma.val, method = 'VAR')$HBIC
      }

    }else{
      print("Use BIC!")
      BIC.new <- BIC(residual, phi = phi.hat.full.new, method = 'VAR' )$BIC
    }
    # print("BIC.new:"); print(BIC.new)
    BIC.diff <- BIC.old - BIC.new
    # print("BIC.diff:");print(BIC.diff)
    BIC.old <- BIC.new
    if(BIC.diff <= 0){
      pts.sel <- sort(pts.sel.old)
      break
    }
    jumps[loc.block] <- 0
    plot(jumps, type = 'o', main= 'l2 norm')
  }

  # print(pts.sel)



  if ( length(pts.sel) == 0){cp.final <- c(); pts.list <-  vector("list",0);}
  if ( length(pts.sel) > 0){
    cp.final <- pts.sel;
    # REMOVE BOUNDARY POINTS and other cleaning
    cp.final <- cp.final[which(cp.final > sum(blocks.size[1:3]))];
    cp.final <- cp.final[which(cp.final < (TT-sum(blocks.size[(length(blocks.size)-2):length(blocks.size)])))];
    cp.final <- sort(cp.final);
    # print(cp.final)


    # if there are multipler change points
    # use the p-means clustering
    if(length(cp.final) > 4){
      print(fviz_nbclust(matrix(cp.final,length(cp.final),1), kmeans, nstart = 25,  method = "gap_stat", k.max = min(12,length(cp.final)-1), nboot = 100)+
              labs(subtitle = "Gap statistic method"))
      cl <- fviz_nbclust(matrix(cp.final,length(cp.final),1), kmeans, nstart = 25,  method = "gap_stat", k.max = min(12,length(cp.final)-1), nboot = 100)+
        labs(subtitle = "Gap statistic method")

      cl.data <- cl$data;
      gap <- cl.data$gap;
      # se <- cl.data$SE.sim;
      # i.cl <- 0;
      # while (i.cl < (length(gap)-1)) {
      #   i.cl <- i.cl + 1;
      #   if( gap[i.cl] > gap[i.cl+1] - se[i.cl+1] ){cl.number <- i.cl; break;}
      # }
      # print("choose the maximum gap stat")
      cl.number <- which.max(gap)
      # print(cl.number)

      cl.final <- kmeans(cp.final, centers = cl.number);
      pts.list <-  vector("list",cl.number);
      loc.new <- cl.final$cluster;
      cl.reorder = c(1:cl.number)[order(cl.final$centers)]
      for (i in c(1:cl.number)) {
        pts.i <- cp.final[which(loc.new==cl.reorder[i])]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }

    #!!!!!!!!!!!!!!!!need  ajustment!!!!!!
    if(length(cp.final) <= 4 & length(cp.final) >1 ){
      print("small number of cp !!!!")
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))

      for (i in 2:length(cp.final)){
        if (cp.final[i]-cp.final[i-1]<= max(3*mean(blocks.size))){
          cl.number <-cl.number-1
          loc.new[i] <-loc.new[i-1]
        }else{
          loc.new[i] <- i
        }
      }

      pts.list <-  vector("list",cl.number);
      #for (i in unique(loc.new)) {
      loc.new.unique <- unique(loc.new)
      for (i in 1:length(loc.new.unique)) {
        pts.i <- cp.final[which(loc.new==loc.new.unique[i])]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }

    if(length(cp.final) == 1 ){
      cl.number <- length(cp.final);
      loc.new <- rep(1,length(cp.final))
      pts.list <-  vector("list",cl.number);
      for (i in unique(loc.new)) {
        pts.i <- cp.final[which(loc.new==i)]
        print(paste0("Cluster: ", i))
        print(pts.i)
        pts.list[[i]] <- pts.i
      }
    }

    if(length(cp.final) == 0 ){
      pts.list <-  vector("list", 0);
    }

  }

  #compute the estimated phi
  phi.par.sum <- vector("list",n.new);
  phi.par.sum[[1]] <- phi.hat.full[,1:(p*q)];
  for(i in 2:n.new){
    phi.par.sum[[i]] <- phi.par.sum[[i-1]] + phi.hat.full[,((i-1)*p*q+1):(i*p*q)];
  }


  #plot(blocks[1:n.new],jumps,main = "JUMPS.FULL", type = "o")

  print("First step DONE!!!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
  return(list(jumps.l2 = jumps.l2, jumps.l1 = jumps.l1,
              pts.list = pts.list, phi.full = phi.par.sum))

}


#' Exhaustive search step
#'
#' @description Perform the exhaustive search to "thin out" redundant break points.
#'
#' @param data_y input data matrix, with each column representing the time series component
#' @param q the AR order
#' @param max.iteration max number of iteration for the fused lasso
#' @param tol tolerance for the fused lasso
#' @param cp.first the selected break points after the first step
#' @param beta.est the estimated parameters by block fused lasso
#' @param blocks the blocks
#' @return A list object, which contains the followings
#' \describe{
#'   \item{cp.final}{a set of selected break point after the exhaustive search step}
#'   \item{phi.hat.list}{the estimated coefficient matrix for each segmentation}
#' }
#' @import graphics
#' @import ggplot2
#' @import stats
var.second.step.search <- function(data_y, q, max.iteration = max.iteration, tol = tol,
                                   cp.first, beta.est, blocks){

  TT <- length(data_y[,1]); p.y <- length(data_y[1,]);
  p <- length(data_y[1,]);

  n.new <- length(blocks) - 1;
  blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj]  );


  cp.search <- cp.first;
  cl.number <- length(cp.first)

  cp.list <- vector("list",cl.number+2);
  cp.list[[1]] <- c(1);
  cp.list[[cl.number+2]] <- c(TT+1);

  cp.index.list <- vector("list",cl.number+2);
  cp.index.list[[1]] <- c(1);
  cp.index.list[[cl.number+2]] <- c(n.new+1);

  for (i in 1:cl.number) {
    cp.list[[i+1]] <- cp.first[[i]]
    cp.index.list[[i+1]] <- match(cp.first[[i]], blocks)
  }


  cp.search <- rep(0,cl.number);
  cp.list.full <- cp.list

  phi.hat.list <- vector("list", cl.number+1)
  for(i in 1:(cl.number)){
    idx <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2);
    phi.hat.list[[i]] <- beta.est[[idx]]



    if(length(cp.list[[i+1]]) > 1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1] +1) :  (cp.list[[i+1]][length(cp.list[[i+1]])] -1 ) )
    }
    if(length(cp.list[[i+1]]) == 1){
      cp.list.full[[i+1]] <- c((cp.list[[i+1]][1]-  (blocks.size[cp.index.list[[i+1]][1] ])+1 ) :  (cp.list[[i+1]][length(cp.list[[i+1]])] +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1  ) )
    }
    #print(cp.list.full[[i+1]] )

    #compare the SSE of first num and last num
    num  = cp.list.full[[i+1]][1]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    phi.hat <- beta.est[[idx.1]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
    if(len.1 == 1){
      temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
    }


    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    phi.hat <- beta.est[[idx.2]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
    if(len.2 == 1){
      temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
    }

    sse1 <- temp.1 + temp.2;


    num  <- cp.list.full[[i+1]][length(cp.list.full[[i+1]])]
    lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
    ub.1 <- num - 1;
    len.1 <- ub.1 - lb.1 + 1;
    idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
    phi.hat <- beta.est[[idx.1]]
    forecast <- sapply(c(lb.1:ub.1), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
    if(len.1 == 1){
      temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
    }else{
      temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
    }

    lb.2 <- num ;
    ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
    len.2 <- ub.2 - lb.2 + 1;
    idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
    phi.hat <- beta.est[[idx.2]]
    forecast <- sapply(c(lb.2:ub.2), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
    if(len.2 == 1){
      temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
    }else{
      temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
    }
    sse2 <- temp.1 + temp.2;


    if(sse1 <= sse2){
      sse.full <- 0;
      ii <- 0
      # print(cp.list.full[[i+1]])
      for(num in cp.list.full[[i+1]]  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        phi.hat <- beta.est[[idx.1]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
        if(len.1 == 1){
          temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
        }else{
          temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
        }

        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        phi.hat <- beta.est[[idx.2]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
        if(len.2 == 1){
          temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
        }
        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full, 0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][min(which(sse.full == min(sse.full)))];
    }
    if(sse1 > sse2){
      sse.full <- 0;
      ii <- 0
      # print(rev(cp.list.full[[i+1]]) )
      for(num in rev(cp.list.full[[i+1]])  ){
        ii <- ii + 1
        lb.1 <- min(cp.list[[i+1]]) - (blocks.size[cp.index.list[[i+1]][1] ]);
        ub.1 <- num - 1;
        len.1 <- ub.1 - lb.1 + 1;
        idx.1 <- floor((min(cp.index.list[[i+1]]) + max(cp.index.list[[i]]))/2) ;
        phi.hat <- beta.est[[idx.1]]
        forecast <- sapply(c(lb.1:ub.1), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
        if(len.1 == 1){
          temp.1 <- sum( (data_y[lb.1:ub.1,]-forecast)^2 );
        }else{
          temp.1 <- sum( (t(data_y[lb.1:ub.1,])-forecast)^2 );
        }

        lb.2 <- num ;
        ub.2 <- max(cp.list[[i+1]]) +   (blocks.size[cp.index.list[[i+1]][length(cp.list[[i+1]])] ]) -1;
        len.2 <- ub.2 - lb.2 + 1;
        idx.2 <- floor((min(cp.index.list[[i+2]]) + max(cp.index.list[[i+1]]))/2) ;
        phi.hat <- beta.est[[idx.2]]
        forecast <- sapply(c(lb.2:ub.2), function(jjj) pred.var(t(data_y), matrix(phi.hat, ncol = p*q), q, jjj, p, 1)  )
        if(len.2 == 1){
          temp.2 <- sum( ( data_y[lb.2:ub.2,]-forecast)^2 );
        }else{
          temp.2 <- sum( (t(data_y[lb.2:ub.2,])-forecast)^2 );
        }

        sse.full[ii] <- temp.1 + temp.2;
        # print(ii)
        # print(sse.full[ii])
        if(ii >= min(20, length(cp.list.full[[i+1]])) && sse.full[ii] >=  quantile(sse.full, 0.20) ){
          break
        }
      }
      cp.search[i] <- cp.list.full[[i+1]][length(cp.list.full[[i+1]]) + 1 - min(which(sse.full == min(sse.full)))];
    }

  }
  idx <- floor((min(cp.index.list[[cl.number+2]]) + max(cp.index.list[[cl.number+1]]))/2);
  phi.hat.list[[cl.number+1]] <- beta.est[[idx]]

  print("cp.final:")
  print(cp.search)
  print("Second step DONE!!!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

  return(list(cp.final = cp.search, phi.hat.list = phi.hat.list))
}




#' Prediction function for VAR (block)
#' @param Y data for prediction
#' @param phi parameter matrix
#' @param q the AR order
#' @param TT the start time point for prediction
#' @param p the number of time series components
#' @param h the length of observation to predict
#' @return prediction matrix
#'
pred.block.var <- function(Y, phi, q, TT, p, h){
  concat.Y <- matrix(0, p, q+h); concat.Y[, 1:q] <- Y[, (TT-q):(TT-1)];
  for ( j in 1:h){
    temp <- matrix(0,p,1);
    for (i in 1:q){temp <- temp +  phi[,((i-1)*p+1):(i*p)]%*%concat.Y[,q+j-i];}
    concat.Y[,q+j] <- temp;
  }
  return(as.matrix(concat.Y[,(q+1):(q+h)]))
}

#' Prediction function for VAR 2
#' @param Y data for prediction
#' @param phi parameter matrix
#' @param q the AR order
#' @param TT the start time point for prediction
#' @param p the number of time series components
#' @param h the length of observation to predict
#' @return prediction matrix
#'
pred.var <- function(Y, phi, q, TT, p, h = 1){
  concat.Y <- matrix(0,p,q+h); concat.Y[,1:q] <- Y[, (TT-q):(TT-1)];
  for ( j in 1:h){
    temp <- matrix(0,p,1);
    for (i in 1:q){temp <- temp +  phi[,((i-1)*p+1):(i*p)]%*%concat.Y[,q+j-i];}
    concat.Y[,q+j] <- temp;
  }
  return(as.matrix(concat.Y[,q+h]))
}


#' prediction function
#' @param X data for prediction
#' @param phi parameter matrix
#' @param j the start time point for prediction
#' @param p.x the dimension of data X
#' @param p.y the dimension of data Y
#' @param h the length of observation to predict
#' @return prediction matrix
#'
pred <- function(X, phi, j, p.x, p.y, h = 1){
  concat.X <- matrix(0, p.x, 1);
  concat.Y <- matrix(0, p.y, 1);
  concat.X[,1] <- as.matrix(X[, j]);
  temp <- matrix(0, p.y, 1);
  if(p.y == 1){
    temp <- temp +  t(as.matrix(phi[, 1:p.x]))%*%concat.X[, 1];
  }else{
    temp <- temp +  as.matrix(phi[, 1:p.x])%*%concat.X[, 1];

  }
  concat.Y[, 1] <- temp;
  return(as.matrix(concat.Y[, 1]))
}


#' Prediction function (block)
#' @param X data for prediction
#' @param phi parameter matrix
#' @param j the start time point for prediction
#' @param p.x the dimension of data X
#' @param p.y the dimension of data Y
#' @param h the length of observation to predict
#' @return prediction matrix
#'
pred.block <- function(X, phi, j, p.x, p.y, h){
  concat.X <- matrix(0, p.x, h);
  concat.Y <- matrix(0, p.y, h);
  concat.X[, 1:h] <- as.matrix(X[, (j):(j+h-1)]);
  for ( i in 1:h){
    temp <- matrix(0, p.y, 1);
    if(p.y == 1){
      temp <- temp +  t(as.matrix(phi[, (1):(p.x)]))%*%concat.X[, i];
    }else{
      temp <- temp +  as.matrix(phi[, (1):(p.x)])%*%concat.X[, i];
    }

    concat.Y[, i] <- temp;
  }
  return(as.matrix(concat.Y[, 1:h]))
}


#' helper function for detection check
#' @param pts the estimated change points
#' @param brk the true change points
#' @return a vector of timepoints
#'
remove.extra.pts <- function(pts, brk){
  m.hat <- length(brk)-1;
  if(length(pts) <= m.hat){
    return (pts)
  }else{
    pts.temp <- rep(0, m.hat);
    for(i in 1:m.hat){
      origin <- brk[i];
      dis <- rep(0, length(pts));
      for(j in 1:length(pts)){
        dis[j] <- abs(origin - pts[j]);
      }
      ll <- min(which.min(dis));
      pts.temp[i] <- pts[ll];
    }
    pts <- pts.temp;
    return(pts)
  }

}



#' BIC threshold for final parameter estimation
#' @param method method name for the model: Constant: Mean-shift Model; MvLR: Multivariate Linear Regression; MLR: Multiple Linear Regression
#' @param beta.final a combined matrix of estimated parameter coefficient matrices for all stationary segementations
#' @param k dimensions of parameter coefficient matrices
#' @param m.hat number of estimated change points
#' @param brk vector of estimated change points
#' @param data_y input data matrix (response), with each column representing the time series component
#' @param data_x input data matrix (predictor), with each column representing the time series component
#' @param b_n the block size
#' @param nlam number of hyperparameters for grid search
#' @return lambda.val.best, the tuning parameter lambda selected by BIC.
#'
BIC.threshold <- function(method, beta.final, k, m.hat, brk, data_y, data_x = NULL,
                          b_n = 2, nlam = 20){
  # print(b_n)
  brk.full <- c(1, brk)
  jj = 1; flag = 0
  if(!is.null(beta.final)){
    lambda.val.best <- c()
    flag = 1
    # print("jj");print(jj)
    for(i in 1:m.hat){
      temp <- unlist(beta.final[, ((i-1)*k+1):(i*k)] )
      lambda.max <-  max(abs(temp))
      if(lambda.max > 0){
        lambda.min <-  min(abs(temp[temp!=0]))
        # print(lambda.max)
        # print(lambda.min)
        if(lambda.max/lambda.min >= 10^4){
          nlam <- 50
        }
        if(lambda.max/lambda.min >= 10^8){
          lambda.min <-  lambda.max*10^(-4)
        }
        delata.lam <- (log(lambda.max)-log(lambda.min))/(nlam -1)
        lambda.val.full <-  sapply(1:(nlam), function(jjj) lambda.min*exp(delata.lam*(nlam-jjj)))
        # print(lambda.val.full)
        mse.res <- c()
        BIC.res <- c()
        for(j in 1:nlam){
          lambda.val = lambda.val.full[j]
          # print("lambda value:")
          # print(lambda.val)

          # print(beta.temp)
          if(method == 'Constant'){
            beta.temp <- beta.final[((i-1)*k+1):(i*k)]
            beta.temp[abs(beta.temp) < lambda.val] <- 0
            data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
            data_x.temp <- matrix(1, nrow = brk.full[i+1] - brk.full[i])
            data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, nrow = 1)
            residual.temp <- data_y.temp - data_y.est
            BIC.res <- c(BIC.res, BIC(t(residual.temp[b_n : (nrow(residual.temp)-b_n), ]), phi = matrix(beta.temp, ncol = 1) )$BIC )
          }
          if(method == "MLR"){
            beta.temp <- beta.final[((i-1)*k+1):(i*k)]
            beta.temp[abs(beta.temp) < lambda.val] <- 0
            data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
            data_x.temp <- data_x[brk.full[i]:(brk.full[i+1]-1), ]
            data_y.est <- as.matrix(data_x.temp)%*%matrix(beta.temp, ncol = 1)
            residual.temp <- data_y.temp - data_y.est
            # print(residual.temp[1:10,])
            BIC.res <- c(BIC.res, BIC(t(residual.temp[b_n : (nrow(residual.temp)-b_n), ]), phi = matrix(beta.temp, nrow = 1) )$BIC )
          }
          if(method == "MvLR"){
            beta.temp <- beta.final[, ((i-1)*k+1):(i*k)]
            beta.temp[abs(beta.temp) < lambda.val] <- 0
            data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
            data_x.temp <- data_x[brk.full[i]:(brk.full[i+1]-1), ]
            # print(dim(data_x.temp))
            # print(dim(t(beta.temp)))
            data_y.est <- as.matrix(data_x.temp)%*%t(as.matrix(beta.temp))
            # print(dim(data_y.est))
            # print(dim(data_y.temp))
            residual.temp <- data_y.temp - data_y.est
            # print(residual.temp[1:10,])
            BIC.res <- c(BIC.res, BIC(t(residual.temp[b_n : (nrow(residual.temp)-b_n),]), phi = as.matrix(beta.temp))$BIC )
          }
          if(method == "VAR"){
            # k <- p*q
            p <- dim(data_y)[2]
            q <- dim(beta.final)[2]/(m.hat*p)

            beta.temp <- beta.final[, ((i-1)*k+1):(i*k)]
            beta.temp[abs(beta.temp) < lambda.val] <- 0


            data_y.temp <- as.matrix(data_y[(brk.full[i]+q ):(brk.full[i+1]-1), ] )
            data_x.temp <- data_y[brk.full[i]:(brk.full[i+1]-1-q+1), ]

            data_y.est <- matrix(0, dim(data_y.temp)[1], dim(data_y.temp)[2])
            # print(dim(data_x.temp))
            h = brk.full[i+1] - brk.full[i] - q
            for ( j in 1:h){
              temp <- matrix(0, p, 1);
              for (qq in 1:q){
                temp <- temp + beta.temp[, ((qq-1)*p+1):(qq*p)]%*%data_x.temp[j + q - qq, ];
              }
              data_y.est[j, ] <- t(temp)
            }






            # data_y.est <- as.matrix(data_x.temp)%*%t(as.matrix(beta.temp))
            # print(dim(data_y.est))
            # print(dim(data_y.temp))
            residual.temp <- data_y.temp - data_y.est
            BIC.res <- c(BIC.res, BIC(t(residual.temp[b_n : (nrow(residual.temp)-b_n),]), phi = as.matrix(beta.temp), method = "VAR")$BIC )
          }
        }
      }else{
        # print("not good, continue!")
        lambda.min <- 0
        lambda.val.full <- c(0)
        # print(lambda.val.full)
        BIC.res <- c(0)
        flag = 0
      }
      # print("BIC results:")
      # print(BIC.res)
      # if(which.min(BIC.res) == 1){
      #   print("not good, continue!")
      #   flag = 0
      #   break
      # }else{
      #   lambda.val.best <- c(lambda.val.best, lambda.val.full[which.min(BIC.res)])
      # }
      lambda.val.best <- c(lambda.val.best, lambda.val.full[which.min(BIC.res)])
    }
  }
  return(lambda.val.best)
}

#' BIC threshold for final parameter estimation (GGM)
#' @param beta.final a combined matrix of estimated parameter coefficient matrices for all stationary segementations
#' @param k dimensions of parameter coefficient matrices
#' @param m.hat number of estimated change points
#' @param brk vector of estimated change points
#' @param data_y input data matrix (response), with each column representing the time series component
#' @param data_x input data matrix (predictor), with each column representing the time series component
#' @param b_n the block size
#' @param nlam number of hyperparameters for grid search
#' @return lambda.val.best, the tuning parameter lambda selected by BIC.
#'
BIC.threshold.ggm <- function(beta.final, k, m.hat, brk, data_y, data_x = NULL, b_n = 2 , nlam = 20){
  brk.full <- c(1, brk)
  flag = 0
  p.y <- k+1
  BIC.all <- vector("list", m.hat);
  if(!is.null(beta.final)){
    lambda.val.best <- c()
    flag = 1
    for(i in 1:m.hat){
      temp <- unlist(beta.final[, ((i-1)*k+1):(i*k)] )
      # lambda.max <-  max(abs(unlist(beta.final[, ((i-1)*k+1):(i*k)] )))
      lambda.max <-  max(abs(temp))
      if(lambda.max > 0){
        lambda.min <-  min(abs(temp[temp!=0]))
        # lambda.min <-  lambda.max*10^(-2)
        print(lambda.max)
        print(lambda.min)
        if(lambda.max/lambda.min >= 10^4){
          nlam <- 50
        }
        if(lambda.max/lambda.min >= 10^4){
          lambda.min <-  lambda.max*10^(-4)
        }
        delata.lam <- (log(lambda.max)-log(lambda.min))/(nlam -1)
        lambda.val.full <-  sapply(1:(nlam), function(jjj) lambda.min*exp(delata.lam*(nlam-jjj)))
        print(lambda.val.full)
        mse.res <- c()
        BIC.res <- c()
        for(j in 1:nlam){
          lambda.val = lambda.val.full[j]
          # print(lambda.val)
          beta.temp <- beta.final[, ((i-1)*k+1):(i*k)]
          beta.temp[abs(beta.temp) < lambda.val] <- 0
          data_y.temp <- as.matrix(data_y[brk.full[i]:(brk.full[i+1]-1), ] )
          data_y.est <- matrix(0, ncol = ncol(data_y.temp), nrow = nrow(data_y.temp))
          for(j.1 in 1:p.y){
            data.x.temp <- data_y.temp[, -j.1]
            # data_y.est[, j.1] <- pred.block(t(data.x.temp), matrix(beta.temp[j.1, ], ncol = k),
            #                                                1, k, 1, nrow(data_y.temp));
            data_y.est[, j.1] <- data.x.temp%*% matrix(beta.temp[j.1, ], ncol= 1)
          }
          residual.temp <- data_y.temp - data_y.est
          # if(drop ==TRUE){
          #   residual.temp <- residual.temp[b_n:(ncol(residual.temp)-b_n) ,]
          # }
          BIS.sum.all <- 0
          for(j.1 in 1: p.y){
            BIS.sum.all <- BIS.sum.all + BIC(t(residual.temp[b_n : (nrow(residual.temp)-b_n), j.1]), phi = matrix(beta.temp[j.1, ], nrow = 1) )$BIC
          }
          BIC.res <- c(BIC.res, BIS.sum.all )
        }
      }else{
        print("not good, continue!")
        lambda.min <- 0
        lambda.val.full <- c(0)
        print(lambda.val.full)
        BIC.res <- c(0)
        flag = 0
      }
      if(which.min(BIC.res) == 1){
        # print("not good, try larger range of lambda?!")
        flag = 0
        lambda.val.best <- c(lambda.val.best, lambda.val.full[which.min(BIC.res)])
      }else{
        lambda.val.best <- c(lambda.val.best, lambda.val.full[which.min(BIC.res)])
      }
      BIC.all[[i]] <- BIC.res
    }

  }

  return(list(lambda.val.best = lambda.val.best))
}
