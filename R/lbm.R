#' @name lbm
#' @aliases lbm
#' @title Log binomial regression model in Exact method
#' @import stats utils
#' @keywords log binomial model boundary risk ratio/relative risk convergence
#'
#' @description If the maximum likelihood (ML) solution lies on the boundary of the
#'              parameter space in log binomial model, a special method is needed
#'              since the standard fitting algorithm may meet numerical difficulties.
#'              Exact method can overcome the difficulties and address the ML solution
#'              when it lies on the boundary of the parameter space.\code{lbm}
#'              implemented the exact method to address the ML solution in the log binomial model.
#'
#' @usage lbm(formula, data,contrasts = NULL,subset,na.action,lfv=0.95,
#'            vce = "oim",rescode=NULL,control=lbm.control(),\dots)
#'
#' @param formula an object of class "\code{\link[stats]{formula}}" (or one that
#'     can be coerced to that class): a symbolic description of the model
#'     to be fitted. The details of model specification are given under \sQuote{Details}.
#' @param data an optional data frame, list or environment (or object
#'     coercible by \code{as.data.frame} to a data frame) containing the variables
#'     in the model. If not found in \code{data}, the variables are taken from
#'     \code{environment(formula)}, typically the environment from which \code{lbm}
#'     is called.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#'     \code{model.matrix.default}.
#' @param subset a specification of the rows to be used: defaults to all rows.
#'     This can be any valid indexing vector (see [\code{.data.frame}] for the
#'     rows of data or if that is not supplied, a data frame made up of the
#'     variables used in formula.
#' @param na.action a function which indicates what should happen when the data
#'     contain \code{NA}s. The default is set by the \code{na.action} setting of
#'     \code{\link[stats]{na.action}}, and is \code{\link[stats]{na.fail}} if that
#'     is unset. The 'factory-fresh' default is \code{\link[stats]{na.omit}}.
#'     Another possible value is \code{NULL}, no action. Value
#'     \code{\link[stats]{na.exclude}} can be useful.
#' @param lfv a testing range option which decides the range of boundary vector candidates included for testing.
#'     The default value is 0.95, which means the covariate vectors with probability
#'     greater than 0.95 will be included into boundary pairing system as boundary vector candidates.
#' @param vce the type of the information matrix used to attain the variance-covariance matrix.
#'     Two options could be selected, observed information matrix (OIM) and expected information
#'     matrix (EIM). The default \code{vce} is \code{"OIM"}. This argument only works in the data with
#'     boundary vector. If there is no boundary vector included in the data, the results
#'     are from \code{glm}. In the \code{glm}, the standard error is calculated
#'     by expected information matrix.
#' @param rescode is an option to code the response variable if it is a factor.
#' @param control The \code{control} argument of \code{lbm} is by default passed to the arguments of \code{\link{lbm.control}}.
#' @param \dots For \code{lbm}, arguments to be used to form the default control argument if it is not supplied directly.
#'
#' @details A typical predictor has the form \code{response ~ terms} where
#'     response is the (numeric) response vector and terms is a series of
#'     terms which specifies a linear predictor for response. A terms
#'     specification of the form \code{first + second} indicates all the
#'     terms in first together with all the terms in second with any
#'     duplicates removed.

#'     A specification of the form \code{first:second} indicates the set of
#'     terms obtained by taking the interactions of all terms in first with
#'     all terms in second. The specification \code{first*second} indicates
#'     the cross of first and second. This is the same as
#'     \code{first + second + first:second}.

#'     The terms in the formula will be re-ordered so that main effects come
#'     first, followed by the interactions, all second-order, all third-order
#'     and so on: to avoid this pass a terms object as the formula.
#'
#' @return
#' \code{lbm} returns an object of class inheriting from \code{"lbm"} which
#' inherits from the class \code{"lbm"}. The function \code{summary} (i.e.,
#' \code{summary.lbm}) can be used to obtain or print a summary of the estimates
#' and the relevant confidence interval. The argument \code{CF.lvl} in
#' \code{summary} represents the level of confidence interval claimed in the
#' model. The default value is \code{CF.lvl=0.95}. Optionally, Risk ratio estimates
#' and their related confidence interval are offered as an argument \code{RR} in
#' the \code{summary}. The default \code{RR=FALSE} is not to display them.
#'
#' An object of class \code{"lbm"} is a list containing at least the following
#' components:
#'
#' \item{coefficients}{a named vector of coefficients}
#' \item{residuals}{the working residuals, that is the residuals in the final
#' iteration of the IWLS fit.}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the
#' linear predictors by the inverse of the log link function.}
#' \item{linear.predictors}{the linear fit on log scale.}
#' \item{deviance}{twice the absolute value of maximized log-likelihood.}
#' \item{aic}{A version of Akaike's An Information Criterion, minus twice the
#' maximized log-likelihood plus twice the number of parameters, computed by
#' the \code{aic} component of the family. For the binomial model, the dispersion
#' is fixed at one and the number of parameters is the number of coefficients.}
#' \item{null.deviance}{The deviance for the null model, comparable with
#' \code{deviance}. The null model will only include an intercept if there is
#' one in the model.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{df.null}{the residual degrees of freedom for the null model.}
#' \item{response}{the response vector used in the mode.l}
#' \item{vcov}{the unscaled (\code{dispersion = 1}) estimated covariance matrix
#' of the estimated coefficients.}
#' \item{vce}{the type of information matrix applied.}
#' \item{call}{the matched call.}
#' \item{na.action}{(where relevant) information returned by \code{stats::model.frame}
#' on the special handling of \code{NA}.}
#' \item{contrasts}{(where relevant) the contrasts used.}
#' \item{formula}{the formula supplied.}
#' \item{factor}{the order of factors used in the response variable.}
#' \item{bvector}{the matrix of boundary vectors.}
#' \item{bv}{logical. Determines whether the model has boundary vectors.}
#'
#' @references
#'     Petersen, M. R. & Deddens, J. A. (2010). Maximum likelihood estimation
#'     of the log-binomial model. \emph{Communications in Statistics - Theory
#'     and Methods}, 39: 5, 874 - 883.
#'
#' @seealso
#' \code{glm}, \code{lm}.
#'
#' @examples
#' ## Two examples are from Petersen, M. R. & Deddens, J. A. (2010).
#'
#' ## Example 1.
#' x<-c(1:10)
#' y<-c(0,0,0,0,1,0,1,1,1,1)
#' data<-data.frame(x,y)
#' a<-lbm(formula=y~x,data=data,vce="eim")
#'
#' ## Example 2.
#' x1<-c(1:11)
#' x2<-x1^2
#' y<-c(10,6,4,3,3,2,3,3,4,6,10)
#' dat<-cbind(x1,x2,y)
#' dat1<-apply(dat, 1, function(t) {
#'   temp<-data.frame(x1=rep(t[1],10),x2=rep(t[2],10),y=0)
#'   temp$y[1:t[3]]<-1
#'   return(temp)
#' })
#' data<-do.call(rbind, dat1)
#' a<-lbm(formula=y~x1+x2,data=data)
#' summary(a)
#'
#' @export lbm

lbm <- function(formula, data,contrasts = NULL,subset,na.action,lfv=0.95,
                vce = "oim",rescode=NULL,control=lbm.control(),...)
{
  scall <- match.call()
  cf<-match.call(expand.dots = FALSE)
  c<-match(c("formula", "data", "subset","na.action"), names(cf), 0L)
  cf <- cf[c(1L, c)]
  cf$drop.unused.levels <- TRUE
  cf[[1L]] <- quote(stats::model.frame)
  data <- eval(cf, parent.frame())
  vce<-substitute(vce)
  if(!is.character(vce)) vce <- deparse(vce)
  #data <- stats::model.frame(formula, data = data)
  na.action <- attr(data, "na.action")
  terms<-attr(data,"terms")
  y<-model.extract(data,"response")
  # names_observed<-attr(data,"names")[1]
  factor_index<-NULL
  if(is.factor(y)) {
    if(is.null(rescode))lvl<-levels(y) else lvl<-as.factor(rescode)
    if(length(lvl)!=2L) stop("Response variable has to be binary.")
    y_temp<-y
    y<-as.numeric(y)
    y[which(y_temp==lvl[1])]<-0
    y[which(y_temp==lvl[2])]<-1
    factor_index<-data.frame(Levels=lvl, Num=c(0,1))
  }
  x<-model.matrix(terms,data,contrasts)
  contrasts<-attr(x,"contrasts")
  start<-lbm.InitValue(x=x,y=y)
  if(any(is.na(start$start)==TRUE))
    stop("Please check the model. Some covariates may be 100 percent explained by others.
        \rIf interaction term is included, please reform the data and generate them and
        \rrerun the function.")
  res<-lbm.NR(x=x,y=y,start=start$start,
              noconstant=F,control=control,...)
  fitvalues<-rep(0,nrow(x))
  indicator<-rep(0,nrow(x))
  x<-cbind(x,indicator=indicator, fitvalues=fitvalues)
  temp_null.deviance <- res$null.deviance
  deviance_temp<-res$deviance
  x[,"indicator"][which(res$fitted.values > lfv & y!=0)]<-1
  x[,"fitvalues"]<-res$fitted.values
  indicator<-x[,"indicator"]
  start<-res$coefficients
  b_vectors <- as.data.frame(x)[which(as.data.frame(x)$indicator==1),]
  # Remove the first column which is intercept
  b_vectors<-subset.data.frame(b_vectors,select=-1)
  x<-subset.matrix(x,select=-c(indicator,fitvalues))
  b_vectors<-b_vectors[order(b_vectors$fitvalues, decreasing = T),]
  colname<-colnames(b_vectors)[!is.na(match(colnames(b_vectors),colnames(x)))]
  b_vectors<-subset.data.frame(b_vectors,select=colname)
  df.null <- nrow(x)-1
  df.residual <- nrow(x) - ncol(x)
  aic <- res$deviance + 2*ncol(x)
  res <- c(res, list(bv=0,aic = aic,vce=vce,
           df.null = df.null,df.residual = df.residual))
  if(nrow(b_vectors)!=0) {
    res.temp <- tryCatch({
      suppressWarnings(lbm.fit(x=x,y=y, bvectors = b_vectors,control=control,
                      num_bvectors = nrow(b_vectors), indicator = indicator,vce=vce,
                      start=start, deviance=deviance_temp,null.dev=temp_null.deviance,...))
      },
      error=function(err) {
        return(ind<-0)
    })
    if(is.list(res.temp))
      res <- c(res.temp, list(vce=vce, bv=1))
  }
  res <- c(res, list(factor=factor_index,formula = formula,
          contrasts=contrasts,na.action=na.action,response=y))
  res$call<-scall
  class(res)<-c("lbm")
  return(res)
}

lbm.fit<-function(x,y, bvectors = NULL,null.dev=NULL,vce,
                  num_bvectors = NULL, indicator = NULL,control,
                  start=NULL, deviance=NULL,...)
{
  num_var<-ncol(x)-1
  if(nrow(bvectors)>=1) {
    names<-row.names(bvectors)
    bvectors_1e<-floor(bvectors*1e+5)
    ind_dup<-cumsum(!duplicated(bvectors_1e))
    bvectors<-cbind(bvectors,ind_dup=ind_dup)
    bvectors_temp <- bvectors[!duplicated(bvectors_1e[,1:(ncol(bvectors_1e)-1)]), ]
    max_num<-min(num_var,nrow(bvectors_temp))
    improved<-0
    for(j in 1:max_num)
    {
      # if(j==5) break
      message("\nSearching for the admissible combinations of",j,"boundary vectors.\n",appendLF=F)
      if(improved==0) ind_matrix<-t(utils::combn(nrow(bvectors_temp),j))
        else ind_matrix<-t(utils::combn(nrow(bvectors_temp)-bv_num,j-bv_num))
      for(i in 1:nrow(ind_matrix))
      {
        message(".",sep="",appendLF=F)
        if(i %% 10==0) message("|",sep="",appendLF=F)
        if(i %% 5==0 & i %% 10!=0) message(" ",sep="",appendLF=F)
        if(i %% 50==0) message("\n",sep="",appendLF=F)
        x_temp<-cbind(x, indicator=indicator)

        if(improved==0) temp_bvectors<-bvectors_temp[ind_matrix[i,],]
          else temp_bvectors<-rbind(bvectors_temp[1:bv_num,],
            bvectors_temp[(bv_num+1):nrow(bvectors_temp),][ind_matrix[i,],])
        x_temp[,"indicator"][!(row.names(x_temp) %in%
            names[which(sapply(bvectors$ind_dup,
            function(t){t %in% temp_bvectors$ind_dup}))])]<-0
        indicator_temp<-x_temp[,"indicator"]
        x_temp<-subset(x_temp, select=-indicator)
        temp_bvectors<-subset(temp_bvectors,select=-ind_dup)
        conv<-tryCatch({
          lbm.reform(x=x_temp,y=y, t = temp_bvectors,
                    num_bvectors = nrow(temp_bvectors),
                    indicator=indicator_temp, start=start,
                    vce=vce,control=control,...)
        },
          error=function(err){
            return(ind<-0)
        })
        #if(conv$conv==1) print(1) else print(0)
        if(is.list(conv)) {
          if(conv$deviance<=deviance) break
        }
        if(i==control$maxit) break
      }
      if(is.list(conv)) {
        if(conv$deviance<=deviance) {
          deviance<-conv$deviance
          res.temp<-conv
          b_vector<-temp_bvectors
          if(nrow(temp_bvectors)!=nrow(bvectors)) {
            temp_fit<-conv$fitted.values[row.names(bvectors_temp)[-which(row.names(bvectors_temp)%in%
                row.names(temp_bvectors))]]
            # Reorder the boundary vectors left in the bvectors_temp by their fitted probabilities.
            temp_fit_name<-names(temp_fit[order(temp_fit,decreasing = T)])
            temp_fit_name<-c(row.names(temp_bvectors),temp_fit_name)
            bvectors_temp<-bvectors_temp[temp_fit_name,]
          }
          # Most recent number of boundary vectors.
          bv_num<-j
          improved<-1
          message("\nMinimum Deviance improved. \n",appendLF=F)
        }
      }
      if(j==max_num & improved==0)
        stop("No admissible pairs of boundary vectors could be found.")
    }
  }
  if(is.list(res.temp)) {
    res.temp<-lbm.makeup(beta=res.temp$beta,
                        vcov=res.temp$vcov,
                        t.repa=res.temp$t.repa)
    res.temp$beta<-res.temp$beta[colnames(x)]
    res.temp$vcov<-res.temp$vcov[colnames(x),colnames(x)]
    eta <- x %*% res.temp$beta
    eta <- as.vector(eta, mode = "numeric")
    mu <- exp(eta)
    ## Working residuals
    residuals <-  (y - mu)/exp(eta)
    deviance<-lbm.NR.dev(y=y,u=mu)
    df.null <- nrow(x) - 1
    df.residual <- nrow(x) - ncol(x)
    aic <- deviance + 2*ncol(x)
    res <- list(coefficients = res.temp$beta,bvector = b_vector,
                linear.predictors = eta, fitted.values = mu,
                null.deviance = null.dev, deviance = deviance,
                df.null = df.null, df.residual = df.residual,
                residuals=residuals,aic = aic,vcov = res.temp$vcov)
    return(res)
  } else {
    stop("No admissible combination of boundary vectors in the model.")
  }
}


#' @name lbm.control
#' @title Auxiliary control for \code{lbm}
#' @description
#' Auxiliary function for \code{lbm} fitting. Only used
#' internally by \code{lbm}.
#'
#' @usage lbm.control(epsilon = 1e-8, maxit = 100)
#' @param epsilon positive convergence tolerance \code{epsilon};
#' @param maxit integer giving the maximal number of iterations.
#' @seealso \code{glm.control}
#' @return A list with components named as the arguments.
#' @export lbm.control
lbm.control <- function(epsilon = 1e-8, maxit = 100)
{
  if(!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if(!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit)
}

lbm.makeup<-function(beta,vcov,t.repa)
{
  t<-cbind(0,t.repa)
  name.vcov<-c("(Intercept)",names(t.repa))
  beta.complete<-rep(0,length(name.vcov))
  names(t)<-names(beta.complete)<-name.vcov
  vcov.complete<-matrix(0,nrow=length(name.vcov),
                        ncol=length(name.vcov),
                        dimnames=list(name.vcov,name.vcov))
  vcov.complete[row.names(vcov),colnames(vcov)]<-as.matrix(vcov)
  beta.complete[names(beta)]<-beta
  # Makeup estimates and variance-covariance matrix.
  for (i in nrow(t.repa):1)
  {
    cov.vector<-apply(vcov.complete,1,function(x) {
                  cov<--sum(x*t[i,])
                  return(cov)
                })
    var.temp<-sum(diag(t[i,])%*%vcov.complete%*%t(t[i,]))
    vcov.complete[i,]<-vcov.complete[,i]<-cov.vector
    vcov.complete[i,i]<-var.temp
    beta.complete[i]<--sum(t[i,]*beta.complete)
  }
  return(list(beta=beta.complete,
        vcov=vcov.complete))
}

#' @method summary lbm
#' @export
summary.lbm<-function(object, CF.lvl=0.95, RR=FALSE,...)
{
  if(0.5>=CF.lvl || CF.lvl>=1)
    stop("CF.lvl should be a value between 0.5 and 1.")
  alpha<-1-CF.lvl
  var.cf <- diag(object$vcov)
  coef<-object$coefficients
  ## calculate coef table
  s.err <- sqrt(var.cf)
  zvalue <- coef/s.err
  if(object$bv==1) {
    if(tolower(object$vce)=="oim"){
      dn<-c("Estimate", "Std.Err")
    } else {
      dn<-c("Estimate", "Std.Err")
    }
  } else {
    dn<-c("Estimate", "Std.Err")
  }
  pvalue <- 2*pnorm(-abs(zvalue))
  coef.table <- cbind(coef, s.err, zvalue, pvalue)
  confint<-apply(coef.table[,1:2], 1, function(t){
    lower<-t[1]+qnorm(alpha/2)*t[2]
    upper<-t[1]+qnorm(1-alpha/2)*t[2]
    dat<-cbind(lower,upper)
  })
  coef.table<-cbind(coef.table,t(confint))
  ## Borrow idea from the Stata to present confidence interval.
  confint_name<-c(paste("[",(1-alpha)*100,"% ","Conf.",sep=""),"Interval]")
  dimnames(coef.table)<-list(names(object$coefficients),
    c(dn, "z value","Pr(>|z|)",confint_name))
  object$coefficients<-coef.table
  if(RR==T) {
    RR.table<-exp(cbind(coef,t(confint)))
    dimnames(RR.table)<-list(row.names(object$coefficients),
      c("Risk Ratio", confint_name))
    object$RR<-RR.table
  }
  ans<-object
  class(ans) <- "summary.lbm"
  return(ans)
}

printCoefmat.lbm<-function (x, digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),signif.legend = signif.stars,
  dig.tst = max(1L, min(5L, digits-1L)), cs.ind = 1:k, tst.ind = k + 1,
  zap.ind = integer(),  P.values = NULL, has.Pvalue = nc >= 4L &&
    length(cn <- colnames(x)) && substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"),
  eps.Pvalue = .Machine$double.eps, na.print = "NA", quote = FALSE, right = TRUE, ...)
{
  if (is.null(d <- dim(x)) || length(d) != 2L)
    stop("'x' must be coefficient matrix/data frame")
  nc <- d[2L]-2
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  } else if (P.values && !has.Pvalue)
      stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  } else xm <- data.matrix(x)
  k <- nc - has.Pvalue - (if (missing(tst.ind))
    1
    else length(tst.ind))
  if (!missing(cs.ind) && length(cs.ind) > k)
    stop("wrong k / cs.ind")
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], digits)
  if (length(cs.ind)) {
    acs <- abs(coef.se <- xm[, cs.ind, drop = FALSE])
    if (any(ia <- is.finite(acs))) {
      digmin <- 1 + if (length(acs <- acs[ia & acs != 0]))
        floor(log10(range(acs[acs != 0], finite = TRUE)))
      else 0
      Cf[, cs.ind] <- format(round(coef.se, max(1L, digits -
          digmin)), digits = digits)
    }
  }
  if (length(tst.ind))
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
      digits = digits)
  if (any(r.ind <- !((1L:d[2L]) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc))))
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue)
    ok[, -nc]
  else ok
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".")
    x1 <- chartr(dec, ".", x1)
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(1L,
      digits - 1L))
  }
  if (any(ina))
    Cf[ina] <- na.print
  ## print(Cf)
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP], digits = dig.tst,
        eps = eps.Pvalue)
      ## signif.stars <- signif.stars && any(pv[okP] < 0.1)
      ## print(signif.stars)
      if (signif.stars) {
        Signif <- symnum(pv, corr = FALSE, na = FALSE,
          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
          symbols = c("***", "**", "*", ".", " "))
        Cf <- cbind(Cf, format(Signif))
        ## print(Cf)
      }
    }
    else signif.stars <- FALSE
  }
  else signif.stars <- FALSE

  print.default(Cf, quote = quote, right = right, na.print = na.print,
    ...)
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(Signif,
      "legend")))
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
        4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}

#' @method print summary.lbm
#' @export
print.summary.lbm <- function (x, digits = max(6L, getOption("digits")-3L),
                              signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n")
  dput(x$call)
  cat("\nFormula:\n", paste(deparse(x$formula), sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nBoundary Vectors:\n")
  print(x$bvector, digits=digits)
  if(!is.null(x$factor)) {
    cat("\nIndex of response variable:\n")
    factor_names<-as.character(x$factor[,1])
    factor_names<-format(factor_names,
      width=max(nchar(factor_names)),
      justify="right")
    cat(factor_names[1],"$ ",x$factor[,2][1],"\n",
      factor_names[2],"$ ",x$factor[,2][2],"\n",sep="")
  }
  cat("\nCoefficients:\n")
  coefs <- x$coefficients
  printCoefmat.lbm(coefs, digits=min(4L,digits-2L),
    signif.stars = signif.stars, na.print = "NA")
  if(!is.null(x$RR)) {
    cat("\nRisk Ratio Estimates:\n")
    print(x$RR, digits=digits)
  }
  cat("\n", apply(cbind(paste(format(c("Null","Residual"), justify="right"),
    "deviance:"), format(unlist(x[c("null.deviance","deviance")]),
      digits = max(5L, digits + 3L)), " on",
    format(unlist(x[c("df.null","df.residual")])),
    " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(5L, digits + 3L)), "\n", sep = "")
  invisible(x)
}

#' @method print lbm
#' @export
print.lbm <- function(x, digits = max(3L, getOption("digits")-1L), ...)
{
  cat("\nCall:\n")
  dput(x$call)
  cat("\nFormula:\n", paste(deparse(x$formula),
    sep = "\n", collapse = "\n"), "\n", sep = "")
  cat("\nBoundary Vectors:\n")
  print(x$bvector, digits=digits-2L)
  if(!is.null(x$factor)) {
    cat("\nIndex of response variable:\n")
    factor_names<-as.character(x$factor[,1])
    factor_names<-format(factor_names,
      width=max(nchar(factor_names)),
      justify="right")
    cat(factor_names[1],"$ ",x$factor[,2][1],"\n",
      factor_names[2],"$ ",x$factor[,2][2],"\n",sep="")
  }
  if(length(coef(x))) {
    cat("\nCoefficients")
    if(is.character(co <- x$contrasts))
      cat("  [contrasts: ",
        apply(cbind(names(co),co), 1L, paste, collapse = "="), "]")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
      print.gap=2, quote = FALSE)
  } else cat("No coefficients\n\n")
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
    x$df.residual, "Residual\n")
  cat("     Null Deviance:",	format(signif(x$null.deviance,
    digits=max(5L, digits + 3L))),
    "\n Residual Deviance:", format(signif(x$deviance,
      digits=max(5L, digits + 3L))),
    "\tAIC:", format(signif(x$aic, digits=max(4L, digits + 3L))),"\n")
  invisible(x)
}
#' @method vcov lbm
#' @export
vcov.lbm <- function(object, ...) object$vcov

#' @method vcov summary.lbm
#' @export
vcov.summary.lbm <- function(object, ...)  object$vcov

#' @method logLik lbm
#' @export
logLik.lbm <- function(object, ...)  -object$deviance/2

#' @method logLik summary.lbm
#' @export
logLik.summary.lbm <- function(object, ...)  -object$deviance/2
