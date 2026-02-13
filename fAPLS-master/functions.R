if (!("expm" %in% rownames(installed.packages())))
  install.packages("expm")
if (!("fda" %in% rownames(installed.packages())))
  install.packages("fda")
if (!("FRegSigCom" %in% rownames(installed.packages())))
  install.packages("FRegSigCom")
if (!("latex2exp" %in% rownames(installed.packages())))
  install.packages("latex2exp")
if (!("MASS" %in% rownames(installed.packages())))
  install.packages("MASS")
if (!("plsdepot" %in% rownames(installed.packages())))
  install.packages("plsdepot")
if (!("plsgenomics" %in% rownames(installed.packages())))
  install.packages("plsgenomics")
library(expm)
library(fda)
library(FRegSigCom)
library(plsdepot)
library(plsgenomics)


#'
#'
#' type =
#' 100: \int f(t)dt
#' 201: \int f(s,t)dt
#' 211: \int f(s,t)g(t)dt
#' 222: \int f(s,w)g(w,t)dw

integral = function(f, g = NULL, domain, type){
  if (type == 100){
    f = as.vector(f)
    len.f = length(f)
    result = sum((f[-1] + f[-len.f]) * diff(domain))/2
    return(result)
  }

  if (type == 201){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = matrix((rowSums(f[, -1] * gap.mat) + rowSums(f[, -ncol.f] * gap.mat))/2, nrow = nrow.f, ncol = 1)
    return(result)
  }

  if (type == 211){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = ((f[, -1] * gap.mat) %*% as.matrix(g[-1]) + (f[, -ncol.f] * gap.mat) %*% as.matrix(g[-ncol.f]))/2
    return(as.matrix(result))
  }

  if (type == 222){
    nrow.f = dim(f)[1]
    ncol.f = dim(f)[2]
    gap.mat = matrix(diff(domain), nrow = nrow.f, ncol = length(domain)-1, byrow = T)
    result = ((f[, -1] * gap.mat) %*% g[-1, ] + (f[, -ncol.f] * gap.mat) %*% g[-ncol.f, ])/2
    return(result)
  }
}

# output basis (as a list) of a p-dimensional Krylov subspace
KS = function(p, r.xx, r.xy, domain.x){
  basis = array(NA, dim = c(dim(r.xy), p))
  for (i in 1:p){
    if (i==1){
      basis[, , 1] = r.xy
    }else
      basis[, , i] = integral(r.xx, basis[, , i-1], domain.x, type = 222)
  }
  return(basis)
}

# Gram-Schmidt orthonormalization w.r.t. ker
# basis.origi is a list
GSortho = function(basis.origi, ker, domain.x, domain.y = NULL){

  basis.ortho = array(NA, dim = dim(basis.origi))

  if (is.null(domain.y)){
    p = dim(basis.origi)[2]
    for (i in 1:p){
      psi.tmp = as.matrix(basis.origi[, i])
      if (i > 1){
        for (j in 1:(i-1)){
          integ1 = integral(ker, psi.tmp, domain.x, type = 211)
          integ2 = integral(integ1 * as.matrix(basis.ortho[, j]), domain = domain.x, type = 100)
          psi.tmp = psi.tmp - as.matrix(basis.ortho[, j]) * integ2
        }
      }
      integ1 = integral(ker, psi.tmp, domain.x, type = 211)
      integ2 = integral(integ1 * psi.tmp, domain = domain.x, type = 100)
      basis.ortho[, i] = psi.tmp / integ2^.5
    }
  }

  if (!is.null(domain.y)){
    p = dim(basis.origi)[3]
    for (i in 1:p){
      psi.tmp = basis.origi[, , i]
      if (i > 1){
        for (j in 1:(i-1)){
          integ1 = integral(ker, psi.tmp, domain.x, type = 222)
          integ2 = integral(integ1 * basis.ortho[, , j], domain = domain.y, type = 201)
          integ3 = integral(integ2, domain = domain.x, type = 100)
          psi.tmp = psi.tmp - basis.ortho[, , j] * integ3
        }
      }
      integ1 = integral(ker, psi.tmp, domain.x, type = 222)
      integ2 = integral(integ1 * psi.tmp, domain = domain.y, type = 201)
      integ3 = integral(integ2, domain = domain.x, type = 100)
      basis.ortho[, , i] = psi.tmp / integ3^.5
    }
  }

  return(basis.ortho)
}

test.ortho = function(basis.ortho, ker, domain.x, domain.y = NULL){
  if (is.null(domain.y) & is.null(ker)){
    p = ncol(basis.ortho)
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        H[i, j] = integral(basis.ortho[, i] * basis.ortho[, j], domain = domain.x, type = 100)
      }
    }
  }

  if (is.null(domain.y) & !is.null(ker)){
    p = ncol(basis.ortho)
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ = integral(ker, basis.ortho[, j], domain.x, type = 211)
        H[i, j] = integral(basis.ortho[, i] * integ, domain = domain.x, type = 100)
      }
    }
  }

  if (!is.null(domain.y) & is.null(ker)){
    p = dim(basis.ortho)[3]
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ = integral(basis.ortho[, , i] * basis.ortho[, , j], domain = domain.y, type = 201)
        H[i, j] = integral(integ, domain = domain.x, type = 100)
      }
    }
  }

  if (!is.null(domain.y) & !is.null(ker)){
    p = dim(basis.ortho)[3]
    H = array(NA, dim = c(p, p))
    for (i in 1:p){
      for (j in 1:p){
        integ1 = integral(ker, basis.ortho[, , j], domain.x, type = 222)
        integ2 = integral(basis.ortho[, , i] * integ1, domain = domain.y, type = 201)
        H[i, j] = integral(integ2, domain = domain.x, type = 100)
      }
    }
  }

  return(H)
}

# fAPLS
fAPLS = function(x.old, y.old, domain.x, domain.y, x.new = NULL, y.new = NULL, ridge.par, p.max, tune = 'CV'){

  if (length(ridge.par) > 1)
    hyper.par = cbind(ridge.par, sample.int(p.max, size = length(ridge.par), replace = T))
  else
    hyper.par = cbind(ridge.par, p.max)

  if (!is.null(x.new)){
    mu.x = colMeans(rbind(x.old, x.new))
  }else{
    mu.x = colMeans(x.old)
  }
  mu.y = colMeans(y.old)
  n = nrow(x.old)
  r.xx = cov(x.old)
  r.xy = cov(x.old, y.old)

  if (tune == 'GCV'){
    for (i in 1:length(ridge.par)){
      basis.origi = KS(p.max, r.xx + ridge.par[i] * diag(ncol(r.xx)), r.xy, domain.x)
      basis.ortho = GSortho(basis.origi, r.xx, domain.x, domain.y)
      gamma.vec = numeric(hyper.par[i, 2])
      for (p in 1:hyper.par[i, 2]){
        integ = integral(r.xy * basis.ortho[, , p], domain = domain.y, type = 201)
        gamma.vec[p] = integral(integ, domain = domain.x, type = 100)
        if (p == 1){
          Beta.tilde.curr = gamma.vec[p] * basis.ortho[, , p]
          y.old.tilde.curr = matrix(mu.y, nrow = n, ncol = length(mu.y), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.tilde.curr, domain.x, type = 222)
          Beta.tilde.opti = Beta.tilde.curr
          y.old.tilde.opti = y.old.tilde.curr
          p.opti = p
          lambda.opti = ridge.par[i]
        }else{
          Beta.tilde.curr = Beta.tilde.curr + gamma.vec[p] * basis.ortho[, , p]
          y.old.tilde.curr = matrix(mu.y, nrow = n, ncol = length(mu.y), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.tilde.curr, domain.x, type = 222)
          if (
            sum(integral((y.old.tilde.curr - y.old)^2, domain = domain.y, type = 201))/(n - p - 1)^2 <
            sum(integral((y.old.tilde.opti - y.old)^2, domain = domain.y, type = 201))/(n - p.opti - 1)^2
          ){
            Beta.tilde.opti = Beta.tilde.curr
            y.old.tilde.opti = y.old.tilde.curr
            p.opti = p
            lambda.opti = ridge.par[i]
          }
        }
      }
    }
  }

  if (tune == 'CV'){
    nfold = 5
    cv.array = array(NA, dim = c(p.max, length(ridge.par), nfold))
    groups <- split(sample(1:n), rep(1:nfold, length = n))
    for (k in 1:nfold){
      x.train = x.old[-groups[[k]], ]
      y.train = y.old[-groups[[k]], ]
      x.vali = x.old[groups[[k]], ]
      y.vali = y.old[groups[[k]], ]
      mu.x.train = colMeans(x.train)
      mu.y.train = colMeans(y.train)

      r.xx.train = cov(x.train)
      r.xy.train = cov(x.train, y.train)

      for (i in 1:length(ridge.par)){
        basis.origi = KS(p.max, r.xx.train + ridge.par[i]*diag(ncol(r.xx.train)), r.xy.train, domain.x)
        basis.ortho = GSortho(basis.origi, r.xx.train, domain.x, domain.y)
        gamma.vec = numeric(hyper.par[i, 2])
        for (p in 1:hyper.par[i, 2]){
          integ = integral(r.xy.train * basis.ortho[, , p], domain = domain.y, type = 201)
          gamma.vec[p] = integral(integ, domain = domain.x, type = 100)
          if (p == 1){
            Beta.tilde.curr = gamma.vec[p] * basis.ortho[, , p]
          }else{
            Beta.tilde.curr = Beta.tilde.curr + gamma.vec[p] * basis.ortho[, , p]
          }
          y.vali.tilde.curr = matrix(mu.y.train, nrow = nrow(x.vali), ncol = length(mu.y.train), byrow = T) +
            integral(sweep(x.vali, 2, mu.x.train), Beta.tilde.curr, domain.x, type = 222)
          cv.array[p, i, k] =
            sum(integral((y.vali.tilde.curr - y.vali)^2, domain = domain.y, type = 201))/
            sum(integral(sweep(y.vali, 2, mu.y.train)^2, domain = domain.y, type = 201))
          # cv.array[p, i, k] =
          #   mean(integral((y.vali.tilde.curr - y.vali)^2, domain = domain.y, type = 201))
        }
      }
    }
    cv.mat = apply(cv.array, c(1, 2), mean)
    loc.opt = which(cv.mat == min(cv.mat, na.rm = T), arr.ind = T)
    lambda.opti = ridge.par[loc.opt[1, 2]]
    p.opti = loc.opt[1, 1]

    basis.origi = KS(p.max, r.xx + lambda.opti*diag(ncol(r.xx)), r.xy, domain.x)
    basis.ortho = GSortho(basis.origi, r.xx, domain.x, domain.y)
    gamma.vec = numeric(p.opti)
    for (p in 1:p.opti){
      integ = integral(r.xy * basis.ortho[, , p], domain = domain.y, type = 201)
      gamma.vec[p] = integral(integ, domain = domain.x, type = 100)
      if (p == 1){
        Beta.tilde.opti = gamma.vec[p] * basis.ortho[, , p]
      }else{
        Beta.tilde.opti = Beta.tilde.opti + gamma.vec[p] * basis.ortho[, , p]
      }
    }
  }

  y.new.tilde.opti = NULL
  reispe.opti = NULL
  if (!is.null(x.new)){
    y.new.tilde.opti = matrix(mu.y, nrow = nrow(x.new), ncol = length(mu.y), byrow = T) +
      integral(sweep(x.new, 2, mu.x), Beta.tilde.opti, domain.x, type = 222)
    if (!is.null(y.new))
      reispe.opti =
        sum(integral((y.new.tilde.opti - y.new)^2, domain = domain.y, type = 201))/
        sum(integral(sweep(y.new, 2, mu.y)^2, domain = domain.y, type = 201))
  }
  return(list(Beta = Beta.tilde.opti, lambda = lambda.opti, ncomp = p.opti, y.pred = y.new.tilde.opti, reispe = reispe.opti))
}

# FPC regression
###Dongjiao added K to present the number of basis
#The original function should be FPCR = function(x.old, y.old, domain.x, domain.y, x.new = NULL, y.new = NULL, FVE.x, FVE.y, tune = T)
FPCR = function(x.old, y.old, domain.x, domain.y, x.new = NULL, y.new = NULL, FVE.x, FVE.y, tune = T, K){
  #expansion.x = expand.basis(x.old, domain.x, basis.name = "bspline")#Dongjiao comment this
  expansion.x = expand.basis(x.old, domain.x, basis.name = "bspline", K)#Dongjiao add this
  x.pca.obj = pca.fd(expansion.x$fdObj, nharm = min(dim(x.old) - 1))
  p.max.x = which.max(cumsum(x.pca.obj$varprop) >= FVE.x)
  eigen.fun.x = t(eval.basis(evalarg = domain.x, basisobj = x.pca.obj$harmonics$basis) %*% x.pca.obj$harmonics$coefs[, 1:p.max.x])

  #expansion.y = expand.basis(y.old, domain.y, basis.name = "bspline")#Dongjiao comment this
  expansion.y = expand.basis(y.old, domain.y, basis.name = "bspline", K)#Dongjiao add this
  y.pca.obj = pca.fd(expansion.y$fdObj, nharm = min(dim(y.old) - 1))
  p.max.y = which.max(cumsum(y.pca.obj$varprop) >= FVE.y)
  eigen.fun.y = t(eval.basis(evalarg = domain.y, basisobj = y.pca.obj$harmonics$basis) %*% y.pca.obj$harmonics$coefs[, 1:p.max.y])

  if (!is.null(x.new)){
    mu.x = colMeans(rbind(x.old, x.new))
  }else{
    mu.x = colMeans(x.old)
  }
  mu.y = colMeans(y.old)

  r.xx = cov(x.old)
  r.xy = cov(x.old, y.old)
  n = nrow(x.old)

  Coef.all = array(NA, dim = c(dim(r.xy), p.max.x, p.max.y))
  for (p.x in 1:p.max.x){
    for (p.y in 1:p.max.y){
      integ = integral(r.xy, eigen.fun.y[p.y, ], domain.y, type = 211)
      Coef.all[, , p.x, p.y] = integral(eigen.fun.x[p.x, ] * integ, domain = domain.x, type = 100)/
        x.pca.obj$values[p.x] *
        as.matrix(eigen.fun.x[p.x, ]) %*% t(as.matrix(eigen.fun.y[p.y, ]))
      if (tune == F & p.x == p.max.x & p.y == p.max.y){
        Beta.hat.opti = rowSums(Coef.all, dims = 2)
        p.opti.x = p.max.x
        p.opti.y = p.max.y
      }
      if (tune == T){
        if (p.x == 1 & p.y == 1){
          Beta.hat.opti = Coef.all[, , 1, 1]
          y.old.hat.opti = matrix(mu.y, ncol = length(mu.y), nrow = nrow(x.old), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.hat.opti, domain.x, type = 222)
          p.opti.x = p.x
          p.opti.y = p.y
        }else{
          Beta.hat.curr = rowSums(Coef.all[, , 1:p.x, 1:p.y], dims = 2)
          y.old.hat.curr = matrix(mu.y, ncol = length(mu.y), nrow = nrow(x.old), byrow = T) +
            integral(sweep(x.old, 2, mu.x), Beta.hat.curr, domain.x, type = 222)
          if (
            sum(integral((y.old.hat.curr - y.old)^2, domain = domain.y, type = 201))/(n - max(p.x, p.y) - 1)^2 <
            sum(integral((y.old.hat.opti - y.old)^2, domain = domain.y, type = 201))/(n - max(p.opti.x, p.opti.y) - 1)^2
          ){
            p.opti.x = p.x
            p.opti.y = p.y
            Beta.hat.opti = Beta.hat.curr
            y.old.hat.opti = y.old.hat.curr
          }
        }
      }
    }
  }

  y.new.hat.opti = NULL
  reispe.opti = NULL
  if (!is.null(x.new)){
    y.new.hat.opti = matrix(mu.y, ncol = length(mu.y), nrow = nrow(x.new), byrow = T) +
      integral(sweep(x.new, 2, mu.x), Beta.hat.opti, domain.x, type = 222)
    if (!is.null(y.new))
      reispe.opti =
        sum(integral((y.new.hat.opti - y.new)^2, domain = domain.y, type = 201))/
        sum(integral(sweep(y.new, 2, mu.y)^2, domain = domain.y, type = 201))
  }
  return(list(Beta = Beta.hat.opti, ncomp = c(p.opti.x, p.opti.y), y.pred = y.new.hat.opti, reispe = reispe.opti))
}

SigComp = function(x.old, y.old, domain.x, domain.y, x.new = NULL, y.new = NULL, p.max){
  K.x = min(4 + ncol(x.old) - 2, nrow(x.old) - 1)
  K.y = min(4 + ncol(y.old) - 2, nrow(y.old) - 1)
  mu.y = colMeans(y.old)
  sigcomp.obj = cv.sigcom(X = list(x.old), Y = y.old, t.x = list(domain.x), t.y = domain.y, s.n.basis = K.x, t.n.basis = K.y, K.cv = 5, upper.comp = p.max)

  Beta.hat = getcoef.sigcom(fit.obj = sigcomp.obj)$beta[[1]]

  y.new.hat = NULL
  reispe = NULL
  if (!is.null(x.new)){
    y.new.hat = pred.sigcom(fit.obj = sigcomp.obj, X.test = list(x.new))
    if (!is.null(y.new))
      reispe =
        sum(integral((y.new.hat - y.new)^2, domain = domain.y, type = 201))/
        sum(integral(sweep(y.new, 2, mu.y)^2, domain = domain.y, type = 201))
  }
  return(list(Beta = Beta.hat, reispe = reispe, ncomp = sigcomp.obj$fitted_model$opt_K))
}

NIPALS = function(x.old, y.old, domain.x, domain.y, x.new = NULL, y.new = NULL, p.max){
  nfold = 5
  n = nrow(x.old)
  cv.mat = matrix(NA, nrow = p.max, ncol = nfold)
  groups <- split(sample(1:n), rep(1:nfold, length = n))

  for (k in 1:nfold){
    x = x.old
    y = y.old[-groups[[k]], ]
    d_time_x = domain.x
    d_time_y = domain.y
    norder = 4
    K_x = min(norder + ncol(x)-2, nrow(x)-1)
    K_y = min(norder + ncol(y)-2, nrow(y)-1)

    ###Basis functions
    Phi_y = Bsplines_FDA(d_time = d_time_y, nbf = K_y, norder = norder)
    Phi_x = Bsplines_FDA(d_time = d_time_x, nbf = K_x, norder = norder)

    ###Coef w.r.t. the given basis
    m_basis_y = create.bspline.basis(rangeval = range(d_time_y), nbasis = K_y, norder = norder)
    m_basis_x = create.bspline.basis(rangeval = range(d_time_x), nbasis = K_x, norder = norder)
    S_X = Data2fd(argvals = d_time_x, y = t(x), basisobj = m_basis_x, dfscale = 1.2)$coefs
    S_Y = Data2fd(argvals = d_time_y, y = t(y), basisobj = m_basis_y, dfscale = 1.2)$coefs

    ###Observed response functions
    y_fun = t(S_Y) %*% t(Phi_y)

    ###Square roots of the inner product matrices
    Inn_prod_y_sqrt = sqrtm(inprod(m_basis_y, m_basis_y))
    Inn_prod_x_sqrt = sqrtm(inprod(m_basis_x, m_basis_x))

    ###Arguments for the PLS regression
    #Response
    Reg_y_train = t(S_Y) %*% Inn_prod_y_sqrt
    #Predictors
    Reg_mat = t(S_X) %*% Inn_prod_x_sqrt
    Reg_mat_train = Reg_mat[-groups[[k]],]
    Reg_mat_test = Reg_mat[groups[[k]],]

    ###Mean of variables
    mean_y = apply(Reg_y_train, 2, mean)
    mean_x = apply(Reg_mat, 2, mean)

    ###Center the response and predictors for training stations
    #Response
    CY = matrix(NA, nrow = nrow(Reg_y_train), ncol = ncol(Reg_y_train))
    for(i in 1:nrow(Reg_y_train))
      CY[i,] = Reg_y_train[i,] - mean_y
    #Predictors
    CX = matrix(NA, nrow = nrow(Reg_mat_train), ncol = ncol(Reg_mat_train))
    for(i in 1:nrow(Reg_mat_train))
      CX[i,] = Reg_mat_train[i,] - mean_x

    ###Center the predictors for test stations
    CX_test = matrix(NA, nrow = nrow(Reg_mat_test), ncol = ncol(Reg_mat_test))
    for(i in 1:nrow(Reg_mat_test))
      CX_test[i,] = Reg_mat_test[i,] - mean_x

    for (p in 1:p.max){
      ###NIPALS
      model_nipals = plsreg2(predictors = CX, responses = CY, comps = p, crosval = F)
      Bhat_nipals = model_nipals$reg.coefs[1:K_x,]

      ###Obtain fitted coefficients
      fitted_coefs_nipals = CX_test %*% Bhat_nipals
      for(i in 1:nrow(fitted_coefs_nipals))
        fitted_coefs_nipals[i,] = fitted_coefs_nipals[i,] + mean_y

      ###Obtain fitted functions
      fitted_functions_nipals = (fitted_coefs_nipals %*% solve(Inn_prod_y_sqrt)) %*% t(Phi_y)
      Beta_hat_nipals = Phi_x %*% solve(Inn_prod_x_sqrt) %*% Bhat_nipals %*% solve(Inn_prod_y_sqrt) %*% t(Phi_y)

      ###reispe for NIPALS
      S_Y_test = Data2fd(argvals = d_time_y, y = t(y.old[groups[[k]], ]), basisobj = m_basis_y, dfscale = 1.2)$coefs
      y_fun_test = t(S_Y_test) %*% t(Phi_y)
      reispe = sum(diag(tcrossprod(t(S_Y_test) %*% Inn_prod_y_sqrt - fitted_coefs_nipals)))/
        sum(diag(tcrossprod(sweep(t(S_Y_test) %*% Inn_prod_y_sqrt, 2, mean_y))))

      cv.mat[p, k] = reispe
    }
  }

  p.opt = which.min(apply(cv.mat, 1, mean))
  x = rbind(x.old, x.new)
  y = y.old
  d_time_x = domain.x
  d_time_y = domain.y
  norder = 4
  K_x = min(norder + ncol(x)-2, nrow(x)-1)
  K_y = min(norder + ncol(y)-2, nrow(y)-1)

  ###Basis functions
  Phi_y = Bsplines_FDA(d_time = d_time_y, nbf = K_y, norder = norder)
  Phi_x = Bsplines_FDA(d_time = d_time_x, nbf = K_x, norder = norder)

  ###Coef w.r.t. the given basis
  m_basis_y = create.bspline.basis(rangeval = range(d_time_y), nbasis = K_y, norder = norder)
  m_basis_x = create.bspline.basis(rangeval = range(d_time_x), nbasis = K_x, norder = norder)
  S_X = Data2fd(argvals = d_time_x, y = t(x), basisobj = m_basis_x, dfscale = 1.2)$coefs
  S_Y = Data2fd(argvals = d_time_y, y = t(y), basisobj = m_basis_y, dfscale = 1.2)$coefs

  ###Observed response functions
  y_fun = t(S_Y) %*% t(Phi_y)

  ###Square roots of the inner product matrices
  Inn_prod_y_sqrt = sqrtm(inprod(m_basis_y, m_basis_y))
  Inn_prod_x_sqrt = sqrtm(inprod(m_basis_x, m_basis_x))

  ###Arguments for the PLS regression
  #Response
  Reg_y_train = t(S_Y) %*% Inn_prod_y_sqrt
  #Predictors
  Reg_mat = t(S_X) %*% Inn_prod_x_sqrt
  Reg_mat_train = Reg_mat[1:nrow(x.old),]
  Reg_mat_test = Reg_mat[-(1:nrow(x.old)),]

  ###Mean of variables
  mean_y = apply(Reg_y_train, 2, mean)
  mean_x = apply(Reg_mat, 2, mean)

  ###Center the response and predictors for training stations
  #Response
  CY = matrix(NA, nrow = nrow(Reg_y_train), ncol = ncol(Reg_y_train))
  for(i in 1:nrow(Reg_y_train))
    CY[i,] = Reg_y_train[i,] - mean_y
  #Predictors
  CX = matrix(NA, nrow = nrow(Reg_mat_train), ncol = ncol(Reg_mat_train))
  for(i in 1:nrow(Reg_mat_train))
    CX[i,] = Reg_mat_train[i,] - mean_x

  ###Center the predictors for test stations
  CX_test = matrix(NA, nrow = nrow(Reg_mat_test), ncol = ncol(Reg_mat_test))
  for(i in 1:nrow(Reg_mat_test))
    CX_test[i,] = Reg_mat_test[i,] - mean_x

  ###NIPALS
  model_nipals = plsreg2(predictors = CX, responses = CY, comps = p.opt, crosval = F)
  Bhat_nipals = model_nipals$reg.coefs[1:K_x,]

  ###Obtain fitted coefficients
  fitted_coefs_nipals = CX_test %*% Bhat_nipals
  for(i in 1:nrow(fitted_coefs_nipals))
    fitted_coefs_nipals[i,] = fitted_coefs_nipals[i,] + mean_y

  ###Obtain fitted functions
  fitted_functions_nipals = (fitted_coefs_nipals %*% solve(Inn_prod_y_sqrt)) %*% t(Phi_y)
  Beta_hat_nipals = Phi_x %*% solve(Inn_prod_x_sqrt) %*% Bhat_nipals %*% solve(Inn_prod_y_sqrt) %*% t(Phi_y)

  ###reispe for NIPALS
  S_Y_test = Data2fd(argvals = d_time_y, y = t(y.new), basisobj = m_basis_y, dfscale = 1.2)$coefs
  y_fun_test = t(S_Y_test) %*% t(Phi_y)
  reispe = sum(diag(tcrossprod(t(S_Y_test) %*% Inn_prod_y_sqrt - fitted_coefs_nipals)))/
    sum(diag(tcrossprod(sweep(t(S_Y_test) %*% Inn_prod_y_sqrt, 2, mean_y))))

  return(list(Beta = Beta_hat_nipals, reispe = reispe, ncomp = p.opt))
}

SIMPLS = function(x.old, y.old, domain.x, domain.y, x.new = NULL, y.new = NULL, p.max){
  nfold = 5
  n = nrow(x.old)
  cv.mat = matrix(NA, nrow = p.max, ncol = nfold)
  groups <- split(sample(1:n), rep(1:nfold, length = n))

  for (k in 1:nfold){
    x = x.old
    y = y.old[-groups[[k]], ]
    d_time_x = domain.x
    d_time_y = domain.y
    norder = 4
    K_x = min(norder + ncol(x)-2, nrow(x)-1)
    K_y = min(norder + ncol(y)-2, nrow(y)-1)

    ###Basis functions
    Phi_y = Bsplines_FDA(d_time = d_time_y, nbf = K_y, norder = norder)
    Phi_x = Bsplines_FDA(d_time = d_time_x, nbf = K_x, norder = norder)

    ###Coef w.r.t. the given basis
    m_basis_y = create.bspline.basis(rangeval = range(d_time_y), nbasis = K_y, norder = norder)
    m_basis_x = create.bspline.basis(rangeval = range(d_time_x), nbasis = K_x, norder = norder)
    S_X = Data2fd(argvals = d_time_x, y = t(x), basisobj = m_basis_x, dfscale = 1.2)$coefs
    S_Y = Data2fd(argvals = d_time_y, y = t(y), basisobj = m_basis_y, dfscale = 1.2)$coefs

    ###Observed response functions
    y_fun = t(S_Y) %*% t(Phi_y)

    ###Square roots of the inner product matrices
    Inn_prod_y_sqrt = sqrtm(inprod(m_basis_y, m_basis_y))
    Inn_prod_x_sqrt = sqrtm(inprod(m_basis_x, m_basis_x))

    ###Arguments for the PLS regression
    #Response
    Reg_y_train = t(S_Y) %*% Inn_prod_y_sqrt
    #Predictors
    Reg_mat = t(S_X) %*% Inn_prod_x_sqrt
    Reg_mat_train = Reg_mat[-groups[[k]],]
    Reg_mat_test = Reg_mat[groups[[k]],]

    ###Mean of variables
    mean_y = apply(Reg_y_train, 2, mean)
    mean_x = apply(Reg_mat, 2, mean)

    ###Center the response and predictors for training stations
    #Response
    CY = matrix(NA, nrow = nrow(Reg_y_train), ncol = ncol(Reg_y_train))
    for(i in 1:nrow(Reg_y_train))
      CY[i,] = Reg_y_train[i,] - mean_y
    #Predictors
    CX = matrix(NA, nrow = nrow(Reg_mat_train), ncol = ncol(Reg_mat_train))
    for(i in 1:nrow(Reg_mat_train))
      CX[i,] = Reg_mat_train[i,] - mean_x

    ###Center the predictors for test stations
    CX_test = matrix(NA, nrow = nrow(Reg_mat_test), ncol = ncol(Reg_mat_test))
    for(i in 1:nrow(Reg_mat_test))
      CX_test[i,] = Reg_mat_test[i,] - mean_x

    ###SIMPLS
    model_simpls = pls.regression(Xtrain = CX, Ytrain = CY, Xtest = CX, ncomp = 1:p.max, unit.weights=FALSE)
    for (p in 1:p.max){
      Bhat_simpls = model_simpls$B[ , , p]

      ###Obtain fitted coefficients
      fitted_coefs_simpls = CX_test %*% Bhat_simpls
      for(i in 1:nrow(fitted_coefs_simpls))
        fitted_coefs_simpls[i,] = fitted_coefs_simpls[i,] + mean_y

      ###Obtain fitted functions
      fitted_functions_simpls = fitted_coefs_simpls %*% solve(Inn_prod_y_sqrt) %*% t(Phi_y)
      Beta_hat_simpls = Phi_x %*% solve(Inn_prod_x_sqrt) %*% Bhat_simpls %*% solve(Inn_prod_y_sqrt) %*% t(Phi_y)

      ###reispe for SIMPLS
      S_Y_test = Data2fd(argvals = d_time_y, y = t(y.old[groups[[k]], ]), basisobj = m_basis_y, dfscale = 1.2)$coefs
      y_fun_test = t(S_Y_test) %*% t(Phi_y)
      reispe = sum(diag(tcrossprod(t(S_Y_test) %*% Inn_prod_y_sqrt - fitted_coefs_simpls)))/
        sum(diag(tcrossprod(sweep(t(S_Y_test) %*% Inn_prod_y_sqrt, 2, mean_y))))

      cv.mat[p, k] = reispe
    }
  }

  p.opt = which.min(apply(cv.mat, 1, mean))
  x = rbind(x.old, x.new)
  y = y.old
  d_time_x = domain.x
  d_time_y = domain.y
  norder = 4
  K_x = min(norder + ncol(x)-2, nrow(x)-1)
  K_y = min(norder + ncol(y)-2, nrow(y)-1)

  ###Basis functions
  Phi_y = Bsplines_FDA(d_time = d_time_y, nbf = K_y, norder = norder)
  Phi_x = Bsplines_FDA(d_time = d_time_x, nbf = K_x, norder = norder)

  ###Coef w.r.t. the given basis
  m_basis_y = create.bspline.basis(rangeval = range(d_time_y), nbasis = K_y, norder = norder)
  m_basis_x = create.bspline.basis(rangeval = range(d_time_x), nbasis = K_x, norder = norder)
  S_X = Data2fd(argvals = d_time_x, y = t(x), basisobj = m_basis_x, dfscale = 1.2)$coefs
  S_Y = Data2fd(argvals = d_time_y, y = t(y), basisobj = m_basis_y, dfscale = 1.2)$coefs

  ###Observed response functions
  y_fun = t(S_Y) %*% t(Phi_y)

  ###Square roots of the inner product matrices
  Inn_prod_y_sqrt = sqrtm(inprod(m_basis_y, m_basis_y))
  Inn_prod_x_sqrt = sqrtm(inprod(m_basis_x, m_basis_x))

  ###Arguments for the PLS regression
  #Response
  Reg_y_train = t(S_Y) %*% Inn_prod_y_sqrt
  #Predictors
  Reg_mat = t(S_X) %*% Inn_prod_x_sqrt
  Reg_mat_train = Reg_mat[1:nrow(x.old),]
  Reg_mat_test = Reg_mat[-(1:nrow(x.old)),]

  ###Mean of variables
  mean_y = apply(Reg_y_train, 2, mean)
  mean_x = apply(Reg_mat, 2, mean)

  ###Center the response and predictors for training stations
  #Response
  CY = matrix(NA, nrow = nrow(Reg_y_train), ncol = ncol(Reg_y_train))
  for(i in 1:nrow(Reg_y_train))
    CY[i,] = Reg_y_train[i,] - mean_y
  #Predictors
  CX = matrix(NA, nrow = nrow(Reg_mat_train), ncol = ncol(Reg_mat_train))
  for(i in 1:nrow(Reg_mat_train))
    CX[i,] = Reg_mat_train[i,] - mean_x

  ###Center the predictors for test stations
  CX_test = matrix(NA, nrow = nrow(Reg_mat_test), ncol = ncol(Reg_mat_test))
  for(i in 1:nrow(Reg_mat_test))
    CX_test[i,] = Reg_mat_test[i,] - mean_x

  ###SIMPLS
  model_simpls = pls.regression(Xtrain = CX, Ytrain = CY, Xtest = CX, ncomp = p.opt, unit.weights=FALSE)
  Bhat_simpls = model_simpls$B

  ###Obtain fitted coefficients
  fitted_coefs_simpls = CX_test %*% Bhat_simpls
  for(i in 1:nrow(fitted_coefs_simpls))
    fitted_coefs_simpls[i,] = fitted_coefs_simpls[i,] + mean_y

  ###Obtain fitted functions
  fitted_functions_simpls = fitted_coefs_simpls %*% solve(Inn_prod_y_sqrt) %*% t(Phi_y)
  Beta_hat_simpls = Phi_x %*% solve(Inn_prod_x_sqrt) %*% Bhat_simpls %*% solve(Inn_prod_y_sqrt) %*% t(Phi_y)

  ###reispe for SIMPLS
  S_Y_test = Data2fd(argvals = d_time_y, y = t(y.new), basisobj = m_basis_y, dfscale = 1.2)$coefs
  y_fun_test = t(S_Y_test) %*% t(Phi_y)
  reispe = sum(diag(tcrossprod(t(S_Y_test) %*% Inn_prod_y_sqrt - fitted_coefs_simpls)))/
    sum(diag(tcrossprod(sweep(t(S_Y_test) %*% Inn_prod_y_sqrt, 2, mean_y))))

  return(list(Beta = Beta_hat_simpls, reispe = reispe, ncomp = p.opt))
}

orthoBasis = function(order, denseGrid, type, normalized = T){
  if (!("orthopolynom" %in% rownames(installed.packages())))
    install.packages("orthopolynom")
  library(orthopolynom)

  if (type == 'shiftedLegendre') {
    poly.value = polynomial.values(slegendre.polynomials(n = max(order), normalized = normalized), x = denseGrid)
    res = matrix(NA, nrow = length(order), ncol = length(denseGrid))
    for (i in 1:length(order)){
      res[i, ] = poly.value[[order[i]+1]]
    }
  }
  return(res)
}

#' expand.basis is to expand 1-argument functional data with respect to a basis.
#' The output martrix C is guaranteed to have full column rank.
#' @param x is a matrix. Each row represents a sample.
#' @param nOrder=4 by default indicates the cubic B-spline.
### Dongjiao add K to present the number of spline basis.
## the original function is expand.basis = function(x, domain.x, nOrder = 4, basis.name)
expand.basis = function(x, domain.x, nOrder = 4, basis.name, K){
  if (!("fda" %in% rownames(installed.packages()))) install.packages("fda")
  if (!all(c(
    exists("create.bspline.basis", mode="function"),
    exists("int2Lfd", mode="function"),
    exists("fdPar", mode="function"),
    exists("Data2fd", mode="function"),
    exists("inprod", mode="function"),
    exists('deriv.fd', mode='function')
  ))) library("fda")

  timePts = domain.x
  #K = min(nOrder+ncol(x)-2, nrow(x)-1) # number of basis functions#dongjiao comment this
  #K=10#dongjiao add this
  #K=20#dongjiao add this
  if (basis.name == "bspline")
    basis = create.bspline.basis(rangeval=range(timePts), nbasis=K, norder = nOrder, names="bspl")
  if (basis.name == "fourier")
    basis = create.fourier.basis(rangeval = c(0,1), nbasis = K)
  D2Lfd = int2Lfd(m=2)

  # tune the value of lambda for smoothing x
  log10lambda = (-20):20
  gcvsave = NULL
  for (i in log10lambda) {
    lambda = 10^i
    D2fdPar = fdPar(basis, D2Lfd, lambda)
    smooth = smooth.basis(y=t(x), argvals=timePts, D2fdPar, dfscale=1.2)
    gcvsave = c(gcvsave, sum(smooth$gcv))
  }

  lambda.opt = 10^log10lambda[which.min(gcvsave)]
  D2fdPar = fdPar(basis, D2Lfd, lambda = lambda.opt)
  fd.smooth = smooth.basis(y=t(x), argvals=timePts, D2fdPar, dfscale=1.2)

  return(list(fdObj = fd.smooth$fd, coef = t(fd.smooth$fd$coefs), basis = basis, D2fdPar = D2fdPar, lambda = lambda))
}

#' compute p.max
#'
#'
## Dongjiao add K to presnt the number of basis
## The original function is pUpper.compu = function(x, domain.x, proportion, basis.name)
pUpper.compu = function(x, domain.x, proportion, basis.name, K){
  mu.x = colMeans(x)
  #expansion.x = expand.basis(x, domain.x, basis.name = basis.name)#Dongjiao comment this
  expansion.x = expand.basis(x, domain.x, basis.name = basis.name, K=K)#Dongjiao added this
  W = inprod(expansion.x$basis, expansion.x$basis)
  # x.fd = smooth.basisPar(argvals=1:ncol(x),
  #                        y = t(sweep(x, 2, mu.x)),
  #                        fdobj=expansion.x$basis,
  #                        Lfdobj=int2Lfd(2),
  #                        lambda=expansion.x$lambda)
  #x.pca.obj = pca.fd(expansion.x$fdObj, nharm = min(dim(x) - 1), centerfns = T)#dongjiao comment this
  x.pca.obj = pca.fd(expansion.x$fdObj, nharm = min(dim(x) - 1, dim(expansion.x$fdObj$coefs)-1), centerfns = T)
  pMax = which.max(cumsum(x.pca.obj$values)/sum(x.pca.obj$values) >= proportion)
  return(pMax)
}

#' integrate error measures

creatErrMat = function(errLst){
  nrowErrMat = min(lengths(errLst))
  errMat = NULL
  if (nrowErrMat > 0){
    for (i in 1:length(errLst)){
      errMat = cbind(errMat, errLst[[i]][1:nrowErrMat])
    }
    colnames(errMat) = names(errLst)
  }
  return(errMat)
}

#'  boxplot of error measures

boxplotErr = function(nrowplot, errMat, type){
  if (!("ggplot2" %in% rownames(installed.packages())))
    install.packages("ggplot2")
  library(ggplot2)
  library(reshape2)

  melton = melt(
    data.frame(errMat[1:nrowplot,],
              Replication = 1:nrowplot),
              id.vars = "Replication")
  bplot = ggplot(melton, aes(x = variable, y = value, colour = variable)) +
    # geom_violin() +
    geom_boxplot(outlier.shape = NA, width = .5) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = '', y = type, title = '') +
    theme_bw() +
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 25),
          axis.text.x = element_text(size = 25, angle = 60, hjust = .8)
          )

  plot(bplot)
  file = ifelse(simu,
                paste0(type, '_',
                       'Simu', simuCase, '_',
                       nrowplot, 'repeats_',
                       propTrain * 100, 'train_',
                       rho * 1e1L, 'rho_',
                       round(SNR), 'SNR_',
                       FVE.x * 1e2L, 'FVEx.pdf'),
                paste0(type, '_',
                       'Real', realCase, '_',
                       nrowplot, 'repeats_',
                       propTrain*100, 'train_',
                       FVE.x * 1e2L, 'FVEx.pdf')
  )
  ggsave(file = file,
         width = 6,
         height = 8,
         units = 'in',
         dpi = 300,
         path = "figure")
}

PlotMultiCurve = function(x, y){
  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  if (!("tidyr" %in% rownames(installed.packages()))) install.packages("tidyr")
  library(ggplot2)

  curves.simu = data.frame(ID=1L:nrow(x), day=x, class=factor(y))
  curves.simu = tidyr::gather(curves.simu, day, value, -c(ID, class))
  curvesplot.simu = ggplot(curves.simu, aes(x=day, y=value, colour=class)) +
    geom_line(aes(group=ID))
  curvesplot.simu +
    ggtitle("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5))
  }

Bsplines_FDA = function(d_time, nbf, norder=4){
  require(fda)
  basis = create.bspline.basis(rangeval = range(d_time), nbasis = nbf, norder)
  Phi = eval.basis(evalarg = d_time, basisobj = basis)
  return(Phi)
}

#'
#'
#'

impute.curve = function(x){
  if (!("spatialEco" %in% rownames(installed.packages()))) install.packages("spatialEco")
  library(spatialEco)

  gcv = NULL
  lambda.vec = 10^(0:6)

  for (lambda in lambda.vec){
    x.smooth = array(NA, dim = dim(x))
    for (i in 1:nrow(x)){
      x.smooth[i, ] = poly.regression(x[i, ], s = lambda)
    }
    gcv = c(gcv, sum(((x-x.smooth)[!is.na(x-x.smooth)])^2))
  }
  lambda = lambda.vec[which.min(gcv)]
  for (i in 1:nrow(x)){
    x.smooth[i, ] = poly.regression(x[i, ], s = lambda, impute = T, na.only = T)
  }

  return(x.smooth)
}

PlotMultiCurve.x = function(x, domain.x, group){

  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  if (!("tidyr" %in% rownames(installed.packages()))) install.packages("tidyr")
  library(ggplot2)

  data.reform = data.frame(ID = 1:nrow(x), class = factor(group), time = x)
  data.reform = tidyr::gather(data.reform, time, value, -c(ID, class))
  data.reform$time = rep(domain.x, each = nrow(x))
  curvesplot.simu =
    ggplot(data.reform, aes(x=time, y=value, colour=class, linetype=class)) +
    geom_line(aes(group=ID), color = 'green')
  curvesplot.simu +
    coord_cartesian(ylim = range(x)+c(-.1, .1)) +
    labs(x='s', y='X(s)') +
    ggtitle("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.title = element_text(hjust = .5),
          text = element_text(size = 25)
    )
}

PlotMultiCurve.y = function(x, domain.y, group){

  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  if (!("tidyr" %in% rownames(installed.packages()))) install.packages("tidyr")
  library(ggplot2)

  data.reform = data.frame(ID = 1:nrow(x), class = factor(group), time = x)
  data.reform = tidyr::gather(data.reform, time, value, -c(ID, class))
  data.reform$time = rep(domain.y, each = nrow(x))
  curvesplot.simu =
    ggplot(data.reform, aes(x=time, y=value, colour=class, linetype=class)) +
    geom_line(aes(group=ID), color='blue')
  curvesplot.simu +
    coord_cartesian(ylim = range(x)+c(-.1, .1)) +
    labs(x='t', y='Y(t)') +
    ggtitle("") +
    theme_bw()+
    theme(legend.position = "none",
          panel.border = element_blank(),
          plot.title = element_text(hjust = .5),
          text = element_text(size = 25)
    )
}

