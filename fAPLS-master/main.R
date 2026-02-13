#' Put main.R & supportFun.R into the identical folder.
#' In the same folder, create two subfolders named figure' & 'Rimage', respectively.

rm(list = ls())
if (!("rstudioapi" %in% rownames(installed.packages()))) 
  install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("functions.R")

#### Global settings
options(warn = -1, digits = 4)
set.seed(1)
FVE.x = .99
FVE.y = .99
RR = 50 # number of replication
propTrain = .8 # proportion of real dataset left for training
realCase = 1 # 1 for DTI; 2 for BG
simu = F # T for simulation, F for real data
simuCase = 1 # 1 for beta = P_2(s)P_2(t); 2 for beta = P_4(s)P_4(t)
rho = .9 # control the autocorrelation for epsilon
sigma = .05^.5 # sigma for epsilon

ridge.par = c(0)
tune.method = 'CV' # 'CV' or 'GCV'

#### Simulated data ####

if (simu == T) {
  N = 300
  
  if (simuCase %in% c(1)){
    tRange = c(0, 1)
    domain.x = seq(from = tRange[1], to = tRange[2], length.out = 101)
    domain.y = seq(from = tRange[1], to = tRange[2], length.out = 100)
    # Beta.true
    Lambda = c(100, 10, 1)
    J = length(Lambda) # number of eigenfunctions
    eigenFun.x = orthoBasis(order = 2:(J+1), denseGrid = domain.x, type = 'shiftedLegendre', normalized = T)
    eigenFun.y = orthoBasis(order = 2:(J+1), denseGrid = domain.y, type = 'shiftedLegendre', normalized = T)
    
    Beta.true = matrix(0, nrow = length(domain.x), ncol = length(domain.y))
    for (j in 1:J){
      Beta.true = Beta.true + eigenFun.x[j, ] %o% eigenFun.y[j, ]
    }

    Sigma.epsilon = matrix(NA, nrow = length(domain.y), ncol = length(domain.y))
    for (i in 1:nrow(Sigma.epsilon)){
      for (j in 1:ncol(Sigma.epsilon)){
        Sigma.epsilon[i, j] = sigma^2 * (rho^(abs(domain.y[i] - domain.y[j])))
      }
    }

    # Create x and y
    x.all = list()
    y.all = list()
    for (R in 1:RR) {
      eigenScore = MASS::mvrnorm(n = N, mu = numeric(J), Sigma = diag(Lambda))
      x.all[[R]] = eigenScore %*% eigenFun.x
      epsilon = MASS::mvrnorm(n = N, mu = numeric(length(domain.y)), Sigma = Sigma.epsilon)
      y.all[[R]] = integral(x.all[[R]], Beta.true, domain = domain.x, type = 222) + epsilon
    }
    
    SNR = (
      mean(integral(y.all[[RR]]^2, domain = domain.y, type = 201))/mean(integral(epsilon^2, domain = domain.y, type = 201)) - 1
    )^.5;SNR
  }
  
  if (simuCase %in% c(2)){
    tRange = c(0, 1)
    domain.x = seq(from = tRange[1], to = tRange[2], length.out = 101)
    domain.y = seq(from = tRange[1], to = tRange[2], length.out = 100)
    Sigma.x1 = matrix(NA, nrow = length(domain.x), ncol = length(domain.x))
    for (i in 1:nrow(Sigma.x1)){
      for (j in 1:ncol(Sigma.x1)){
        Sigma.x1[i, j] = exp(-(10 * (domain.x[i] - domain.x[j]))^2)
      }
    }
    Sigma.x2 = matrix(NA, nrow = length(domain.x), ncol = length(domain.x))
    for (i in 1:nrow(Sigma.x2)){
      for (j in 1:ncol(Sigma.x2)){
        Sigma.x2[i, j] = (1 + 20*abs(domain.x[i]-domain.x[j]) + 1/3*(20*(domain.x[i]-domain.x[j]))^2) *
          exp(-20*abs(domain.x[i]-domain.x[j]))
      }
    }
    Sigma.y2 = matrix(NA, nrow = length(domain.y), ncol = length(domain.y))
    for (i in 1:nrow(Sigma.y2)){
      for (j in 1:ncol(Sigma.y2)){
        Sigma.y2[i, j] = (1 + 20*abs(domain.y[i]-domain.y[j]) + 1/3*(20*(domain.y[i]-domain.y[j]))^2) *
          exp(-20*abs(domain.y[i]-domain.y[j]))
      }
    }
    Sigma.epsilon = matrix(NA, nrow = length(domain.y), ncol = length(domain.y))
    for (i in 1:nrow(Sigma.epsilon)){
      for (j in 1:ncol(Sigma.epsilon)){
        Sigma.epsilon[i, j] = sigma^2 * (rho^(abs(domain.y[i] - domain.y[j])))
      }
    }
    
    Beta.0.true = MASS::mvrnorm(n = 1, mu = numeric(length(domain.y)), Sigma = Sigma.y2)
    Beta.true = 5 * (
      MASS::mvrnorm(n = 1, mu = numeric(length(domain.x)), Sigma = Sigma.x2) %o%
      MASS::mvrnorm(n = 1, mu = numeric(length(domain.y)), Sigma = Sigma.y2)
    ) + 10 * (
      MASS::mvrnorm(n = 1, mu = numeric(length(domain.x)), Sigma = Sigma.x2) %o%
        MASS::mvrnorm(n = 1, mu = numeric(length(domain.y)), Sigma = Sigma.y2)
    ) + 15 * (
      MASS::mvrnorm(n = 1, mu = numeric(length(domain.x)), Sigma = Sigma.x2) %o%
        MASS::mvrnorm(n = 1, mu = numeric(length(domain.y)), Sigma = Sigma.y2)
    )
    
    # Create x and y
    x.all = list()
    y.all = list()
    for (R in 1:RR){
      x.all[[R]] = MASS::mvrnorm(n = N, mu = numeric(length(domain.x)), Sigma = Sigma.x1)
      epsilon = MASS::mvrnorm(n = N, mu = numeric(length(domain.y)), Sigma = Sigma.epsilon)
      y.all[[R]] = matrix(Beta.0.true, nrow = N, ncol = length(domain.y), byrow = T) + 
        integral(x.all[[R]], Beta.true, domain = domain.x, type = 222) + epsilon
    }
    
    SNR = (
      mean(integral(y.all[[RR]]^2, domain = domain.y, type = 201))/mean(integral(epsilon^2, domain = domain.y, type = 201)) - 1
    )^.5;SNR
  }
  
  if (simuCase %in% c(3)){
    tRange = c(0, 1)
    domain.x = seq(from = tRange[1], to = tRange[2], length.out = 101)
    domain.y = seq(from = tRange[1], to = tRange[2], length.out = 100)
    
    Beta.true = sin(pi * domain.x) %o% cos(2 * pi * domain.y)
    x.fun = function(domain){
      m = 10
      x = numeric(length(domain))
      zeta = rnorm(2*m)
      for (i in 1:length(domain)){
        x[i] = sum((1/rep(1:m, times = 2)^2) * zeta * c(sin((1:m) * pi * domain[i]), cos((1:m) * pi * domain[i])))
      }
      return(x)
    }
    
    Beta.0.true = 2 * exp(-(domain.y-1)^2)
    Sigma.epsilon = matrix(NA, nrow = length(domain.y), ncol = length(domain.y))
    for (i in 1:nrow(Sigma.epsilon)){
      for (j in 1:ncol(Sigma.epsilon)){
        Sigma.epsilon[i, j] = sigma^2 * (rho^(abs(domain.y[i] - domain.y[j])))
      }
    }
    
    # Create x and y
    x.all = list()
    y.all = list()
    for (R in 1:RR){
      x.all[[R]] = matrix(NA, nrow = N, ncol = length(domain.x))
      for (i in 1:N){
        x.all[[R]][i, ] = x.fun(domain.x)
      }
      epsilon = MASS::mvrnorm(n = N, mu = numeric(length(domain.y)), Sigma = Sigma.epsilon)
      y.all[[R]] = matrix(Beta.0.true, nrow = N, ncol = length(domain.y), byrow = T) + 
        integral(x.all[[R]], Beta.true, domain = domain.x, type = 222) + epsilon
    }
    
    SNR = (
      mean(integral(y.all[[RR]]^2, domain = domain.y, type = 201))/mean(integral(epsilon^2, domain = domain.y, type = 201)) - 1
    )^.5;SNR
  }
  
}


#### Real data ####
source("functions.R")
if (simu == F){
  if (realCase == 1) {
    # cca vs. rcst (DTI in R package refund)
    
    if (!("refund" %in% rownames(installed.packages()))) install.packages("refund")
    Beta.true = NULL
    x = impute.curve(refund::DTI$cca)
    y = impute.curve(refund::DTI$rcst)
    domain.x = (1:ncol(x))/ncol(x)
    domain.y = (1:ncol(y))/ncol(y)
    N = nrow(x)
    PlotMultiCurve.x(x[1:100,], domain.x, numeric(100))
    PlotMultiCurve.y(y[1:100,], domain.y, numeric(100))
    
    p.max = pUpper.compu(x, domain.x, FVE.x, basis.name = "bspline")
  }
  
  if (realCase == 2){
    # Boy's gait (gait in R package fda)
    Beta.true = NULL
    x = t(fda::gait[,,"Hip Angle"])
    y = t(fda::gait[,,"Knee Angle"])
    domain.x = seq(from=0.025, to=0.975, by=0.05)
    domain.y = seq(from=0.025, to=0.975, by=0.05)
    N = nrow(x)
    PlotMultiCurve.x(x, domain.x, numeric(N))
    PlotMultiCurve.y(y, domain.y, numeric(N))
    
    p.max = pUpper.compu(x, domain.x, FVE.x, basis.name = "fourier")
  }
}

# Dangerous! clear existing result
res.fAPLS = list(); time.fAPLS = numeric(); reispe.fAPLS = numeric(); reisee.fAPLS = numeric(); ncomp.fAPLS = numeric()
res.SigComp = list(); time.SigComp = numeric(); reispe.SigComp = numeric(); reisee.SigComp = numeric(); ncomp.SigComp = numeric()
res.NIPALS = list(); time.NIPALS = numeric(); reispe.NIPALS = numeric(); reisee.NIPALS = numeric(); ncomp.NIPALS = numeric()
res.SIMPLS = list(); time.SIMPLS = numeric(); reispe.SIMPLS = numeric(); reisee.SIMPLS = numeric(); ncomp.SIMPLS = numeric()

# Check current progress
check.fAPLS = length(time.fAPLS)
check.SigComp = length(time.SigComp)
check.NIPALS = length(time.NIPALS)
check.SIMPLS = length(time.SIMPLS)

################ Replica ##################
source("functions.R")
for (R in 1:RR) {
  
  if (simu == T) {
    
    sampIdx = sample(1:N, round(N * propTrain))
    x.old = x.all[[R]][sampIdx, ]
    y.old = y.all[[R]][sampIdx, ]
    x.new = x.all[[R]][-sampIdx, ]
    y.new = y.all[[R]][-sampIdx, ]
    p.max = pUpper.compu(x.all[[R]], domain.x, FVE.x, basis.name = "bspline")

    # regression
    if (R > check.NIPALS){
      ptm0 = proc.time()[3]
      res.NIPALS[[R]] = NIPALS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.NIPALS[R] = ptm1 - ptm0
      
      reispe.NIPALS[R] = res.NIPALS[[R]]$reispe
      ncomp.NIPALS[R] = res.NIPALS[[R]]$ncomp
      
      if (!is.null(Beta.true)){
        integ1 = integral((res.NIPALS[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
        integ2 = integral((res.NIPALS[[R]]$Beta)^2, domain = domain.y, type = 201)
        reisee.NIPALS[R] = integral(integ1, domain = domain.x, type = 100)/
          integral(integ2, domain = domain.x, type = 100)
      }
    }
    
    if (R > check.SIMPLS){
      ptm0 = proc.time()[3]
      res.SIMPLS[[R]] = SIMPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SIMPLS[R] = ptm1 - ptm0
      
      reispe.SIMPLS[R] = res.SIMPLS[[R]]$reispe
      ncomp.SIMPLS[R] = res.SIMPLS[[R]]$ncomp
      
      if (!is.null(Beta.true)){
        integ1 = integral((res.SIMPLS[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
        integ2 = integral((res.SIMPLS[[R]]$Beta)^2, domain = domain.y, type = 201)
        reisee.SIMPLS[R] = integral(integ1, domain = domain.x, type = 100)/
          integral(integ2, domain = domain.x, type = 100)
      }
    }
  
    if (R > check.fAPLS){
      ptm0 = proc.time()[3]
      res.fAPLS[[R]] = fAPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, ridge.par, p.max, tune = tune.method)
      ptm1 = proc.time()[3]
      time.fAPLS[R] = ptm1 - ptm0
      
      reispe.fAPLS[R] = res.fAPLS[[R]]$reispe
      ncomp.fAPLS[R] = res.fAPLS[[R]]$ncomp
      
      if (!is.null(Beta.true)){
        integ1 = integral((res.fAPLS[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
        integ2 = integral((res.fAPLS[[R]]$Beta)^2, domain = domain.y, type = 201)
        reisee.fAPLS[R] = integral(integ1, domain = domain.x, type = 100)/
          integral(integ2, domain = domain.x, type = 100)
      }
    }
    
    if (R > check.SigComp){
      ptm0 = proc.time()[3]
      res.SigComp[[R]] = SigComp(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SigComp[R] = ptm1 - ptm0
      
      reispe.SigComp[R] = res.SigComp[[R]]$reispe
      ncomp.SigComp[R] = res.SigComp[[R]]$ncomp
      
      if (!is.null(Beta.true)){
        integ1 = integral((res.SigComp[[R]]$Beta - Beta.true)^2, domain = domain.y, type = 201)
        integ2 = integral((res.SigComp[[R]]$Beta)^2, domain = domain.y, type = 201)
        reisee.SigComp[R] = integral(integ1, domain = domain.x, type = 100)/
          integral(integ2, domain = domain.x, type = 100)
      }
    }
  }
  
  if (simu == F) {
    
    sampIdx = sample(1:N, round(N * propTrain))
    x.old = x[sampIdx, ]
    y.old = y[sampIdx, ]
    x.new = x[-sampIdx, ]
    y.new = y[-sampIdx, ]
    
    # regression
    if (R > check.NIPALS){
      ptm0 = proc.time()[3]
      res.NIPALS[[R]] = NIPALS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.NIPALS[R] = ptm1 - ptm0

      reispe.NIPALS[R] = res.NIPALS[[R]]$reispe
      ncomp.NIPALS[R] = res.NIPALS[[R]]$ncomp
    }
    
    if (R > check.SIMPLS){
      ptm0 = proc.time()[3]
      res.SIMPLS[[R]] = SIMPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SIMPLS[R] = ptm1 - ptm0
      
      reispe.SIMPLS[R] = res.SIMPLS[[R]]$reispe
      ncomp.SIMPLS[R] = res.SIMPLS[[R]]$ncomp
    }
    
    if (R > check.fAPLS){
      ptm0 = proc.time()[3]
      res.fAPLS[[R]] = fAPLS(x.old, y.old, domain.x, domain.y, x.new, y.new, ridge.par, p.max, tune = tune.method)
      ptm1 = proc.time()[3]
      time.fAPLS[R] = ptm1 - ptm0
      
      reispe.fAPLS[R] = res.fAPLS[[R]]$reispe
      ncomp.fAPLS[R] = res.fAPLS[[R]]$ncomp
    }
    
    if (R > check.SigComp){
      ptm0 = proc.time()[3]
      res.SigComp[[R]] = SigComp(x.old, y.old, domain.x, domain.y, x.new, y.new, p.max)
      ptm1 = proc.time()[3]
      time.SigComp[R] = ptm1 - ptm0
      
      reispe.SigComp[R] = res.SigComp[[R]]$reispe
      ncomp.SigComp[R] = res.SigComp[[R]]$ncomp
    }
  }
  
  if (R == 1){
    file = ifelse(simu,
                  paste0('Rimage/',
                         'Simu', simuCase, '_',
                         RR, 'repeats_',
                         propTrain * 100, 'train_',
                         rho * 10, 'rho_',
                         round(SNR), 'SNR_',
                         FVE.x * 1e2L, 'FVEx.RData'),
                  paste0('Rimage/',
                         'Real', realCase, '_', 
                         RR, 'repeats_',
                         propTrain * 100, 'train_',
                         FVE.x * 1e2L, 'FVEx.RData')
    )
  }
  if (R %% 50 == 0){
    cat(R, '\n')
  }else 
    cat(R)
  
  # save the R space
  if (R == RR){
    if (simu == T){
      # shrink the R space for simulation
      rm(x.all, y.all)
    }
    save.image(file = file)
  }
}

################ Plots and summaries of results ##################
source("functions.R")

# ReISEPE matrix and boxplots
reispeLst = list(
  reispe.fAPLS
  ,reispe.SigComp
  ,reispe.NIPALS
  ,reispe.SIMPLS
)
names(reispeLst) = c(
  'fAPLS'
  ,'SigComp'
  ,'NIPALS'
  ,'SIMPLS'
)
reispeMat = creatErrMat(reispeLst)
round(cbind(
  apply(reispeMat, 2, mean) * 100,
  apply(reispeMat, 2, sd) * 100
), digits = 2)
boxplotErr(RR, reispeMat, type = 'ReISPE')

# ReISEE matrix and boxplots
if (simu == T & !is.null(Beta.true)){
  reiseeLst = list(
    reisee.fAPLS
    ,reisee.SigComp
    ,reisee.NIPALS
    ,reisee.SIMPLS
  )
  names(reiseeLst) = names(reispeLst)
  reiseeMat = creatErrMat(reiseeLst)
  round(cbind(
    apply(reiseeMat, 2, mean) * 100,
    apply(reiseeMat, 2, sd) * 100
  ), digits = 2)
  boxplotErr(RR, reiseeMat, type = 'ReISEE')
}

# Number of components
ncompLst = list(
  ncomp.fAPLS
  ,ncomp.SigComp
  ,ncomp.NIPALS
  ,ncomp.SIMPLS
)
names(ncompLst) = names(reispeLst)
ncompMat = creatErrMat(ncompLst)
round(cbind(
  apply(ncompMat, 2, mean),
  apply(ncompMat, 2, sd)
), digits = 1)

# Running time
timeLst = list(
  time.fAPLS
  ,time.SigComp
  ,time.NIPALS
  ,time.SIMPLS
)
names(timeLst) = names(reispeLst)
timeMat = creatErrMat(timeLst)
round(cbind(
  apply(timeMat, 2, sum)
), digits = 1)
