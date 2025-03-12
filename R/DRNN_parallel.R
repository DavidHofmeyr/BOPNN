one_hot_scale <- function(X, convert = c()){
  if(length(convert)==0){
    idxs <- 1:ncol(X)
    scls <- apply(X, 2, function(z) ifelse(sd(z)==0, 0, 1/sd(z)))
    out <- as.matrix(X*matrix(scls, nrow(X), ncol(X), byrow = T))
  } 
  else{
    out <- c()
    idxs <- c()
    scls <- c()
    for(j in 1:ncol(X)){
      if(j%in%convert){
        u <- unique(X[,j])
        if(length(u) > 1){
          for(d in 1:(length(u)-1)){
            out <- cbind(out, (X[,j]==u[d])*1)
            idxs <- c(idxs, j)
            scls <- c(scls, 1)
          }
        }
        else{
          out <- cbind(out, numeric(nrow(X)))
          idxs <- c(idxs, j)
          scls <- c(scls, 1)
        }
      }
      else{
        out <- cbind(out, X[,j]*ifelse(sd(X[,j])==0, 0, 1/sd(X[,j])))
        idxs <- c(idxs, j)
        scls <- c(scls, ifelse(sd(X[,j])==0, 0, 1/sd(X[,j])))
      }
    }
  }
  list(out, idxs, scls)
}


BOPNN <- function(data_df, nn = 3, ndim = NULL, dim_prop1 = .5, ndim_samp = NULL, dim_prop0 = min(.75, 5/sqrt(ncol(data_df))), n_prop = 0.63, nbag = 100, coarseness = 6, vimp = FALSE, eps = 0, cores = 1, OOB = FALSE, case_weights = NULL, plotMat = FALSE){
  if(!('y'%in%names(data_df))) stop('data_df must contain a response named "y"')
  n <- nrow(data_df)
  
  nnp <- nn
  
  x0 <- data_df[,names(data_df)!='y']
  y0 <- data_df$y
  
  if(is.null(ndim_samp)) ndim_samp <- max(2, floor(ncol(x0)*dim_prop0))
  if(is.null(ndim)) ndim <- max(2, floor(ndim_samp*dim_prop1))
  n_samp <- floor(n*n_prop)
  
  convert <- which(unlist(lapply(x0, function(l){
    if(is.factor(l)) TRUE
    else if(typeof(l)=='integer' & length(unique(l)) <= coarseness) TRUE
    else FALSE
  })))
  
  scaled <- one_hot_scale(x0, convert)
  
  X <- scaled[[1]]
  varx <- apply(X, 2, var)
  if(min(varx)<1e-20) X[,varx < 1e-20] <- rnorm(length(X[,varx < 1e-20]))*1e-5
  idxs <- scaled[[2]]
  scales <- scaled[[3]]
  
  d <- ncol(x0)
  
  ylevs <- levels(as.factor(y0))
  
  y <- as.numeric(as.factor(y0))
  
  Y <- sapply(1:max(y), function(k) y==k)*1

  if(cores > 1){
    cluster <- makeCluster(cores)
    clusterExport(cluster, c(ls(), 'nn2'), envir = environment())
    
    mods <- parLapply(cl = cluster, 1:nbag,
                      function(x){
                        if(is.null(case_weights)) smp <- sample(1:n, n_samp)
                        else smp <- sample(1:n, n_samp, prob = case_weights, replace = TRUE)
                        dsmp <- sample(1:d, ndim_samp)
                        duse <- which(idxs%in%dsmp)
                        
                        idin <- 1:n
                        idout <- 1:n
                        for(k in 1:max(y)){
                          relix <- smp[y[smp]==k]
                          relix2 <- smp[y[smp]!=k]
                          kuse2 <- min(nnp, length(relix2))
                          ny <- length(relix)
                          nny <- length(relix2)
                          if(ny==1){
                            idin[relix] <- relix
                            idout[relix] <- relix2[nn2(X[relix2,duse], query = matrix(X[relix,duse], nrow = 1), k = kuse2, eps = eps)$nn.idx[,kuse2]]
                          }
                          else if(ny>0){
                            kuse <- min(nnp+1,ny)
                            idin[relix] <- relix[nn2(X[relix,duse], k = kuse, eps = eps)$nn.idx[,kuse]]
                            if(nny > 1) idout[relix] <- relix2[nn2(X[relix2,duse], query = X[relix,duse], k = kuse2, eps = eps)$nn.idx[,kuse2]]
                            else if(nny == 1) idout[relix] <- relix2[nn2(matrix(X[relix2,duse], nrow = 1), query = X[relix,duse], k = kuse2, eps = eps)$nn.idx[,kuse2]]
                          }
                        }
                        E <- eigen(t(X[smp,duse]-X[idin[smp],duse])%*%(X[smp,duse]-X[idin[smp],duse]) + diag(length(duse))*1e-5)
                        sqinv <- t(E$vectors)/sqrt(E$values)
                        E2 <- eigen(sqinv%*%t(X[smp,duse]-X[idout[smp],duse])%*%(X[smp,duse]-X[idout[smp],duse])%*%t(sqinv))
                        dr <- t(sqinv)%*%E2$vectors[,1:ndim]
                        dr <- dr/matrix(sqrt(colSums(dr^2)), nrow(dr), ncol(dr), byrow = TRUE)
                        if(OOB){
                          nns <- nn2(X[smp,duse]%*%dr, query = X[-smp,duse]%*%dr, k = nn, eps = eps)$nn.idx
                          pr_oob <- Y*0
                          for(nnn in 1:nn) pr_oob[-smp,] <- pr_oob[-smp,] + Y[smp[nns[,nnn]],]/nn
                          return(list(dr = dr, duse = duse, ixs = smp, pr_oob = pr_oob, evals = E2$values[1:ndim], oob = t(apply(pr_oob, 1, function(z) z==max(z)))))
                        }
                        else return(list(dr = dr, duse = duse, ixs = smp, evals = E2$values[1:ndim]))
                      })
    stopCluster(cluster)
    
  }
  else{
    mods <- lapply(1:nbag,
                      function(x){
                        smp <- sample(1:n, n_samp)
                        dsmp <- sample(1:d, ndim_samp)
                        duse <- which(idxs%in%dsmp)
                        
                        idin <- 1:n
                        idout <- 1:n
                        for(k in 1:max(y)){
                          relix <- smp[y[smp]==k]
                          relix2 <- smp[y[smp]!=k]
                          kuse2 <- min(nnp, length(relix2))
                          ny <- length(relix)
                          nny <- length(relix2)
                          if(ny==1){
                            idin[relix] <- relix
                            idout[relix] <- relix2[nn2(X[relix2,duse], query = matrix(X[relix,duse], nrow = 1), k = kuse2, eps = eps)$nn.idx[,kuse2]]
                          }
                          else if(ny>0){
                            kuse <- min(nnp+1,ny)
                            idin[relix] <- relix[nn2(X[relix,duse], k = kuse, eps = eps)$nn.idx[,kuse]]
                            if(nny > 1) idout[relix] <- relix2[nn2(X[relix2,duse], query = X[relix,duse], k = kuse2, eps = eps)$nn.idx[,kuse2]]
                            else if(nny==1) idout[relix] <- relix2[nn2(matrix(X[relix2,duse], nrow = 1), query = X[relix,duse], k = kuse2, eps = eps)$nn.idx[,kuse2]]
                          }
                        }
                        E <- eigen(t(X[smp,duse]-X[idin[smp],duse])%*%(X[smp,duse]-X[idin[smp],duse]) + diag(length(duse))*1e-5)
                        sqinv <- t(E$vectors)/sqrt(E$values)
                        E2 <- eigen(sqinv%*%t(X[smp,duse]-X[idout[smp],duse])%*%(X[smp,duse]-X[idout[smp],duse])%*%t(sqinv))
                        dr <- t(sqinv)%*%E2$vectors[,1:ndim]
                        dr <- dr/matrix(sqrt(colSums(dr^2)), nrow(dr), ncol(dr), byrow = TRUE)
                        if(OOB){
                          nns <- nn2(X[smp,duse]%*%dr, query = X[-smp,duse]%*%dr, k = nn, eps = eps)$nn.idx
                          pr_oob <- Y*0
                          for(nnn in 1:nn) pr_oob[-smp,] <- pr_oob[-smp,] + Y[smp[nns[,nnn]],]/nn
                          return(list(dr = dr, duse = duse, ixs = smp, pr_oob = pr_oob, evals = E2$values[1:ndim], oob = t(apply(pr_oob, 1, function(z) z==max(z)))))
                        }
                        else return(list(dr = dr, duse = duse, ixs = smp, evals = E2$values[1:ndim]))
                      })
  }
  if(vimp){
    vvals <- numeric(d)
    for(b in 1:nbag){
      Vb <- mods[[b]]$dr*matrix(sqrt(mods[[b]]$evals), nrow(mods[[b]]$dr), ncol(mods[[b]]$dr), byrow = TRUE)
      Delta <- rowSums(Vb^2)
      for(dd in 1:d){
        idi <- which(idxs[mods[[b]]$duse]==dd)
        if(length(idi) > 0) vvals[dd] <- vvals[dd] + mean(Delta[idi])
      }
    }
  }
  else vvals <- numeric(d)
  if(plotMat){
    X0 <- X*0
    for(b in 1:nbag) X0 <- X0 + X[,mods[[b]]$duse]%*%mods[[b]]$dr%*%t(mods[[b]]$dr)%*%t(sapply(mods[[b]]$duse, function(dd) (1:ncol(X))==dd))
    dr_all <- eigen(cov(X0))$vectors
  }
  else dr_all <- NULL
  if(OOB){
    Cl <- Y*0
    P <- Cl
    for(b in 1:nbag){
      Cl <- Cl + mods[[b]]$oob/nbag
      P <- P + mods[[b]]$pr_oob/nbag
    } 
    return(structure(list(mods = mods, probs_oob = P/rowSums(P), predicted.oob = ylevs[apply(Cl, 1, which.max)], oob.uncertainty = -apply(Cl/rowSums(Cl), 1, function(p) sum(p*log(p+1e-100))), ooberr = mean(apply(Cl, 1, which.max)!=y), nn = nn, convert = convert, idxs = idxs, X = X, x0 = x0, ylevs = ylevs, y = y, scales = scales, vimp = vvals/nbag, dr_all = dr_all), class = 'BOPNN'))
  }
  else return(structure(list(mods = mods, nn = nn, convert = convert, idxs = idxs, X = X, x0 = x0, ylevs = ylevs, y = y, scales = scales, vimp = vvals/nbag, dr_all = dr_all), class = 'BOPNN'))
}


add_plotMat <- function(model){
  X0 <- model$X*0
  for(b in 1:length(model$mods)) X0 <- X0 + model$X[,model$mods[[b]]$duse]%*%model$mods[[b]]$dr%*%t(model$mods[[b]]$dr)%*%t(sapply(model$mods[[b]]$duse, function(dd) (1:ncol(model$X))==dd))
  model$dr_all <- eigen(cov(X0))$vectors
  return(model)
}

add_vimps <- function(model){
  vvals <- numeric(ncol(model$x0))
  for(b in 1:length(model$mods)){
    Vb <- model$mods[[b]]$dr*matrix(sqrt(model$mods[[b]]$evals), nrow(model$mods[[b]]$dr), ncol(model$mods[[b]]$dr), byrow = TRUE)
    Delta <- rowSums(Vb^2)
    for(dd in 1:length(vvals)){
      idi <- which(model$idxs[model$mods[[b]]$duse]==dd)
      if(length(idi) > 0) vvals[dd] <- vvals[dd] + mean(Delta[idi])
    }
  }
  model$vimp <- vvals/length(model$mods)
  return(model)
}

predict.BOPNN <- function(model, test_df){
  ntest <- nrow(test_df)
  Xtest <- matrix(nrow = ntest, ncol = length(model$idxs))
  for(j in 1:ncol(model$x0)){
    if(j%in%model$convert){
      u <- unique(model$x0[,j])
      for(ui in 1:(length(u)-1)){
        Xtest[,which(model$idxs == j)[ui]] <- (test_df[,which(names(test_df)!='y')[j]] == u[ui])*1
      }
    }
    else Xtest[,which(model$idxs==j)] <- test_df[,which(names(test_df)!='y')[j]]
  }
  Xtest <- Xtest*matrix(model$scales, ntest, ncol(Xtest), byrow = TRUE)
  Y <- sapply(1:max(model$y), function(k) model$y==k)*1
  nbag <- length(model$mods)
  preds <- matrix(0, ntest, max(model$y))
  for(b in 1:nbag){
    p0 <- preds*0
    nns <- nn2(model$X[model$mods[[b]]$ixs,model$mods[[b]]$duse]%*%model$mods[[b]]$dr, query = Xtest[,model$mods[[b]]$duse]%*%model$mods[[b]]$dr, k = model$nn, eps = model$eps)$nn.idx
    for(nnn in 1:model$nn) p0 <- p0 + Y[model$mods[[b]]$ixs[nns[,nnn]],]/model$nn
    preds <- preds + t(apply(p0, 1, function(zz) zz==max(zz)))
  }
  model$ylevs[apply(preds, 1, which.max)]
}


plot.BOPNN <- function(model, new_data = NULL, labels = NULL, dims = c(1, 2), ...){
  if(is.null(model$dr_all)) stop('You need to first compute the ensemble projection. \n Use command "model <- add_plotMat(model)" and then "plot(model)".')
  if(is.null(new_data)) Xnew <- model$X
  else{
    Xnew <- matrix(nrow = nrow(new_data), ncol = length(model$idxs))
    for(j in 1:ncol(model$x0)){
      if(j%in%model$convert){
        u <- unique(model$x0[,j])
        for(ui in 1:(length(u)-1)){
          Xnew[,which(model$idxs == j)[ui]] <- (new_data[,which(names(new_data)!='y')[j]] == u[ui])*1
        }
      }
      else Xnew[,which(model$idxs==j)] <- new_data[,which(names(new_data)!='y')[j]]
    }
    Xnew <- Xnew*matrix(model$scales, nrow(Xnew), ncol(Xnew), byrow = TRUE)
  }
  if(is.null(labels)) plot(Xnew%*%model$dr_all[,dims], xlab = paste0("Discriminant Projection ", dims[1]), ylab = paste0("Discriminant Projection ", dims[2]), ...)
  else plot(Xnew%*%model$dr_all[,dims], xlab = paste0("Discriminant Projection ", dims[1]), ylab = paste0("Discriminant Projection ", dims[2]), col = labels, pch = labels, ...)
}


BOPNN_tune_parallel <- function(data_df, ranges = list(), beta_params = list(), nsamp = 30, eps = 0, cores = 1, case_weights = NULL){
  if(is.null(ranges$nn)) ranges$nn <- c(1, 5)
  if(is.null(ranges$dim_prop0)) ranges$dim_prop0 <- ncol(data_df)^(-.5)*c(1, min(10, ncol(data_df)^.5))
  if(is.null(ranges$dim_prop1)) ranges$dim_prop1 <- c(1/2, .99)
  if(is.null(ranges$n_prop)) ranges$n_prop <- c(.63, .63)
  
  if(is.null(beta_params$nn)) beta_params$nn <- c(1, 1)
  if(is.null(beta_params$dim_prop0)) beta_params$dim_prop0 <- c(1, 1)
  if(is.null(beta_params$dim_prop1)) beta_params$dim_prop1 <- c(1, 1)
  if(is.null(beta_params$n_prop)) beta_params$n_prop <- c(1, 1)
  
  cluster <- makeCluster(cores)
  clusterExport(cluster, c(ls(), 'nn2', 'BOPNN', 'one_hot_scale'), envir = environment())
  
  mods <- parLapply(cl = cluster, 1:nsamp,
                    function(x){
                      nn <- round(rbeta(1, beta_params$nn[1], beta_params$nn[2])*(ranges$nn[2]-ranges$nn[1])+ranges$nn[1])
                      dim_prop0 <- rbeta(1, beta_params$dim_prop0[1], beta_params$dim_prop0[2])*(ranges$dim_prop0[2]-ranges$dim_prop0[1])+ranges$dim_prop0[1]
                      dim_prop1 <- rbeta(1, beta_params$dim_prop1[1], beta_params$dim_prop1[2])*(ranges$dim_prop1[2]-ranges$dim_prop1[1])+ranges$dim_prop1[1]
                      n_prop <- rbeta(1, beta_params$n_prop[1], beta_params$n_prop[2])*(ranges$n_prop[2]-ranges$n_prop[1])+ranges$n_prop[1]
                      mod <- try(BOPNN(data_df, nn = nn, dim_prop1 = dim_prop1, dim_prop0 = dim_prop0, n_prop = n_prop, eps = eps, OOB = TRUE, cores = 1, case_weights = case_weights))
                      if(inherits(mod, 'try-error')) return(list(ooberr = 1))
                      else return(mod)
                    })
  stopCluster(cluster)
  prfs <- unlist(lapply(mods, function(l) l$ooberr))
  mods[[which.min(prfs)]]
}

