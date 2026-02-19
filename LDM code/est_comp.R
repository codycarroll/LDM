rm(list = ls())
source("./simulate_data.R")
library(fdapace)
library(fdasrvf)




getcomp = function(data, j, smooth = FALSE,bw=0.01){
  times = data$t %>% unique
  m = times %>% length
  comps = names(data)[2+(1:p)]
  compj = data[[comps[j]]] %>% matrix(nrow = m) %>% t 
  
  if(smooth == TRUE){
    compj = apply(compj, 1, function(x){Lwls1D(bw=bw, kernel_type = "gauss",
                                               xin = times, yin = x, xout = times)}) %>% t
  }
  return(compj)
}


est_components = function(data, method=NULL, smooth = FALSE, bw=0.01,
                         latent_method = c("heuristic", "global")){

  latent_method = match.arg(latent_method)

  if (is.null(method)) {
    method = "WFDA"
  }
  
  n = data$id %>% unique %>% length
  p = ncol(data) - 2
  comps = names(data)[2+(1:p)]
  domain <- range(data$t)
  workGrid <- ((data$t %>% unique) - domain[1])/(domain[2]-domain[1]) 
  id <- data$id %>% unique
  m = workGrid %>% length
  
  
  #amplitude factors
  A_ij.hat = list()
  for(i in 1:n){
    A_i = numeric(p)
    unwarped = list()
    for(j in 1:p){
      compj.mat = getcomp(data, j, smooth = smooth,bw=bw)
      #get max for ith subject and jth comp
      A_i[j] = max(compj.mat[i,])
    }
    A_ij.hat[[i]] = A_i
  }
  
  
  #estimate subject level warping functions and aligned mean curves per component
  H.hat.j = list()
  aligned.est = list()
  
  warp.obj = list()
  if(method == "WFDA"){
    for(j in 1:p){
      compj.mat = getcomp(data, j, smooth = smooth,bw=bw)
      scaledj.mat = (compj.mat/apply(compj.mat, 1, max)) #eliminate amp var
      
      wfda.in.j = MakeFPCAInputs(tVec = workGrid, yVec = scaledj.mat)
      #use subsetting if large sample to speed up
      if(n>=100){
        WFDA.j = WFDA(Ly = wfda.in.j$Ly , Lt = wfda.in.j$Lt,
                      optns = list(nknots = 4, subsetProp = 0.5))
      }else{
        WFDA.j = WFDA(Ly = wfda.in.j$Ly , Lt = wfda.in.j$Lt,
                      optns = list(nknots = 4))
      }

      H.hat.j[[j]] = WFDA.j$hInv
      aligned.est[[j]] = WFDA.j$aligned %>% colMeans
      warp.obj[[j]] = WFDA.j
    }
  }
  
  if(method == "SRVF"){
    for(j in 1:p){
      compj.mat = getcomp(data, j, smooth = smooth,bw=bw)
      scaledj.mat = (compj.mat/apply(compj.mat, 1, max)) 
      
      warp.j = time_warping(f = scaledj.mat %>% t, time = workGrid)
      #invert the warping functions to get H hat
      H.hat.j[[j]] = apply(warp.j$warping_functions , 2, function(x){approx(x = x, y = workGrid, xout = workGrid)$y}) %>% t
      aligned.est[[j]] = warp.j$fn %>% t %>% colMeans
      warp.obj[[j]] = warp.j
    }
  }
  
  
  #average over components
  H.hat = apply(simplify2array(H.hat.j), 1:2, mean)
  
  #aligned curves = component tempos with orig scale
  aligned.mat = aligned.est %>% unlist %>% matrix(ncol= p) %>% t
  
  #align component tempos to estimate component transports (Psi_j)
  aligned = time_warping(f = aligned.mat %>% t, time = workGrid)

  #estimate latent curve
  if(latent_method == "global"){
    # Global alignment method (Sec. 3.2): for each subject, randomly select one component,
    # normalize, align those n curves, and average the aligned curves as lambda hat
    J_i = sample(1:p, n, replace = TRUE)
    Z_star_mat = matrix(0, nrow = n, ncol = m)
    for(i in 1:n){
      compj_i = getcomp(data, J_i[i], smooth = smooth, bw = bw)
      Z_star_mat[i,] = compj_i[i,] / max(compj_i[i,])
    }
    if(method == "WFDA"){
      wfda.in.Z = MakeFPCAInputs(tVec = workGrid, yVec = Z_star_mat)
      if(n >= 100){
        WFDA.Z = WFDA(Ly = wfda.in.Z$Ly, Lt = wfda.in.Z$Lt,
                      optns = list(nknots = 4, subsetProp = 0.5))
      } else {
        WFDA.Z = WFDA(Ly = wfda.in.Z$Ly, Lt = wfda.in.Z$Lt,
                      optns = list(nknots = 4))
      }
      L.hat = WFDA.Z$aligned %>% colMeans
    } else {  # SRVF
      aligned.Z = time_warping(f = Z_star_mat %>% t, time = workGrid)
      L.hat = aligned.Z$fmean
    }
    L.hat = L.hat / max(L.hat)
  } else {
    # Heuristic: use the Frechet mean of the p aligned component tempos as lambda hat
    L.hat = aligned$fmean / max(aligned$fmean)
  }
  
  #component transports
  Psi.inv.hat = list()
  Psi.hat = list()
  for(j in 1:p){
    if(latent_method == "global"){
      # Align each component tempo to the global L.hat to get Psi_j^{-1}
      # pair_align_functions(f1, f2) finds gam s.t. f2(gam(t)) ~ f1(t)
      Psi.inv.hat[[j]] = pair_align_functions(f1 = L.hat,
                                              f2 = aligned.est[[j]],
                                              time = workGrid)$gam
    } else {
      Psi.inv.hat[[j]] = aligned$warping_functions[,j] #psi hat inverses
    }
    Psi.hat[[j]] = approx(x = Psi.inv.hat[[j]], y = workGrid, xout = workGrid)$y
  }
  
  #mixed effect warps 
  G_ij.hat = list()
  G_ij_inv.hat = list()
  
  for(i in 1:n){
    G_i.hat = list()
    G_i_inv.hat = list()
    
    for(j in 1:p){
      G_i.hat[[j]] = approx(x = workGrid, y = Psi.hat[[j]], xout = H.hat[i,])$y
      
      #line below does no dim reduction
      #G_i.hat[[j]] = approx(x = workGrid, y = Psi.hat[[j]], xout = H.hat.j[[j]][i,])$y
      G_i_inv.hat[[j]] = approx(x = G_i.hat[[j]], y = workGrid, xout = workGrid)$y
      
    }
    G_ij.hat[[i]] = G_i.hat
    G_ij_inv.hat[[i]] = G_i_inv.hat
    
  }
  

  #Xij hat
  fit_ij = function(i, j){
    LPH.hat = approx(x = workGrid, y = L.hat, xout = G_ij.hat[[i]][[j]])$y
    X.hat = A_ij.hat[[i]][j] * LPH.hat
    return(X.hat)
  }
  
  
  X.hat = list()
  for(i in 1:n){
    X.hat_df = data.frame(id = id[i],
                          t = workGrid)
    for(j in 1:p){
      newcolname = comps[j]
      X.hat_df = X.hat_df %>% 
        mutate(!!sym(newcolname) := fit_ij(i,j))
    }
    X.hat[[i]] = X.hat_df
  }
  
  
  #XCTs
  pairs = t(combn(p,2))
  L = nrow(pairs)
  T_jk.hat = list()
  
  for(l in 1:L){
    j = pairs[l,1]
    k = pairs[l,2]
    
    T_jk.hat[[l]] = approx(x = workGrid, y = Psi.inv.hat[[j]], xout = Psi.hat[[k]])$y
    T_jk.hat[[l]][m] = 1 #ensure numerical approx retains cdf structure
  }
  
  
  
  
  
  output = list(
    id = id,
    workGrid = workGrid,
    A_ij.hat = A_ij.hat,
    L.hat = L.hat,
    Psi.hat = Psi.hat,
    Psi.inv.hat = Psi.inv.hat,
    H.hat = H.hat, 
    G_ij.hat = G_ij.hat,
    X.hat = X.hat, 
    T_jk.hat = T_jk.hat,
    warp.obj = warp.obj
  )
  
  return(output)
}



#AUC decomposition

decomposition=function(dataset, tcol, varcol){
  
  #begin fn; 
  names(dataset)[tcol] <- c("t")
  varnames <- names(dataset)[varcol]
  M <- length(unique(dataset$t))
  domain <- range(dataset$t)
  workGrid <- ((dataset$t %>% unique) - domain[1])/(domain[2]-domain[1]) 
  id <- dataset$id %>% unique
  
  #split by id and define n & p 
  Lmods <- dlply(dataset, .(id), function(x)return(x))
  n <- length(Lmods)
  p <- length(varcol)
  
  
  afq_decomp=function(Lmods){
    
    
    #initialize outputs
    a_ik = list()
    F_ik = list()
    Q_ik = list()
    output=list()
    
    for(i in 1:n){
      temp=Lmods[[i]]
      Qtemp=temp 
      cols=ncol(temp)
      a_i = numeric(p)
      
      for(j in 1:p){
        a_i[j]=fastTrapz(x = workGrid, y = temp[,j+2])
        temp[,cols+j]=cumtrapz(x = workGrid, y = temp[,j+2])/a_i[j]
        Qtemp[,cols+j]=approx(x = temp[,cols+j], y = workGrid, xout = seq(0,1,length.out = M))$y
        Qtemp[M,cols+j] = 1
      }
      
      a_ik[[i]]=a_i
      F_ik[[i]]=temp[,c(1,2,(cols+1):(cols+p))]
      Q_ik[[i]]=Qtemp[,c(1,2,(cols+1):(cols+p))]
      
      names(a_ik[[i]])=varnames
      names(F_ik[[i]])[varcol]=varnames
      names(Q_ik[[i]])[varcol]=varnames
    }
    
    output=list(a_ik=a_ik, F_ik=F_ik, Q_ik=Q_ik)
    
    return(output)
  }
  afq.obj  = afq_decomp(Lmods = Lmods)
  
  gamma_est=function(afq.obj){
    #p = afq.obj$a_ik[[1]] %>% length
    #M = afq.obj$F_ik[[1]]$t %>% length
    gamma_hat = list()
    gamma_inv = list()
    gammaprime_hat = list()
    
    for(j in 1:p){
      
      Qmatj = matrix(0, nrow = M, ncol = n)
      
      for( i in 1:n ){
        Qmatj[,i] = afq.obj$Q_ik[[i]][,2+j]
      }
      
      gamma_inv[[j]] = Qmatj %>% rowMeans
      gamma_hat[[j]] = approx(x = Qmatj %>% rowMeans, y = workGrid, xout = workGrid)$y
      gammaprime_hat[[j]] = Lwls1D(xin = workGrid, yin = gamma_hat[[j]], xout = workGrid, 
                                   bw = 0.03, nder = 1, kernel_type = "gauss")
      
    }
    
    gamma = list(hat = gamma_hat, inv = gamma_inv, primehat = gammaprime_hat)
    return(gamma)
  }
  gamma = gamma_est(afq.obj = afq.obj)
  
  latent = function(gamma){
    latentinv = apply(simplify2array(gamma$inv), 1, mean)
    latent = approx(x = latentinv, y = workGrid, xout = workGrid)$y
    return(latent)
  }
  
  latent = latent(gamma)
  
  psi_est = function(gamma){
    psi = list()
    latentinv = apply(simplify2array(gamma$inv), 1, mean)
    for(j in 1:p){
      psi[[j]] = approx(x = workGrid, y = latentinv, xout = gamma$hat[[j]])$y
    }
    
    return(psi)
  }
  
  psi = psi_est(gamma)
  
  Hik_est = function(afq.obj, gamma){
    #p = afq.obj$a_ik[[1]] %>% length
    Fiks = afq.obj$F_ik
    Hiks_hat = Fiks
    Hiks_inv = Fiks
    for(i in 1:n){
      temp=Hiks_hat[[i]]
      tempinv=Hiks_inv[[i]]
      Fi=Fiks[[i]]
      for(j in 1:p){
        Ginvj=gamma$inv[[j]]
        Fik=Fi[,j+2]
        temp[,j+2]=approx(x = workGrid, y = Ginvj, xout = Fik)$y
        temp[M,j+2]=1
        
        tempinv[,j+2]=approx(x = temp[,j+2], y = workGrid, xout = workGrid)$y
        tempinv[M,j+2]=1
      }
      Hiks_hat[[i]] = temp
      Hiks_inv[[i]] = tempinv
    }
    Hiks = list(hat = Hiks_hat, inv = Hiks_inv)
    return(Hiks)
  }
  Hiks = Hik_est(afq.obj = afq.obj, gamma = gamma)
  
  Hibar_est=function(Hiks){
    Hibar = list()
    Hibarinv = list() 
    
    for( i in 1:n ){
      
      Hibarinv[[i]] = rowMeans(Hiks$inv[[i]][,varcol])
      Hibar[[i]] = approx(x = rowMeans(Hiks$inv[[i]][,varcol]), y = workGrid, xout = workGrid)$y
      
    }
    
    Hi = list(bar = Hibar, barinv = Hibarinv)
    return(Hi)
  }
  Hi = Hibar_est(Hiks)
  
  Tijk_est = function(afq.obj){
    xcts = list()
    pairs = t(combn(p,2))
    L = nrow(pairs)
    Tjk_hat = list()
    Tjk_inv = list()
    
    
    for(i in 1:n){
      Tijk = matrix(0, nrow = M, ncol = L)
      for(l in 1:L){
        j = pairs[l,1]
        k = pairs[l,2]
        Tijk[,l] = approx(x = workGrid, y = afq.obj$Q_ik[[i]][,varcol[j]], xout = afq.obj$F_ik[[i]][,varcol[k]])$y
        Tijk[M,l] = 1
      }
      xcts[[i]] = Tijk
      colnames(xcts[[i]]) = paste0(pairs[,1],"to",pairs[,2])
    }
    Tjk_hat = apply(simplify2array(xcts), 1:2, mean)
    
    Tijk.obj = list(Tijk = xcts, Tjk_hat = Tjk_hat)
    return(Tijk.obj)
  }
  
  Tijk.obj = Tijk_est(afq.obj)
  
  fit_ij = function(i, j, bw, scaled = FALSE){
    Xi_ij = approx(x = workGrid, y = psi[[j]], xout = Hi$bar[[i]])$y
    LPHhat = approx(x = workGrid, y = latent, xout = Xi_ij)$y
    fhat = Lwls1D(bw = bw, kernel_type = "gauss", nder = 1,
                  xin = workGrid, yin = LPHhat, xout = workGrid)
    if(scaled == TRUE){
      Xhat = fhat         
    }else{
      Xhat = afq.obj$a_ik[[i]][j] * fhat
    }
    return(Xhat)
  }
  
  fitall = function(scaled = FALSE){
    Xhat = list()
    for(i in 1:n){
      Xhat_df = data.frame(id = id[i],
                           t = workGrid)
      for(j in 1:p){
        newcolname = varnames[j]
        Xhat_df = Xhat_df %>% 
          mutate(!!sym(newcolname) := fit_ij(i,j, bw = 0.03, scaled))
      }
      Xhat[[i]] = Xhat_df
    }
    
    return(bind_rows(Xhat))
  }
  Xhat = fitall()
  Xhatscaled = fitall(scaled = TRUE)
  
  output = list(id = id,
                aik = afq.obj$a_ik, 
                Fik = afq.obj$F_ik, 
                Qik = afq.obj$Q_ik,
                Tijk = Tijk.obj$Tijk,
                Tjk_hat = Tijk.obj$Tjk_hat,
                Gammahat = gamma$hat,
                Gammainv = gamma$inv,
                Gammaprimehat = gamma$primehat, 
                Lambda = latent,
                psi = psi,
                Hik = Hiks$hat,
                Hikinv = Hiks$inv,
                Hibar = Hi$bar,
                Hibarinv = Hi$barinv,
                Xhat = Xhat, 
                Xhatscaled = Xhatscaled,
                workGrid = workGrid
  )
  
  return(output)
}



###shifts 

padcurves = function(data, cushion){
  #define before/after grids
  beforeworkGrid = seq(-cushion,-0.01,by= 0.01)
  workGrid = seq(0,1,by= 0.01)
  afterworkGrid = seq(1.01,1+cushion,by= 0.01)
  paddedworkGrid = c(beforeworkGrid, workGrid, afterworkGrid)
  
  #before/after data frames
  #add extra copies of first/last timepoint to allow for shifting
  
  before = data %>% slice(rep(1, each=length(beforeworkGrid )))
  before$t = beforeworkGrid
  after= data %>% slice(rep(m, each=length(afterworkGrid )))
  after$t = afterworkGrid
  
  #combine
  padded = bind_rows(before, data, after)
  padded
}

normalize = function(data){
  for(j in 1:p){
    data[,2+j] = data[,2+j]/max(data[,2+j])
  }
  data
}

fastTrapz <- function(x, y){
  
  return( fdapace:::trapzRcpp( X = x, Y = y ) )
  
}

L2dist<-function(c,dataset,pair, I=c(0,1)){
  
  #normalize curves
  dataset = ddply(dataset, .(id), normalize)
  
  # start = Sys.time()
  #get unshifted subsets
  spurtmods=subset(dataset,dataset$t>= I[1] & dataset$t < I[2])
  Lmods <- dlply(spurtmods, .(id), function(x)return(x))
  ##get shifted subsets
  shiftmods=subset(dataset, round(dataset$t+c,2) >= I[1] & round(dataset$t+c,2) < I[2])
  shiftLmods = dlply(shiftmods, .(id), function(x)return(x))
  #pair -> j,k 
  j = pair[1]
  k = pair[2]
  ##compute L2 distance
  L2=numeric(length(Lmods))
  for(i in 1:length(Lmods)){
    L2[i]=fastTrapz(x = Lmods[[i]]$t, y = ((Lmods[[i]])[,2+j]-(shiftLmods[[i]])[,2+k])^2)
  }
  distance=sum(L2)
  distance
  # timing = Sys.time() - start
  # ret = list(distance = distance, timing = timing)
  # return(ret)
}


pairWFDA=function(pair, paddata){
  mods=data
  s_jk=round(optim(par = 0,  fn = L2dist, method = "Brent", dataset=paddata, pair = pair, I=c(0,1), lower = -.4, upper = .4)$par,4) 
  #s_jk=round(optim(par = 0,  fn = L2dist, method = "CG", dataset=paddata, pair = pair, I=c(0,1))$par,4) 
  cat("The optimal pairwise shift between modalities",paste0(pair[1]),"and",paste0(pair[2]),"over all subjects in the dataset is", paste0(s_jk), "\n")
  return(s_jk)
}



smoothcomps = function(data,bw=0.01){
  for(j in 1:p){
    data[,j+2] = Lwls1D(bw=bw, kernel_type = "gauss", xin = data$t, yin = data[,j+2], xout = data$t)
  }
  data
}




est_shifts=function(data, bw=0.01){
  
  #PAD
  paddata = ddply(data, .(id), function(df){padcurves(df, cushion = 1)})
  
  #SMOOTH
  paddata = ddply(paddata, .(id), function(df){smoothcomps(df, bw=bw)})
  
  
  #pairwise shifts
  allpairs=t(combn(p, 2))
  optShifts=apply(allpairs, 1,  FUN=pairWFDA, paddata=paddata)
  s=c(optShifts,0)
  A1=matrix(0,nrow = (p*(p-1)/2), ncol = p)
  for(i in 1:(p*(p-1)/2)){
    A1[i,allpairs[i,]]=c(1,-1)
  }
  A=rbind(A1,1)
  
  #global shifts
  theta=round(1/p*t(A)%*%s,3)
  theta
  
  psi_est_shift=list()
  for(j in 1:p){
    psi_est_shift[[j]] = workGrid-theta[j]
  }
  
  output = list(tau=s,
                theta=theta,
                psi_est_shift=psi_est_shift)
  
  return(output)
}

#check that it works... yes
set.seed(1234321)
simdata = simulate_data(n = 100, sigma_amp = 5, sigma_warp = .5, sigma_error = 5)
data = simdata$simdata
est_shifts(data = data)

