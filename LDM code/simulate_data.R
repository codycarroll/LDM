


library(ggplot2); theme_set(theme_bw());theme_update(plot.title = element_text(hjust = 0.5))

library(fdapace)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(pracma)
library(plyr)
library(readxl)
library(MASS)
library(ggpubr)
library(fdadensity)
library(seriation)
library(officer)
library(rvg)
#source(file = "../../subroutines.R")



###### SET WORKGRID AND PARAMETERS
m = 101
p = 4
n = 100
t = workGrid = seq(0, 1, length.out = m)

Lshape = function(t){
  20-(5)*cos(4*pi*t)+(3)*sin(pi*t^2)+15*(t)^2 
}

ameans = c(100,100,100,100) 

L = function(t){
  (20-(5)*cos(4*pi*t)+(3)*sin(pi*t^2)+15*(t)^2) / max(Lshape(t))
}


psi_true = list()
psi_inv = list()


alpha = c(1,2)
beta = c(1/2,2)

lambda = .5

psi = list() 
psiinv = list()

for(j in 1:(p/2)){
  t = seq(0, 1, length.out = m)
  psi[[j]] = lambda*pbeta(t, alpha[j], beta[j]) + (1-lambda)*t
  psiinv[[j]] = approx(psi[[j]], t, t)$y
  
  psiinv[[j+2]] = 2*t - psiinv[[j]]
  psi[[j+2]] = approx(psiinv[[j+2]], t, t)$y
  
}


psi_true=list()
psi_inv=list()

order = c(1,3,2,4)
for(j in 1:p){
  psi_true[[j]] = psi[[order[j]]]
  psi_inv[[j]] = psiinv[[order[j]]]
  
}

Tjk_true=list()
l=combn(p,2) %>% t
for(i in 1:nrow(l)){
  j=l[i,1]
  k=l[i,2]
  Tjk_true[[i]]=approx(x = t, y = psi_inv[[j]], xout = psi_true[[k]])$y
}



n=50; sigma_amp=5; sigma_warp=.5; sigma_error = 0; sigma_dist = 0.3
simulate_data = function(n, sigma_amp, sigma_warp, sigma_error = 0, sigma_dist = 0){
  
  #generate warping
  hinvmat = matrix(0, nrow = n, ncol = m)
  hmat = matrix(0, nrow = n, ncol = m)
  z = rnorm(n = n, 0, sigma_warp^2)
  for(i in 1:n){
    hinvmat[i,] = (exp(t*z[i])-1)/(exp(z[i])-1)
    hmat[i,] =  approx(hinvmat[i,], t, t)$y
  }
  
  #generate random amplitudes
  amat = matrix(0, nrow = n , ncol = p)
  for(j in 1:p){
    amat[,j] = pmax(rnorm(n, mean = ameans[j], sd = sigma_amp^2), 5)
  }
  
  #function to make curve shape w warping
  makecurve = function(i,j){
    if(sigma_dist==0){
      dist=t
    }else{
      d = rnorm(n = 1, 0, sigma_dist^2)
      dinv = (exp(t*d)-1)/(exp(d)-1)
      dist =  approx(dinv, t, t)$y
    }
    hd = approx(t,  hmat[i,], xout = dist)$y
    PoH = approx(t, psi_true[[j]], xout = hd)$y
    L(PoH)
  }
  
  
  simdata = matrix(0, nrow = n*m, ncol = p+2)
  
  simdata[,1] = rep(1:n,each = m)
  simdata[,2] = rep(t, times = n)
  
  for(j in 1:p){
    complist = list()
    for(i in 1:n){
      complist[[i]] =  amat[i, j] * makecurve(i, j) + rnorm(n = m, 
                                                            mean = 0, 
                                                            sd = sigma_error)
      complist[[i]] = pmax(complist[[i]], .1) #force positive; avoid measurement error leading to negatives
    }
    simdata[,j+2] = unlist(complist) 
  }
  
  simdata = as.data.frame(simdata)
  names(simdata) = c("id", "t", "A","B","C","D")
  
  return(list(simdata = simdata, hmat = hmat, amat = amat))
}

#comp barycenter curves with no warping 
makecurvenowarp = function(j){
  PoI = approx(t, psi_true[[j]], xout = t)$y
  L(PoI)
}




###IMAGE 
# 
# 
# create_pptx <- function(plt = last_plot(), path = file.choose(), w, h){
#   if(!file.exists(path)) {
#     out <- read_pptx()
#   } else {
#     out <- read_pptx(path)
#   }
#   out %>%
#     add_slide(layout = "Title and Content", master = "Office Theme") %>%
#     ph_with(value = dml(ggobj = plt), location = ph_location(width = w, height = h)) %>%
#     print(target = path)
# }
# 
# filefolder = "/Users/codycarroll/repos/xctransport/text/ppt/simfig_warps.pptx"
# 
# 
# set.seed(1234321)
# simdata = simulate_data(n = 100,sigma_amp = 4, sigma_warp = 0.01, sigma_error = 0, sigma_dist = 0.01)
# # 
# simdatanowarp = matrix(0, nrow = m, ncol = p+1)
# 
# simdatanowarp[,1] = t
# 
# for(j in 1:p){
#   simdatanowarp[,j+1] = ameans[j] * makecurvenowarp(j)
# }
# 
# simdatanowarp = as.data.frame(simdatanowarp)
# names(simdatanowarp) = c("t", "A","B","C","D")
# # 
# ggdata=melt(simdata$simdata, id.vars = c("id","t"))
# ggdatanowarp=melt(simdatanowarp, id.vars = "t")
# ggcolors = c("red", "darkorange", "turquoise2", "blue")
# 
# 
# ggsim = ggplot()+
#   geom_line(data = ggdata,
#             aes(x=t, y = value, group=id, color = variable),
#             size=.6, alpha = 0.2) +
#   geom_line(data = ggdatanowarp,
#             aes(x=t, y = value),
#             color = "black") +
#   xlab("t") +
#   ylab("")+
#   facet_wrap(~variable, ncol = 1) +
#   scale_color_manual(values=ggcolors) +
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank()) +
#   scale_x_continuous(expand = c(0.01,0.01))
# ggsim
# 
# 
# 
# latentdf = data.frame(t=t, L=L(t))
# 
# 
# gglatent = ggplot() +
#   geom_line(data = latentdf, aes(x = t, y = L) , color = "black", size = 1.25)+
#   ylab(" ") +
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank())+
#   #ggtitle("Latent Curve") +
#   scale_x_continuous(expand = c(0.01,0.01))
# 
# 
# #gamma
# 
# simdatanowarpscale = matrix(0, nrow = m, ncol = p+1)
# 
# simdatanowarpscale[,1] = t
# 
# for(j in 1:p){
#   simdatanowarpscale[,j+1] = makecurvenowarp(j)
# }
# 
# simdatanowarpscale = as.data.frame(simdatanowarpscale)
# names(simdatanowarpscale) = c("t", "A","B","C","D")
# ggdatanowarpscale=melt(simdatanowarpscale, id.vars = "t")
# 
# gggamma = ggplot(ggdatanowarpscale, aes(x=t, y = value, color = variable))+
#   xlab("t") +
#   ylab("")+
#   ylim(c(0.5,1.5)) +
#   geom_line(size=1.15, alpha = .8) +
#   facet_wrap(~variable) +
#   scale_color_manual(values=ggcolors) +
#   facet_wrap(~variable, ncol = 1) +
#   theme(legend.position = "none",
#   axis.text.y = element_blank(),
#   axis.ticks = element_blank())
# #psi 
# 
# simdata = simulate_data(n = 100,
#                         sigma_amp=5,
#                         sigma_warp=.5,
#                         sigma_error=0,
#                         sigma_dist=.5)
# 
# #set grid, id, and components 
# workGrid = seq(0,1,length.out = m)
# data = simdata$simdata
# comps = LETTERS[1:p]
# id = data$id
# 
# #estimate the components, XCT
# test = est_components(data = simdata$simdata, smooth = TRUE)
# 
# psi_df = data.frame(t = t,
#                     `1` = psi_true[[1]],
#                     `2` = psi_true[[2]],
#                     `3`= psi_true[[3]],
#                     `4`= psi_true[[4]], check.names=FALSE)
# 
# psi_df = psi_df %>%
#   melt(id.vars = 't')
# 
# 
# ggpsi = ggplot(data = psi_df,
#                aes(x = t, y = value,
#                    group = variable, color = variable)) +
#   geom_line(size = 1.15) +
#   ylab(" ") +
#   scale_color_manual(values=ggcolors) +
#   theme(legend.position = "none")
# ggpsi
# 
# T_df = data.frame(t = t,
#                   AB = Tjk_true[[1]],
#                   AC = Tjk_true[[2]],
#                   AD= Tjk_true[[3]],
#                   BC= Tjk_true[[4]],
#                   BD= Tjk_true[[5]],
#                   CD= Tjk_true[[6]],
#                   check.names=FALSE)
# 
# T_df = T_df %>%
#   melt(id.vars = 't')
# 
# pairs = t(combn(p,2))
# L = nrow(pairs)
# T_jk.hat = list()
# 
# for(l in 1:L){
#   j = pairs[l,1]
#   k = pairs[l,2]
#   
#   T_jk.hat[[l]] = approx(x = workGrid, y = Psi.inv.hat[[k]], xout = Psi.hat[[j]])$y
#   T_jk.hat[[l]][m] = 1 #ensure numerical approx retains cdf structure
# }
# 
# 
# xcts_warp_df = data.frame(t = workGrid,
#                           AB = T_jk.hat[[1]],
#                           AC = T_jk.hat[[2]],
#                           AD = T_jk.hat[[3]],
#                           BC = T_jk.hat[[4]],
#                           BD = T_jk.hat[[5]],
#                           CD = T_jk.hat[[6]]
# )
# 
# 
# 
# ggT = ggplot() +
#     geom_line(data = T_df,
#               aes(x = t, y = value,
#                   group = variable, color = variable))+
#     geom_line(data = xcts_warp_df %>% melt(id.vars = 't'),
#               mapping = aes(x = t, y = value,
#                             group = variable, color = variable), linetype = "dashed")+
#     
#   geom_line(size = 1.15) +
#   ylab(" ") +
#   #scale_color_manual(values=ggcolors) +
#   theme(legend.position = "none")
# ggT
#   
# 
# #H 
# hdf = data.frame(id = as.factor(rep(1:n, each = m)), t = rep(t, times = n), warp = c(t(simdata$hmat)))
# 
# ggwarp = ggplot() +
#   geom_line(data = hdf, aes(x = t, y = warp, group = id),
#             size=.6, alpha = 0.3, color = "black") +
#   ylab("") +
#   theme(legend.position = "none",
#         axis.text.y = element_blank(),
#         axis.ticks = element_blank())+
#   scale_x_continuous(expand = c(0.01,0.01)) 
# 
# ggwarp
# 
# 
# 
# create_pptx(plt = gglatent, path = filefolder, w= 2.5, h = 2)
# create_pptx(plt = ggpsi, path = filefolder, w= 2.5, h = 2)
# create_pptx(plt = gggamma, path = filefolder, w= 2.25, h = 6)
# create_pptx(plt = ggwarp, path = filefolder, w= 2.5, h = 2)
# create_pptx(plt = ggsim, path = filefolder, w = 2.25, h = 6)
# 
# 
#T_24
# approx(x = t, y = psi_inv[[1]], xout = psi_true[[2]])$y %>% plot
# 
# plot(t, psi_true[[2]])


