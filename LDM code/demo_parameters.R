
library(reshape2)
library(ggplot2)
library(magrittr)
set.seed(1234321)
simdata = simulate_data(n = 100, 
                        sigma_amp = 5, 
                        sigma_warp = .5,
                        sigma_error = 5)

sa = c(1,2,5)
sw = c(0.01,0.5,1, 2)
se = c(0,1,5,10)

simdatanowarp = matrix(0, nrow = m, ncol = p+1)

simdatanowarp[,1] = t


for(j in 1:p){
  simdatanowarp[,j+1] = ameans[j] * makecurvenowarp(j)
}

simdatanowarp = as.data.frame(simdatanowarp)
names(simdatanowarp) = c("t", "A","B","C","D")

ggdata=melt(simdata$simdata, id.vars = c("id","t"))
ggdatanowarp=melt(simdatanowarp, id.vars = "t")
ggcolors = c("red", "darkorange", "turquoise2", "blue")


ggsim = ggplot()+
  geom_line(data = ggdata,
            aes(x=t, y = value, group=id, color = variable),
            size=.6, alpha = 0.2) +
  geom_line(data = ggdatanowarp,
            aes(x=t, y = value),
            color = "black", size = .75) +
  xlab("t") +
  ylab("")+
  facet_wrap(~variable, ncol = 1) +
  scale_color_manual(values=ggcolors) +
  theme(legend.position = "none",
        #axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0.01,0.01))
ggsim


ggsimindiv = ggplot()+
  geom_line(data = ggdata %>% subset(id==1),
            aes(x=t, y = value, group=variable, color = variable),
            size=.6, alpha = 0.2) +
  # geom_line(data = ggdatanowarp,
  #           aes(x=t, y = value),
  #           color = "black") +
  xlab("t") +
  ylab("")+
  #facet_wrap(~variable, ncol = 1) +
  scale_color_manual(values=ggcolors) +
  theme(legend.position = "none",
        #axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0.01,0.01))
ggsimindiv

