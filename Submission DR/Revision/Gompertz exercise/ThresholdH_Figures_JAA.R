################################################################################################
# Threshold age of H
################################################################################################

library(tidyverse)
library(readr)
library(data.table)
setwd("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Threshold of H/R")

theme_set(
  theme_bw(base_size = 20)+
    theme(axis.title = element_text(size = 20),
          axis.ticks = element_line(color = "grey", size = .5),
          axis.ticks.length = unit(.5, "line"),
          # panel.grid = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y = element_text(vjust = 3),
          panel.border = element_rect(color = "grey90", fill = NA, size = .5),
          line = element_line(lineend = "round"))
)

# Function to calculate threshold age
age.H<- function(Age, ax, dx, lx, ex){ 
  
  ax<- ax
  dx<- dx /100000
  lx<- lx /100000
  ex<- ex
  Age <-  Age
  
  Age[Age == "110+"] <-110
  Age <- as.integer(as.character(Age))
  #e-dagger at age x
  ex.bar    <- ex +ax*(ex - c(ex[-1], ex[length(ex)]))
  ex.bar.dx <- ex.bar * dx
  e.dag.x   <-  rev(cumsum(rev(ex.bar.dx)))/ lx
  
  #Keyfitz's entropy
   Hx <- e.dag.x / ex
  
  #Cumulative hazard
   cumhaz <- -log(lx)
   
  #g(x)
   gx <- cumhaz + Hx - 1- Hx[1]
   
   #wx
   
   wx <- (dx / Hx) * gx
  
   #a.dagger
   a.dag <- ex * (1- cumhaz)
   
   #get threshold age of H
    
   return(data.frame(Age,ex, lx, e.dag.x, Hx, cumhaz, gx, wx, a.dag))}

get.a <- function(Age,gx,ex, e.dag.x, a.dag){
  e0 <-  ex[1]
  f1 <- approxfun(Age,gx,method = "linear",rule=2 )
  a.H <- uniroot(function(x) f1(x),c(0,110))$root
  
  
  f2 <- approxfun(Age,e.dag.x,method = "linear",rule=2 )
  f3 <- approxfun(Age, a.dag ,method = "linear",rule=2 )
  
  a.dagger <- uniroot(function(x) f2(x) - f3(x),c(0,110))$root
  
  result <-  data.frame(e0,a.H, a.dagger)
  return(result)
}



Dat <- data.table(read_table2("C:/Users/jmartinez/OneDrive - Syddansk Universitet/Threshold of H/Data/SWE.fltper_1x1.txt", skip = 3))

pr <- Dat[ ,age.H(Age = Age, ax = ax, dx = dx, lx = lx, ex = ex), by = list(Year)]


pr.a <-  subset(pr, Year >=1800)


# get thresholds for all years
as <- pr.a[ ,get.a(Age = Age, gx = gx, ex = ex,  e.dag.x = e.dag.x, a.dag =a.dag), by = list(Year)]

pr.t <- merge(pr.a, as, by = "Year")

ggplot(as)+
  geom_line(aes(Year, e0), size =1)+
  geom_line(aes(Year, a.H), colour= "red")+
  geom_line(aes(Year, a.dagger), colour= "blue")+
  theme_bw()+
  coord_cartesian(ylim = c(20,90))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(strip.background = element_rect(fill="white"))+
  theme(strip.background = element_rect(fill="white"))+
  theme(legend.title=element_blank(), 
        panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.spacing = unit(2, "lines"),
        axis.title = element_text( size = 20),
        axis.ticks = element_line(color = "grey", size = .5),
        axis.ticks.length = unit(.5, "line"),
        # panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(vjust = 3),
        panel.border = element_rect(color = "grey90", fill = NA, size = .5),
        line = element_line(lineend = "round"))+
  xlab("Year") +
  ylab("Threshold ages and life expectancy at birth (in years)")

ggsave("fig/2_Threshold_SWE.pdf", width = 6, height = 4, device = cairo_pdf)


pr.b <- subset(pr.t, Year==2005)

ggplot(pr.b)+
  geom_line(aes(Age, gx), colour = "blue")+
  geom_line(aes(Age, a.dag / 10), colour = "red")+
  geom_line(aes(Age, e.dag.x / 10), colour = "gray50")+
  geom_vline(xintercept = pr.b$a.H[1], linetype = "dotted", colour = "blue", size = 1)+ #threshold aH
  geom_vline(xintercept = pr.b$a.dagger[1], linetype = "dotted", colour = "red", size = 1)+ #threshold aH
  geom_vline(xintercept = pr.b$ex[1], linetype = "dashed", colour = "black")+ #threshold aH
  geom_hline(yintercept = 0, colour = "black")+
  theme_bw()+
  coord_cartesian(ylim = c(-1,8), xlim = c(0,110))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+
  theme(strip.background = element_rect(fill="white"))+
  theme(legend.title=element_blank(), 
        panel.background = element_blank(),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        panel.spacing = unit(2, "lines"),
        axis.title = element_text( size = 20),
        axis.ticks = element_line(color = "grey", size = .5),
        axis.ticks.length = unit(.5, "line"),
        # panel.grid = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(vjust = 3),
        panel.border = element_rect(color = "grey90", fill = NA, size = .5),
        line = element_line(lineend = "round"))+
  xlab("Age") +
  ylab("Functions to determine threshold ages")

ggsave("fig/1_Ages_ITA.pdf", width = 6, height = 4, device = cairo_pdf)



#### Compare entropy with edag ####

source("FUN.R")
Dat$Age <- as.character(Dat$Age)
Dat$Age[Dat$Age == "110+"] <-110
Dat$Age <- as.numeric(Dat$Age)

Var <- as.data.frame(Dat[,Var.fun(x = Age, ax = ax, dx = dx, lx = lx, ex = ex), by = list(Year)])

library(Hmisc)
library(ggpmisc)

ggplot(Var, aes(Year, e.dagger))+
  geom_line()+
  stat_summary(fun.data=mean_cl_normal)+
  geom_smooth(data = subset(Var, Year<= 1900),method='lm',formula=y~x, colour = "red")+
  geom_smooth(data = subset(Var, Year>= 1900 & Year<=1950),method='lm',formula=y~x, colour = "blue")+
  geom_smooth(data = subset(Var, Year>= 1950),method='lm',formula=y~x, colour = "green")+
  stat_poly_eq(data = subset(Var, Year<= 1900),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "red")+
  stat_poly_eq(data = subset(Var, Year>= 1900 & Year<=1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "blue")+
  stat_poly_eq(data = subset(Var, Year>= 1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "green")
  



ggplot(Var, aes(Year, H))+
  geom_line()+
  stat_summary(fun.data=mean_cl_normal)+
  geom_smooth(data = subset(Var, Year<= 1900),method='lm',formula=y~x, colour = "red")+
  geom_smooth(data = subset(Var, Year>= 1900 & Year<=1950),method='lm',formula=y~x, colour = "blue")+
  geom_smooth(data = subset(Var, Year>= 1950),method='lm',formula=y~x, colour = "darkgreen")+
  stat_poly_eq(data = subset(Var, Year<= 1900),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "red")+
  stat_poly_eq(data = subset(Var, Year>= 1900 & Year<=1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "blue")+
  stat_poly_eq(data = subset(Var, Year>= 1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "darkgreen")+
  ylab("Lifetable entropy")


ggplot(Var, aes(Year, e.dagger))+
  geom_line()+
  stat_summary(fun.data=mean_cl_normal)+
  geom_smooth(data = subset(Var, Year<= 1900),method='lm',formula=y~x, colour = "red")+
  geom_smooth(data = subset(Var, Year>= 1900 & Year<=1950),method='lm',formula=y~x, colour = "blue")+
  geom_smooth(data = subset(Var, Year>= 1950),method='lm',formula=y~x, colour = "darkgreen")+
  stat_poly_eq(data = subset(Var, Year<= 1900),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "red")+
  stat_poly_eq(data = subset(Var, Year>= 1900 & Year<=1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "blue")+
  stat_poly_eq(data = subset(Var, Year>= 1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "darkgreen")+
  ylab("e-dagger")

ggplot(Var, aes(Year, er))+
  geom_line()+
  stat_summary(fun.data=mean_cl_normal)+
  geom_smooth(data = subset(Var, Year<= 1900),method='lm',formula=y~x, colour = "red")+
  geom_smooth(data = subset(Var, Year>= 1900 & Year<=1950),method='lm',formula=y~x, colour = "blue")+
  geom_smooth(data = subset(Var, Year>= 1950),method='lm',formula=y~x, colour = "darkgreen")+
  stat_poly_eq(data = subset(Var, Year<= 1900),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "red")+
  stat_poly_eq(data = subset(Var, Year>= 1900 & Year<=1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "blue")+
  stat_poly_eq(data = subset(Var, Year>= 1950),parse=T, aes(label = ..eq.label..), formula=y~x, colour = "darkgreen")+
  ylab("Life expectancy at birth")

