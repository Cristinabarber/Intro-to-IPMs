# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1 - Load up packages and read in data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library("ggplot2")
library("plyr")
library("dplyr")
library("bindr")
library("stringr")
library(MASS)
library("bbmle")

#tables
BN<-read.table(file="C:/Users/cristinabarberal/Documents/Courses/DB_brazilnut/DB_brazilnut.csv",header=T,sep=",")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2 - Fit the vital rate functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Growth linear model

grow_data <- subset(BN, !is.na(BN$DBH2))
mod_grow <- lm(grow_data$DBH2 ~ grow_data$DBH+grow_data$Logging, data = grow_data)
summary(mod_grow)


#Survival binomial model
surv_data <- subset(BN, !is.na(BN$DBH))
mod_surv <- glm(surv_data$Survival ~ surv_data$DBH + surv_data$Logging, family = binomial, data = surv_data)
summary(mod_surv)

#Recruitment:

fruit_data <- subset(BN, BN$Survival==1 & BN$Reproductive==1)
mod_fruit <- glm(fruit_data$Fruits ~ fruit_data$DBH + fruit_data$BN_harvest + fruit_data$Logging, family = "poisson", data = fruit_data)
summary(mod_fruit)


recruit_data <- subset(BN, BN$Survival==1)
mod_recruit <- glm(recruit_data$Recruitment ~ recruit_data$DBH2 + recruit_data$BN_harvest, family = binomial, data = recruit_data)
summary(mod_recruit)

reproduct_data <- subset(BN, BN$Survival==1 & BN$DBH>40)
mod_reproduct <- glm(reproduct_data$Reproductive ~ reproduct_data$DBH + reproduct_data$Logging+ reproduct_data$BN_harvest, family = binomial, data = reproduct_data)
summary(mod_reproduct)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3 - Store the parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # survival
  surv      = as.numeric(coef(mod_surv))
  # growth 
  grow      =  as.numeric(coef(mod_grow))
  grow_sd   =  as.numeric(summary(mod_grow)$sigma)
  # reproduce or not
  repr      =  as.numeric(coef(mod_reproduct))
  # recruit or not
  recr      =  as.numeric(coef(mod_recruit))
  # fruits production
  fruit      =  as.numeric(coef(mod_fruit))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4 - Define the IPM vital rate functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Growth function, 
# probability distribution function of size z1 next time
g_z1z <- function(z1, z){
  mean <- grow[1] + grow[2] * z
  sd <- grow_sd
  p_den_grow <- dnorm(z1, mean = mean, sd = sd)
  return(p_den_grow)
}

# Survival function, logistic regression
s_z <- function(z, log){
  linear_p <- surv[1] + surv[2] * z + surv[3] * log 
  p <- plogis(linear_p)                           
  return(p)
}

# Reproduction function, logistic regression
pb_z <- function(z, log,harv){
  if (z < 40) {
    p <- 0
  } else {
    linear_p <- repr[1] + 
      repr[2] * z + repr[3] * log + repr[4] * harv
    p <- plogis(linear_p)
  }
  return(p)
}

# Recruitment function, logistic regression
pr_z <- function(z1) {
  linear_p <- recr[1]+ recr[2] * z1
  p <- p <- plogis(linear_p)
  return(p)
}

# Fruit production function, given you are size z 
fruit_z <- function(z, log, harv){
  linear_p <- fruit[1] + fruit[2] * z + fruit[3] * log+ fruit[4] * harv
  p <- log(linear_p)
  return(p)  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 5 - Define the IPM component kernel functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define the survival-growth kernel
P_z1z <- function (z1, z, log) {
  return( s_z(z, log) * g_z1z(z1, z) )
}

# Define the reproduction kernel
F_z1z <- function (z, z1, log, harv) {
  return( s_z(z, log) * pr_z(z1) * pb_z(z, log,harv) *fruit_z(z, log,harv)) 
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 6 - Define functions to numerically implement the IPM
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

L<-0 #lower
U<-200 #upper

m<-200

h<-(U-L)/m

meshpts<-L+(1:m)*h-h/2

#Do a matrix with no values
Pnew<-matrix(NA,nrow=m,ncol=m)


#Fill the matrix with the values 
for(i in 1:m) {
  for(j in 1:m){
    Pnew[i,j]<-P_z1z(z1=meshpts[j], z=meshpts[i], log=10)+ F_z1z(z=meshpts[i], z1=meshpts[j],log=0,harv=0)
  }
}



colnames(Pnew)<-meshpts
rownames(Pnew)<-meshpts

# ----- Define a function for plotting a matrix ----- #
# ----- Define a function for plotting a matrix ----- #
# ----- Define a function for plotting a matrix ----- #
myimage <- function(x,sp,mainy=F,lam=1,colpal=rainbow, ...){
  MIN <- min(x,na.rm=T)
  MAX <- max(x,na.rm=T)
  yLabels <- round(as.numeric(rownames(x)),2)
  xLabels <- round(as.numeric(colnames(x)),2)
  title <-c()
  
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rev(colpal(256)[1:200])
  ColorLevels <- seq(MIN, MAX, length=length(ColorRamp))
  
  # Reverse Y axis
  #reverse <- nrow(x) : 1
  yLabels <- yLabels #[reverse]
  x <- x #[reverse,]
  
  # Data Map
  #t(x)
  par(mar = c(5.3,5.3,2.5,2))
  if(mainy==F){
    image(1:length(xLabels), 1:length(yLabels),x, col=ColorRamp,
          axes=FALSE, zlim=c(MIN,MAX),main=mainy,
          xlab="Size at time t",ylab="Size at time t+1")
  }
  else{
    image(1:length(xLabels), 1:length(yLabels),x, col=ColorRamp,
          axes=FALSE, zlim=c(MIN,MAX),main=mainy,
          xlab="Size at time t",ylab="Size at time t+1")
  }
  abline(0,1,lwd=3)
  
  box()
  
  
  if( !is.null(title) ){
    title(main=title)
  }
  xR<-xLabels[seq(from=1,to=length(xLabels),by=30)]
  yR<-yLabels[seq(from=1,to=length(xLabels),by=30)]
  
  #for some reason, this is going from 1:192. I guess its by pixels, not by 
  #rownames
  
  axis(BELOW<-1, at=seq(from=1,to=dim(x)[1],length=length(xR)), labels=xR, cex.axis=1.09,las=2)
  axis(LEFT <-2, at=seq(from=1,to=dim(x)[1],length=length(xR)), labels=yR, las= HORIZONTAL<-1,
       cex.axis=1.09)
  
  
  
  # Color Scale
  par(mar = c(3,2.5,2.5,2))
  image(1, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
  par(mfrow=c(1,1))
}



myimage(x=Pnew,sp="Brazil nut", mainy="Brazil nut")

eigen(Pnew)$values[1]

curve(dnorm(x,mean=11,sd=1),from=10,to=20)