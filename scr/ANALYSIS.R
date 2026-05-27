rm(list=ls())
# Load necessary libraries
library(rpart)
library(rpart.plot)
library(mgcv)

# DATA IMPORTATION
df_standardized <- read.csv("data/STANDARDIZED_DATA.csv", 
                            sep=',', dec=".", check.names=F, fileEncoding = "LATIN1")
# List of variables to analyze
data <- df_standardized[, c(3,4,7,11:13)]

## Outliers detection : 2nd methode
library(robustbase)

combinaisons <- list(c(1,3), c(2:3), c(1:2))
for(i in 1:3){
  X <- data[,combinaisons[[i]]]
  
  mcd <- covMcd(X)
  
  # distance de Mahalanobis robuste
  D2 <- mahalanobis(
    X,
    center = mcd$center,
    cov = mcd$cov
  )
  
  cutoff <- qchisq(
    0.99,
    df = ncol(X)
  )
  
  print(which(D2 > cutoff))
}
outliers <- c(11, 12, 36)

## J'enleve les outliers
data <- data[-outliers,]
data_1 = read.csv("data/RAW_DATA.csv", 
                  sep=',', dec=".", 
                  check.names=F, fileEncoding = "LATIN1")[-c(outliers),c(3,4,7,11,16,17)]

# data <- data
var1 = scale(data$`Total Distance`)
var2 = scale(data$`High Intensity Running (> 15 km/h)`)
var3 = scale(data$`Sprinting (> 25 km/h)`)
VAR <- cbind(var1, var2, var3)

new_data <- data.frame(data)
interact.1 <- scale(new_data$MFO*new_data$VO2max)
interact.2 <- scale(new_data$MFO*new_data$PmaxH)
interact.3 <- scale(new_data$VO2max*new_data$PmaxH)

data <- data.frame(MFO = scale(new_data$MFO), VO2max = scale(new_data$VO2max), PmaxH = scale(new_data$PmaxH),
                   interact.1, interact.2, interact.3)
n = dim(data)[1]

formula <- c(variable~MFO, variable~VO2max, variable~PmaxH, 
             variable~interact.1, variable~interact.2, variable~interact.3)
set.seed(1)
B = 1000
TD = HIR = SPRINT = NULL
TD.value = HIR.value = SPRINT.value = NULL
for (i in 1:6){
  td=hir=sprint=NULL
  td.value = hir.value = sprint.value =NULL
  for (cnt in 1:B){
    boostrap=sample(1:n, n, replace = T)
    # TD
    data$variable = VAR[,1]
    td = c(td,lm(formula[[i]], data=data[boostrap,])$coefficients[-1])
    # HIR
    data$variable = VAR[,2]
    hir = c(hir,lm(formula[[i]], data=data[boostrap,])$coefficients[-1])
    # SPRINT
    data$variable = VAR[,3]
    sprint = c(sprint,lm(formula[[i]], data=data[boostrap,])$coefficients[-1])
  }
  TD = cbind(TD,td)
  HIR = cbind(HIR, hir)
  SPRINT = cbind(SPRINT,sprint)
}

### TD
png("figures/figure1.png", width = 4000, height = 4000, res = 600, type = "cairo")
layout(matrix(c(1, 1, 2, 3), nrow=2, byrow = T))

par(mar=c(5,5,2,2))
plot(NA, xlim=c(0.5,(6+.5)), ylim=c(-1,1),
     axes=F, xlab="", ylab="95% confidence interval of\nthe correlation coefficients", main="Total distance (m)")
abline(h=0, col="black", lty='dashed')
axis(2)
bornes_inf = apply(TD,2,quantile,probs=.025)
bornes_sup = apply(TD,2,quantile,probs=.975)
medianes = apply(TD,2,quantile,probs=.5)
X.Y=data.frame(x1=1:6, y1=bornes_inf,
               x2=1:6, y2=bornes_sup)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]

segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
points(x = 1:6, y = medianes, lwd = 3, cex = 2 , pch = 18, col = 'grey')

axis(1, at=1:6, labels= rep("", 6))
text(x=1.1:6.1, y=-1.25, c("MFO", "VO2max", "PmaxH", "MFO:VO2max","MFO:PmaxH","VO2max:PmaxH"), 
     xpd=NA, srt=45, pos=2)
text(c(2,4), 0.9, c("**", "*"), cex=2)

apply(TD, 2, function(x) quantile(x, p=0.025))
par(mar=c(2,2,2,2))
summary(lm(var1 ~ interact.1))
summary(lm(var1 ~ data$VO2max))
model <- gam(Total.Distance ~ s(VO2max, by=MFO), data=new_data)
vis.gam(model, theta = 315, n.grid = 20, lwd = 0.2, color = "gray",
        zlab = "Total distance")

model.tree <- rpart(data_1$`Total Distance` ~ ., maxdepth = 1, data = data_1[,c(1:3)])
model.tree
par(mar=c(4,4,3,2))
boxplot(data_1$`Total Distance` ~ new_data$VO2max>=55.565, 
        xlab=expression(VO2*"max" ~ "(" * mL * "." * min^{-1} * "." * kg^{-1} * ")"), 
        ylab="Total Distance (m)", axes=F, col=c("gray"))
axis(2)
axis(1, at=1:2, labels=c("<55.5",">55.5"))
segments(x0=1,x1=2,y0=10500, xpd=NA, lwd=2)
text(1.5, 10700, xpd=NA, "***", cex=2)
t.test(new_data$Total.Distance ~ new_data$VO2max>=55.565)
dev.off()

### HIR
png("figures/figure2.png", width = 4000, height = 4000, res = 600, type = "cairo")
layout(matrix(c(1,1,
                2,3), nrow=2, byrow = T))
par(mar=c(5,5,2,2))
plot(NA, xlim=c(0.5,(6+.5)), ylim=c(-1,1),
     axes=F, xlab="", ylab="95% confidence interval of\nthe correlation coefficients", 
     main=expression("High-Intensity Running" ~ "(speed > 15" ~ km * "." * h^{-1} ~ " in meters)"))
abline(h=0, col="black", lty='dashed')
axis(2)
bornes_inf = apply(HIR,2,quantile,probs=.025)
bornes_sup = apply(HIR,2,quantile,probs=.975)
medianes = apply(HIR,2,quantile,probs=.5)
X.Y=data.frame(x1=1:6, y1=bornes_inf,
               x2=1:6, y2=bornes_sup)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]
segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
points(x=1:6,y=medianes, lwd=3, cex=2 , pch=18, col='gray')
# points(x=1,y=medianes[1], lwd=3, cex=2 , pch=15, col='orange')
axis(1, at=1:6, labels= rep("", 6))
text(x=1.1:6.1, y=-1.25, c("MFO", "VO2max", "PmaxH", "MFO:VO2max","MFO:PmaxH","VO2max:PmaxH"), 
     xpd=NA, srt=45, pos=2)
text(c(1, 4, 5), 0.9, c("***", "**", "**"), cex = 2)

par(mar=c(4,5,3,2))
plot(data_1$`High Intensity Running (> 15 km/h)` ~ MFO, data = data_1,
     ylab=expression("High-Intensity Running" ~ "(speed > 15" ~ km * "." * h^{-1} * " in m)"), 
     xlab=expression("MFO" ~ "(" * g * "." * min^{-1} * ")"))

abline(model<-lm(data_1$`High Intensity Running (> 15 km/h)` ~ MFO, data = data_1))
cor.test(data_1$`High Intensity Running (> 15 km/h)`, data_1$MFO)
text(0.6, 1200, expression(paste(r, ": 0.4, ", r^2, ": 0.16, p=0.006")), cex = 0.8, pos = 4)

par(mar=c(4,4,3,2))
boxplot(data_1$`High Intensity Running (> 15 km/h)` ~ data_1$MFO>=0.728, ylim = c(1200, 2600),
        xlab = expression("MFO" ~ "(" * g * "." * min^{-1} * ")"), 
        ylab = expression("High-Intensity Running (speed > 15" ~ km * "." * h^{-1} ~ "in m)"), 
        axes = F, col = c("gray"))
axis(2)
axis(1, at = 1:2, labels=c("<0.73",">0.73"))
segments(x0 = 1, x1 = 2, y0 = 2600, xpd = NA, lwd = 2)
text(1.5, 2700, xpd = NA, "**", cex = 2)
t.test(var2 ~ new_data$MFO>=0.728)
dev.off()

# SPRINTING
png("figures/figure3.png", width = 4000, height = 4000, res = 600, type = "cairo")
layout(matrix(c(1,1,1,1,
                2,2,3,4), nrow=2, byrow = T))
par(mar=c(6,6,2,2))
plot(NA, xlim=c(0.5,(6+.5)), ylim=c(-1,1),
     axes=F, xlab="", ylab="95% confidence interval of\nthe correlation coefficients", 
     main=expression("Sprinting (speed > 25" ~ km * "." * h^{-1} ~ "in meters)"))
abline(h=0, col="black", lty='dashed')
axis(2)
bornes_inf = apply(SPRINT,2,quantile,probs=.025)
bornes_sup = apply(SPRINT,2,quantile,probs=.975)
medianes = apply(SPRINT,2,quantile,probs=.5)
X.Y=data.frame(x1=1:6, y1=bornes_inf,
               x2=1:6, y2=bornes_sup)
X.Y=X.Y[which(is.na(X.Y[,2])==F),]
segments(X.Y$x1, X.Y$y1, X.Y$x2, X.Y$y2, lwd=2)
points(x=1:6,y=medianes, lwd=3, cex=2 , pch=18, col='grey')

axis(1, at=1:6, labels= rep("", 6))
text(x=1.1:6.1, y=-1.25, c("MFO", "VO2max", "PmaxH", "MFO:VO2max","MFO:PmaxH","VO2max:PmaxH"), 
     xpd=NA, srt=45, pos=2)
summary(lm(var3 ~ interact.2))
summary(lm(var3 ~ new_data$MFO))

par(mar=c(2,2,2,2))
model <- gam(Sprinting....25.km.h. ~ s(PmaxH, by=MFO, sp=2), data=new_data)
vis.gam(model, theta = 320, n.grid = 20, lwd = 0.2, color = "gray",
        zlab = "Sprinting")

par(mar=c(4,4,2,2))
model.tree <- rpart(data_1$`Sprinting (> 25 km/h)` ~ ., maxdepth = 1, data = data_1[, c(1:3)])
model.tree

boxplot(data_1$`Sprinting (> 25 km/h)` ~ data_1$MFO>=0.698, 
        xlab = expression("MFO" ~ "(" * g * "." * min^{-1} * ")"), 
        ylab = expression("Sprinting (speed > 25" ~ km * "." * h^{-1} ~ "in meters)"), axes=F, col=c("grey"))
axis(2)
axis(1, at=1:2, labels=c("<0.7", ">0.7"))
segments(x0=1,x1=2,y0=290, xpd=NA, lwd=2)
text(c(1.5), 300, xpd = NA, "*", cex = 1.5)
wilcox.test(data_1$`Sprinting (> 25 km/h)` ~ data_1$MFO>=0.698)

par(mar=c(4,4,2,2))
boxplot(data_1$`Sprinting (> 25 km/h)` ~ data_1$PmaxH>=18.88, 
        xlab = expression("PmaxH" ~ "(" * W * "." * kg^{-1} * ")"), 
        ylab = expression("Sprinting (speed > 25" ~ km * "." * h^{-1} ~ "in meters)"), axes=F, col=c("grey"))
axis(2)
axis(1, at = 1:2, labels = c("<19",">19"))
segments(x0 = 1, x1 = 2, y0 = 290, xpd = NA, lwd = 2)
text(c(1.5), 300, xpd=NA, "p=0.10", cex = 1)
wilcox.test(data_1$`Sprinting (> 25 km/h)` ~ data_1$PmaxH>=18.88)
dev.off()
