# setwd("C:/fat_lab")
setwd(paste(Sys.getenv("HOME"),"/Dropbox/PostDoc/2025_fatscaling/fatscaling/data",sep=""))

# NAVARRETE CAPTIVE v WILD SCALING COMPARISON

# data_orig <- read.csv("masskg_to_fatmasskg.csv")
data_n <- read.csv("Navarrete_fat.csv")
data_nkg <- data_n[,1:2]/1000
# colnames(data_orig) <- c("masskg", "fatkg")
colnames(data_nkg) <- c("masskg", "fatkg")

#IMPORTANT
##########################
#Wild versus captive nkg data:
wildpos <- which(data_n[,3]=="Wild")
captivepos <- which(data_n[,3]=="Captive")
wildmodel <- lm(log(data_nkg[wildpos,2]) ~ log(data_nkg[wildpos,1]))
summary(wildmodel)
captivemodel <- lm(log(data_nkg[captivepos,2]) ~ log(data_nkg[captivepos,1]))
summary(captivemodel)
plot(log(data_nkg[wildpos,1]),log(data_nkg[wildpos,2]),pch=16,col="black")
abline(wildmodel)
points(log(data_nkg[captivepos,1]),log(data_nkg[captivepos,2]),pch=16,col="blue")
abline(captivemodel,col="blue")

#ARE THE SLOPES DIFFERENT?
## vector of responses
y     <- c(log(data_nkg[wildpos,2]),log(data_nkg[captivepos,2]))          
# whatever your y-variable is
## predictor already logged
log_M <- c(log(data_nkg[wildpos , 1]),
           log(data_nkg[captivepos , 1]))
group <- factor(c( rep("wild",    length(wildpos)),
                   rep("captive", length(captivepos)) ))
dat <- data.frame(y, log_M, group)
# full model: different intercepts AND (critically) different slopes
mod_full <- lm(y ~ log_M * group, data = dat)
summary(mod_full)
##########################

