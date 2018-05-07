library(unmarked)
# Simulate data using the multinomial-Poisson model with a
# repeated constant-interval removal design. for FOUR primary periods - 2 days, 2 nights
n <- 40 # number of sites
T <- 4 # number of primary periods
J <-  2# number of secondary periods
lam <- 3
phi <- 0.5
p <- 0.3

### 
#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam) # Local population size
N <- matrix(NA, n, T) # Individuals available for detection
for(i in 1:n) {
N[i,] <- rbinom(T, M[i], phi)
y[i,,1] <- rbinom(T, N[i,], p) # Observe some
Nleft1 <- N[i,] - y[i,,1] # Remove them
y[i,,2] <- rbinom(T, Nleft1, p) # ...
Nleft2 <- N[i,] - y[i,,2]
}
y.ijt1 <- cbind(y[,1,],y[,2,],y[,3,],y[,4,])

######### If you have 2 secondary periods, add these to dataset
### y[i,,2] <- rbinom(T, Nleft1, p) # ...
### Nleft2 <- N[i,] - y[i,,1] # Remove them
#########


# Simulate data using the multinomial-Poisson model with a
# repeated constant-interval removal design. for TWO primary periods - 2 nights
n <- 40 # number of sites
T <- 2 # number of primary periods
J <-  2# number of secondary periods
lam <- 3
phi <- 0.5
p <- 0.3

### 
#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam) # Local population size
N <- matrix(NA, n, T) # Individuals available for detection
for(i in 1:n) {
N[i,] <- rbinom(T, M[i], phi)
y[i,,1] <- rbinom(T, N[i,], p) # Observe some
Nleft1 <- N[i,] - y[i,,1] # Remove them
y[i,,2] <- rbinom(T, Nleft1, p) 
}
y.ijt2 <- cbind(y[,1,],y[,2,])

####################################
####################################
####################################
####################################
## Make a dataframe
umf2 <- unmarkedFrameGMM(y=y.ijt2, numPrimary=T, type="removal")
#?unmarkedFrameGMM
str(umf2)

### Run model with null covariates
(m2 <- gmultmix(~1, ~1, ~1, data=umf2, K=50))


### Get numbers
backTransform(m1, type="lambda") # Individuals per plot
backTransform(m1, type="phi") # Probability of being avilable
(p <- backTransform(m1, type="det")) # Probability of detection
p <- coef(p)
# Multinomial cell probabilities under removal design
c(p, (1-p) * p, (1-p)^2 * p)
# Or more generally:
head(getP(m1))
# Empirical Bayes estimates of super-population size
38 gpcount
re <- ranef(m1)
plot(re, layout=c(5,5), xlim=c(-1,20), subset=site%in%1:25)

############
Nhat <- function(fm) {
    N <- sum(predict(fm, type="state")$Predicted, na.rm=TRUE)
    }  

chisq <- function(fm) {
umf <- getData(fm)
y <- getY(umf)
y[y>1] <- 1
sr <- fm@sitesRemoved
if(length(sr)>0)
y <- y[-sr,,drop=FALSE]
fv <- fitted(fm, na.rm=TRUE)
y[is.na(fv)] <- NA
sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}


fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

##########################################
####  Check model GOF
##########################################
topModel<-m1
topModel<-m2

summary(topModel)
(pb.plot <- parboot(topModel, chisq, nsim=10, report=1)) ##shows distribution graphic
plot(pb.plot) ##shows distribution graphic
pb.plot
