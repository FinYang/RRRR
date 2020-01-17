data <- RRR_sim()
y <- data$y
x <- data$x
z <- data$z
res <- MM(y=data$y, x=data$x, z=data$z, itr = 100)
res <- SSUM(y=data$y, x=data$x, z=data$z)
res <- SAA(y=data$y, x=data$x, z=data$z, addon = 100)
res <- SAA_MM(y=data$y, x=data$x, z=data$z, addon = 100)
res <- RRR(y=data$y, x=data$x, z=data$z)
res <- SAA_GMLE(y=data$y, x=data$x, z=data$z)
data$spec$A
res$A

data$spec$B
res$B

cbind(data$spec$D, data$spec$A_0)
res$D

data$spec$Sigma
res$Sigma



data <- RRR_sim()
res <- RRR(y=data$y, x=data$x, z=data$z, mu = T)
res
data
data <- RRR_sim(D=NULL)
res <- RRR(y=data$y, x=data$x,  mu = T)
res
data
data <- RRR_sim(mu = NULL)
res <- RRR(y=data$y, x=data$x, z=data$z, mu = F)
res
data
data <- RRR_sim(D=NULL, mu=NULL)
res <- RRR(y=data$y, x=data$x,  mu = F)
res
data
