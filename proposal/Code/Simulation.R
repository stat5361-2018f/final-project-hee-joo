library(mvtnorm)
set.seed(20181215)
re = 100
cumulative.beta <- matrix(0,4,100)
cumulative.phi <- matrix(0,6,100)
for(i in 1:re)
{
##First blcok Data generation and Estimates ##
beta = matrix(c(0.1, 0.3, 0.2, 0.3),2,2)
phi = matrix(c(0.1, 0.1,0.2, 0.2,0.1,0.2),2,3)
n = 1000
DATA = DATA_GENERATION(n, beta, phi)
DATASETIING = DATAFORMAT(DATA$Dat, DATA$Y, n)
#initial K and Beta, And Phi
K = solve(var(Y))
Beta = c(0.1, 0.1, 0.1, 0.1)
Phi = c(0.1,0.1,0.1,0.1,0.1,0.1)
FRIST = NR(Beta, Phi, K, DATASETIING$Design, DATASETIING$Design2, DATASETIING$Y, maxit = 10000, tol = 1e-5 )
First.esti = Cumul(n, DATASETIING$Y,DATASETIING$Design, DATASETIING$Design2, FRIST$Beta, FRIST$Phi)

##Second blcok Data generation and Estimates##
beta = matrix(c(0.1, 0.3, 0.2, 0.3),2,2)
phi = matrix(c(0.1, 0.1,0.2, 0.2,0.1,0.2),2,3)
n = 1000
DATA = DATA_GENERATION(n, beta, phi)
DATASETIING = DATAFORMAT(DATA$Dat, DATA$Y, n)
#initial K and Beta, And Phi
K = solve(var(Y))
Beta = c(0.1, 0.1, 0.1, 0.1)
Gam = c(0.1,0.1,0.1,0.1,0.1,0.1)
Second = NR(Beta, Gam, K, DATASETIING$Design, DATASETIING$Design2, DATASETIING$Y, maxit = 10000, tol = 1e-5 )
Second.esti = Cumul(n, DATASETIING$Y,DATASETIING$Design, DATASETIING$Design2, Second$Beta, Second$Phi)
#cumulative estimates
cumulative.beta[,i] = solve(First.esti$A11 + Second.esti$A11)%*%(First.esti$A11%*%c(FRIST$Beta[1:2],FRIST$Beta[3:4])+
                  Second.esti$A11%*%c(Second$Beta[1:2],Second$Beta[3:4])) 
cumulative.phi[,i] = solve(First.esti$A22 + Second.esti$A22)%*%(First.esti$A22%*% c(FRIST$Phi[,1],FRIST$Phi[,2],FRIST$Phi[,3])+
                  Second.esti$A22%*%c(Second$Phi[,1],Second$Phi[,2],Second$Phi[,3])) 
}

apply(cumulative.beta,1, mean)
apply(cumulative.beta,1, sd)
apply(cumulative.phi,1, mean)
apply(cumulative.phi,1, sd)


