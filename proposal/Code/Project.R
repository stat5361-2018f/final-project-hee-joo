DATA_GENERATION = function(n, beta,  phi)
{   
   X = rnorm(n, 0, 1)
   Dat = cbind(1, X)
   Y = matrix(0,n,2)
   for(i in 1:n)
   {   
      Mu =  Dat %*% beta
      tau = Dat %*% phi[,2] 
      diag.phi1 = exp(Dat %*% phi[,1])
      diag.phi2 = exp(Dat %*% phi[,3])
      C = diag(c(diag.phi1[i],diag.phi2[i]))
      C[lower.tri(C)] = tau[i]
      Covari = solve( C %*% t(C))
      Y[i,] = rmvnorm(1,mean = Mu[i,], sigma = Covari)
   }
   return(list( Dat= Dat, Y=Y))
}   
DATAFORMAT = function(Dat,Y,n)
{   
   I = diag(1,2); I2 = diag(1,3)
   Design = array(0, dim = c(2,4,n))
   Design2 = array(0, dim = c(3,6,n))
   for(i in 1:n)
   {   
      Design[,,i] = kronecker(I,t(Dat[i,]))
      Design2[,,i] = kronecker(I2,t(Dat[i,]))
   }
   return(list( Design= Design, Design2=Design2, Y = Y))
}
NR = function(Beta, Gam, K, Design, Design2, Y, maxit = 10000, tol = 1e-5 )
{
   current.Beta = Beta
   current.Gam = Gam
   current.K = K
   iter <- U1 <- U2 <- I11 <- I22 <- 0
   Ubeta <- Ibeta <- Uphi <- Iphi <- 0 
   while(iter <= maxit)
   {
      for(i in 1:dim(Y)[1])
      {
         di1 = Design[2,3:4,i] %*% current.Gam[1:2]
         di2 = Design[2,3:4,i] %*% current.Gam[5:6]
         C =diag(c(exp(di1),exp(di2))); C[lower.tri(C)] = Design[2,3:4,i] %*% current.Gam[3:4] 
         current.K = C %*% t(C)
         R = Y[i,] - Design[,,i] %*% current.Beta
         Ubeta = t(Design[,,i]) %*% current.K %*% R + Ubeta
         Ibeta = t(Design[,,i]) %*% current.K %*% Design[,,i] + Ibeta
         S = R %*% t(R)
         G = c( (S[1,1]*C[1,1]^2 + S[2,1]*C[2,1]*C[1,1]), S[2,1]*C[1,1] + S[2,2]*C[2,1],C[2,2]^2*S[2,2])
         Uphi = -t(Design2[,,i]) %*% (G - c(1,0,1))/2  +Uphi
         H123 = c(2*S[1,1]*C[1,1]^2 + S[2,1]*C[2,1]*C[1,1], S[2,2], 2*S[2,2]*C[2,2]^2)
         H = diag(H123); H[1,2] <- H[2,1] <- S[2,1]*C[1,1]
         Iphi = t(Design2[,,i])%*% (H) %*% Design2[,,i] +Iphi
      }  
      U1 = Ubeta
      I11 = Ibeta
      U2 = Uphi
      I22 =  Iphi
      current.est = c(current.Beta, current.Gam)
      new.est = current.est + c(solve(I11, tol = 1e-20) %*% U1, solve(I22, tol = 1e-20) %*% U2)
      if( sum((new.est -current.est)^2) < tol) break
      current.Beta = new.est[1:length(Beta)]
      current.Gam = new.est[-c(1:length(Beta))]
      print(sum((new.est -current.est)^2))
   }
   current.Beta = matrix(new.est[1:length(Beta)],2,2)
   current.Gam = matrix(new.est[-c(1:length(Beta))],2,3)
   return(list(Beta=current.Beta, Phi = current.Gam))
}   
Cumul = function(n, Y,Design, Design2, current.beta, current.phi, maxit = 10000, tol = 1e-5)
{
   Ibeta<-Iphi<-0
   for(i in 1:n)
   {
      di1 = Design[2,3:4,i] %*% current.phi[1:2]
      di2 = Design[2,3:4,i] %*% current.phi[5:6]
      C =diag(c(exp(di1),exp(di2))); C[lower.tri(C)] = Design[2,3:4,i] %*% current.phi[3:4] 
      current.K = C %*% t(C)
      R = Y[i,] - Design[,,i] %*% c(current.beta[,1],current.beta[,2])
      Ibeta = t(Design[,,i]) %*% current.K %*% Design[,,i] + Ibeta
      S = R %*% t(R)
      H123 = c(2*S[1,1]*C[1,1]^2 + S[2,1]*C[2,1]*C[1,1], S[2,2], 2*S[2,2]*C[2,2]^2)
      H = diag(H123); H[1,2] <- H[2,1] <- S[2,1]*C[1,1]
      Iphi = t(Design2[,,i])%*% (H) %*% Design2[,,i] +Iphi
   }  
   I11 = Ibeta
   I22 =  Iphi
   return(list(current.beta = current.beta, current.phi= current.phi,A11 = I11, A22 = I22 ))
}

