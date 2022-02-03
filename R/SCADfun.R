# Some functions about fSCAD penalty, refer to the code of Lin et al. (2017) with some modifications
slos_temp = function(Y, U, Maxiter, lambda, gamma, beta.basis, absTol, Cutoff)
{
  sum_m = length(Y)

  # beta b-spline basis
  rng = fda::getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  L = beta.basis$nbasis
  M = length(breaks) - 1
  norder = L-M+1
  d = L-M

  L2NNer = sqrt(M/(rng[2]-rng[1]))

  # roughness penalty matrix V
  V = fda::bsplinepen(beta.basis)
  VV = sum_m*gamma*V

  # calculate W
  W = slos.compute.weights(beta.basis)

  # initial estimate of b
  bHat = solve(t(U)%*%U+VV)%*%t(U)%*%Y

  bTilde = bHat

  if(lambda > 0)
  {
    changeThres = absTol
    bTilde = slosLQA(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a=3.7)
    bZero = (abs(bTilde) < Cutoff)
    bTilde[bZero] = 0

    bNonZero = !bZero

  }

  bTilde

}

#W_j
slos.compute.weights = function(basis)
{
  L       = basis$nbasis
  rng     = fda::getbasisrange(basis)
  breaks  = c(rng[1],basis$params,rng[2])
  M       = length(breaks) - 1
  norder  = L-M+1
  W       = array(0,dim=c(norder,norder,M))

  for (j in 1:M)
  {
    temp = fda::inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    W[,,j] = temp[j:(j+norder-1),j:(j+norder-1)]
  }
  W
}

# Step 2 in Lin et al. (2017)
slosLQA = function(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a)
{
  betaNormj = c(0,M)
  bZeroMat = rep(FALSE,L)
  betaNorm = Inf
  sum_m = length(Y)

  bZeroMat[c(1, L)] <- TRUE

  d <- L - M

  it = 1
  while(it <= Maxiter)
  {
    betaNormOld = betaNorm
    betaNorm = sqrt(sum(bHat^2))

    change = (betaNormOld-betaNorm)^2
    if(change < absTol) break

    lqaW = NULL
    lqaWk = matrix(0,L,L)
    for(j in 1:M)
    {
      index = j:(j+d)
      betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j]%*%bHat[j:(j+d)]
      cjk = Dpfunc(sqrt(betaNormj[j])*L2NNer,lambda,a)

      if(cjk != 0)
      {
        if(betaNormj[j] < absTol){
          bZeroMat[index] = TRUE
        }else{
          lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/sqrt(betaNormj[j]))*W[,,j]
        }
      }
    }

    lqaW = lqaWk
    lqaW = lqaW / 2 #W^{(0)}

    bZeroVec = bZeroMat
    bNonZeroVec = !bZeroVec

    UtU = t(U[,bNonZeroVec])%*%U[,bNonZeroVec]
    Ut = t(U[,bNonZeroVec])
    Vp = sum_m*gamma*V[bNonZeroVec,bNonZeroVec]

    theta = solve(UtU+Vp+sum_m*lqaW[bNonZeroVec,bNonZeroVec,drop=F],Ut %*% Y)
    bHat = matrix(0,length(bNonZeroVec),1)
    bHat[bNonZeroVec] = theta

    # print(it)

    it = it + 1
  }
  bHat
}

# SCAD function
Dpfunc = function(u,lambda,a)
{
  if(u<=lambda) Dpval = lambda
  else if(u<a*lambda) Dpval = -(u-a*lambda)/(a-1)
  else Dpval = 0
  Dpval
}
