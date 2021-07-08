##################################################################################################################
### Functions for Direct Filter Approach
##################################################################################################################
# Most of the functions are copied from the R-Package MDFA by Marc Wildi

########################################################
### Grunds?tzliche Funktionen des Direct Filter Approach
########################################################
### Alle Funktionen von Marc Wildi Package: MDFA

### Periodogram

# Plot
periodogram.Plot <- function(x,y = NULL,a = FALSE,b = FALSE){
  if(b){
   Name1 <- "Periodogram"  
  } else {
    Name1 <- "Filter"
  }
  
  if(a){
    len<-nrow(matrix(y))
    per1<-0:(len/2)
    DFT1<-per1
    
    # Discrete Fourier Transformation
    for (k in 0:(len/2)){
      cexp <- exp(-1.i*(1:len)*2*pi*k/len)
      DFT1[k+1]<-sum(cexp*y*sqrt(1/(2*pi*len)))
    }
    
    if (abs(as.integer(len/2)-len/2)<0.1)
      DFT1[k+1]<-DFT1[k+1]/sqrt(2)
    per1<-abs(DFT1)^2
  } 
  
  len1<-nrow(matrix(x))
  per<-0:(len1/2)
  DFT<-per
  
  # Discrete Fourier Transformation
  for (k in 0:(len1/2)){
    cexp <- exp(-1.i*(1:len1)*2*pi*k/len1)
    DFT[k+1]<-sum(cexp*x*sqrt(1/(2*pi*len1)))
  }
  
  if (abs(as.integer(len1/2)-len1/2)<0.1)
    DFT[k+1]<-DFT[k+1]/sqrt(2)
  per<-abs(DFT)^2
  
  # Plot
  p <- plot_ly() %>%
    layout(yaxis = list(
      title = "",
      showticklabels = FALSE,
      showgrid = FALSE),
      xaxis = list(
        ticktext = list("0","pi/6","2pi/6","3pi/6",
                        "4pi/6","5pi/6","pi"), 
        tickvals = as.list((1+0:6*len1/12)/length(per)*pi),
        tickmode = "array"
      )) %>% 
    add_lines(x = seq(1,length(per))/length(per)*pi,y = per, line = list(color = "black"), name = Name1) 
  
  if(a){
    p <- p %>%
      add_lines(x = seq(1,length(per1))/length(per1)*pi,y = per1, line = list(color = "red"),yaxis = "y2", name = "Underlying Data",
                opacity = 0.5) %>%
      layout(yaxis2 = list(
        overlaying = "y",
        side = "right",
        title = "",
        showticklabels = FALSE,
        showgrid = FALSE))}
  p
}

## Univariate
# Periodogram
per<-function(x,plot_T){
  
  # Initialisation
  len<-length(x)
  per<-0:(len/2)
  DFT<-per
  
  # Discrete Fourier Transformation
  for (k in 0:(len/2)){
    cexp <- exp(-1.i*(1:len)*2*pi*k/len)
    DFT[k+1]<-sum(cexp*x*sqrt(1/(2*pi*len)))
  }
  # Frequency zero receives weight 1/sqrt(2)
  #   The periodogram in frequency zero appears once only whereas all other frequencies are doubled
  
  # This is omitted now in order to comply with MDFA
  #   We now change the periodogram in the dfa estimation routines
  #  DFT[1]<-DFT[1]/sqrt(2)
  # Weighths wk: if length of data sample is even then DFT in frequency pi is scaled by 1/sqrt(2) (Periodogram in pi is weighted by 1/2)
  if (abs(as.integer(len/2)-len/2)<0.1)
    DFT[k+1]<-DFT[k+1]/sqrt(2)
  per<-abs(DFT)^2
  if (plot_T)
  {
    par(mfrow=c(1,1))
    plot(per,type="l",axes=F,xlab="Frequency",ylab="Periodogram",
         main="Periodogram")
    axis(1,at=1+0:6*len/12,labels=c("0","pi/6","2pi/6","3pi/6",
                                    "4pi/6","5pi/6","pi"))
    axis(2)
    box()
  }
  return(list(DFT=DFT,per=per))
}

## Multivariate
# Periodogram
periodogram_bp <- function(x, dd, n.pg) {
  # preparations
  n.fit  <- length(x)
  xx     <- as.vector(x[((n.fit - n.pg + 1) : n.fit)])
  npg2   <- n.pg / 2
  perall <- fourtrans <- 0 * 0 : npg2
  
  # case without a seasonal component
  if (dd < 3) {
    for (j in 0 : npg2) {
      fourtrans[j + 1] <- xx %*% exp((1 : (2* npg2)) * 1.i * j * pi / npg2) #length(xx)
      fourtrans[j + 1] <- fourtrans[j + 1] / sqrt(pi * n.pg)
      perall[j + 1] <- abs(fourtrans[j + 1])^2
    }
  }
  
  # case with a seasonal component, special treatment for pi/6
  if (dd >= 3) {
    for (j in (1 : npg2)[(-npg2 / 6) * (1 : 6)]) {
      fourtrans[j + 1] <- xx %*% exp((1 : (2 * npg2)) * 1.i * j * pi / npg2)
      term2 <- abs(1 - exp(j * 1.i * pi / npg2))^2
      term3 <- abs(1 - exp(12 * j * 1.i * pi / npg2))^2
      perall[j + 1] <- abs(fourtrans[j + 1]) / (term2 * term3)
    }
    perall[(npg2 / 6) * (1 : 6) + 1] <- max(perall) * 100000
  }
  # Weights wk: if length of data sample is even then DFT in frequency pi is scaled by 1/sqrt(2) (Periodogram in pi is weighted by 1/2)
  if (abs(as.integer(n.pg/2)-n.pg/2)<0.1)
  {
    fourtrans[npg2+1]<-fourtrans[npg2+1]/sqrt(2)
    perall[npg2+1]<-perall[npg2+1]/(2)
  }
  
  # output
  return(list(perall = perall, fourtrans = fourtrans))
}

# Weight Function
spec_comp <- function(insamp, x, d) {
  # non-stationarity
  if(d == 1) {
    weight_func <- periodogram_bp(diff(x[1 : insamp, 1]), 1,
                                  insamp - 1)$fourtrans
    # explaining variables
    if(length(weight_func) > 1) {
      for(j in 2 : ncol(x)) {
        # since the data is integrated one uses the pseudo-periodogram:
        # diff(data) and d <- 1
        per <- periodogram_bp(diff(x[1 : insamp, j]), 1, insamp - 1)$fourtrans
        weight_func <- cbind(weight_func, per)
      }
    }
  } else {
    weight_func <- periodogram_bp(x[1 : insamp, 1], 0, insamp)$fourtrans
    # explaining variables
    if(length(weight_func) > 1) {
      for(j in 2 : ncol(x)) {
        per <- periodogram_bp(x[1 : insamp, j], 0, insamp)$fourtrans
        weight_func <- cbind(weight_func, per)
      }
    }
  }
  colnames(weight_func) <- colnames(x)
  
  return(list(weight_func = weight_func))
}

### Univariate Direct Filter Approach

dfa_analytic <- function(L,lambda,periodogram,Lag,Gamma,eta,cutoff,i1,i2)
{
  # Frequency 0 appears once only whereas all other frequencies are doubles: therefore we halve periodogram in frequency zero
  periodogram[1]<-periodogram[1]/2
  # Impose meaningful parameter restrictions
  lambda<-abs(lambda)
  eta<-abs(eta)
  K<-length(periodogram)-1
  # Define the amplitude weighting function weight_h (W(omega_k))
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0){
    weight_h<-periodogram*(c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(eta)))
  } else{
    weight_h<-periodogram*rep(1,K+1)
  }
  # First order filter restriction: assigning a `large' weight to frequency zero
  if (i1)
    weight_h[1]<-max(1.e+10,weight_h[1])
  
  X<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)*sqrt(weight_h)
  X_y<-exp(-1.i*Lag*pi*(0:(K))/(K))*rep(1,K+1)
  if (i2)
  {
    # Second order restriction: time shift in frequency zero vanishes
    for (l in 2:(L-1))
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                    cos((L-1-Lag)*pi*(0:(K))/(K))+
                    sqrt(1+Gamma*lambda)*1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                                                sin((L-1-Lag)*pi*(0:(K))/(K))))*sqrt(weight_h))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*
                        cos((L-1-Lag)*pi*(0:(K))/(K))+
                        1.i*(sin((l-1-Lag)*pi*(0:(K))/(K))-((l-1)/(L-1))*sin((L-1-Lag)*pi*(0:(K))/(K)))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
    # MA-Filterweights
    b<-as.vector(solve(xtx,tol = 1e-25)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
    # the last weight is a function of the previous ones through the second order restriction
    b<-c(b,-sum(b*(0:(length(b)-1)))/(length(b)))
  } else
  {
    for (l in 2:L)
    {
      X<-cbind(X,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                    sqrt(1+Gamma*lambda)*1.i*sin((l-1-Lag)*pi*(0:(K))/(K)))*sqrt(weight_h))
      X_y<-cbind(X_y,(cos((l-1-Lag)*pi*(0:(K))/(K))+
                        1.i*sin((l-1-Lag)*pi*(0:(K))/(K))))
    }
    xtx<-t(Re(X))%*%Re(X)+t(Im(X))%*%Im(X)
    # MA-Filterweights
    b<-as.vector(solve(xtx,tol = 1e-25)%*%(t(Re(X_y))%*%(Gamma*weight_h)))
  }
  # Transferfunction
  trffkt<-1:(K+1)
  trffkt[1]<-sum(b)
  for (k in 1:(K))#k<-1
  {
    trffkt[k+1]<-(b%*%exp(1.i*k*(0:(length(b)-1))*pi/(K)))
  }
  return(list(b=b,trffkt=trffkt))
}

### Multivariate Direct Filter Approach

mdfa_analytic <- function(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                        lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                        b0_H0,c_eta,weight_structure,white_noise,
                        synchronicity,lag_mat,troikaner)
{
  
  # Enforce meaningful parameter values
  lambda<-abs(lambda)
  eta<-abs(eta)
  cutoff<-min(abs(cutoff),pi)
  
  # Resolution of discrete frequency-grid
  K<-nrow(weight_func)-1
  weight_target<-weight_func[,1]
  weight_func<-weight_func*exp(-1.i*Arg(weight_target))
  
  white_noise_synchronicity_obj<-white_noise_synchronicity(weight_func,white_noise,synchronicity)
  
  weight_func<-white_noise_synchronicity_obj$weight_func
  
  if (i2&i1&any(weight_constraint==0))
  {
    print(rep("!",100))
    print("Warning: i2<-T is not meaningful when i1<-T and weight_constraint=0 (bandpass)")
    print(rep("!",100))
  }
  if (!(length(b0_H0)==L*(dim(weight_func)[2]-1))&length(b0_H0)>0)
    print(paste("length of b0_H0 vector is ",length(b0_H0),": it should be ",L*(dim(weight_func)[2]-1)," instead",sep=""))
  # The function spect_mat_comp rotates all DFTs such that first column is real!#Lag<-0
  
  spec_mat<-spec_mat_comp(weight_func,L,Lag,c_eta,lag_mat)$spec_mat
  #spec_mat[,2]  spec_math-spec_mat
  structure_func_obj<-structure_func(weight_func,spec_mat)
  
  
  Gamma_structure<-structure_func_obj$Gamma_structure
  spec_mat_structure<-structure_func_obj$spec_mat_structure
  Gamma_structure_diff<-structure_func_obj$Gamma_structure_diff
  spec_mat_structure_diff<-structure_func_obj$spec_mat_structure_diff
  
  # weighting of amplitude function in stopband
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0)
  {
    # Replicate results in McElroy-Wildi (trilemma paper)
    if (lin_eta)
    {
      eta_vec<-c(rep(1,omega_Gamma),1+rep(eta,K-omega_Gamma+1))
    } else
    {
      if (c_eta)
      {
        eta_vec<-c(rep(1,omega_Gamma+1),(1:(K-omega_Gamma))*pi/K +1)^(eta/10)
      } else
      {
        eta_vec<-c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(eta/2))
      }
    }
    weight_h<-weight_func*eta_vec
  } else
  {
    eta_vec<-rep(1,K+1)
    weight_h<-weight_func* eta_vec
  }
  # Frequency zero receives half weight
  weight_h[1,]<-weight_h[1,]*ifelse(c_eta,1,1/sqrt(2))
  # DFT target variable
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error). Note that target of strutural
  # additional structural elements and of original MDFA are the same i.e. rotation must be the same and weight_target are identical too.
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))
  
  # If structure (forecasting) is imposed then we have to undo the customization weighting due to eta (no emphasis of stopband)
  if (sum(abs(weight_structure))>0)
  {
    weight_target_structure<-weight_target_structure_diff<-weight_target/eta_vec
  } else
  {
    weight_target_structure<-weight_target_structure_diff<-rep(0,K+1)
  }
  # DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])
  
  # The spectral matrix is inflated in stopband: effect of eta
  spec_mat<-t(t(spec_mat)*eta_vec) #dim(spec_mat)  as.matrix(Im(spec_mat[150,]))
  
  # Compute design matrix and regularization matrix
  mat_obj<-mat_func(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean,b0_H0,c_eta,lag_mat)
  
  des_mat<-mat_obj$des_mat
  reg_mat<-mat_obj$reg_mat
  reg_xtxy<-mat_obj$reg_xtxy
  w_eight<-mat_obj$w_eight
  
  # Solve estimation problem
  mat_x<-des_mat%*%spec_mat#spec_mat[,1]<-1
  X_new<-t(Re(mat_x))+sqrt(1+Gamma*lambda)*1.i*t(Im(mat_x))
  # xtx can be written either in Re(t(Conj(X_new))%*%X_new) or as below:
  xtx<-t(Re(X_new))%*%Re(X_new)+t(Im(X_new))%*%Im(X_new)
  # If abs(weight_structure)>0 then additional structure is instilled into MDFA-criterion: one-step ahead MSE
  if (sum(abs(weight_structure))>0)
  {
    X_new_structure<-t(abs(weight_structure[1])*des_mat%*%spec_mat_structure)
    xtx_structure<-t(Re(X_new_structure))%*%Re(X_new_structure)+t(Im(X_new_structure))%*%Im(X_new_structure)
    
    X_new_structure_diff<-t(abs(weight_structure[2])*des_mat%*%spec_mat_structure_diff)
    xtx_structure_diff<-t(Re(X_new_structure_diff))%*%Re(X_new_structure_diff)+t(Im(X_new_structure_diff))%*%Im(X_new_structure_diff)
    
    Gamma_structure_diff<-abs(weight_structure)[2]*Gamma_structure
    Gamma_structure<-abs(weight_structure)[1]*Gamma_structure
    
  } else
  {
    X_new_structure<-x_new_structure_diff<-0*des_mat%*%spec_mat_structure
    xtx_structure<-xtx_structure_diff<-0*t(Re(X_new_structure))%*%Re(X_new_structure)+t(Im(X_new_structure))%*%Im(X_new_structure)
  }
  # The filter restrictions (i1<-T and/or i2<-T) appear as constants on the right hand-side of the equation:
  xtxy<-t(Re(t(w_eight)%*%spec_mat)%*%Re(t(spec_mat)%*%t(des_mat))+
            Im(t(w_eight)%*%t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%Im(t(t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%t(des_mat)))
  # scaler makes scales of regularization and unconstrained optimization `similar'
  scaler<-mean(diag(xtx))
  
  
  
  
  if (sum(abs(weight_structure))>0)
  {
    X_inv<-solve(xtx+xtx_structure+xtx_structure_diff+scaler*reg_mat,tol = 1e-25)
    bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)+
                             ((t(Re(X_new_structure)*weight_target_structure))%*%Gamma_structure)+
                             ((t(Re(X_new_structure_diff)*weight_target_structure_diff))%*%Gamma_structure_diff)-
                             xtxy-scaler*reg_xtxy))
  } else
  {
    X_inv<-solve(xtx+scaler*reg_mat,tol = 1e-25)
    bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)-xtxy-scaler*reg_xtxy)) #weight_target[1]<-1
  }
  
  b<-matrix(nrow=L,ncol=length(weight_h_exp[1,]))
  # Reconstruct original parameters from the set of possibly contrained ones
  bhh<-t(des_mat)%*%bh
  for (k in 1:L)
  {
    b[k,]<-bhh[(k)+(0:(length(weight_h_exp[1,])-1))*L]
  }
  
  
  
  
  weight_cm<-matrix(w_eight,ncol=(length(weight_h_exp[1,])))
  # Add level and/or time-shift constraints (if i1<-F and i2<-F then this matrix is zero)
  b<-b+weight_cm
  
  
  # Transferfunction
  trffkt<-matrix(nrow=K+1,ncol=length(weight_h_exp[1,]))
  trffkth<-trffkt
  trffkt[1,]<-apply(b,2,sum)
  trffkth[1,]<-trffkt[1,]
  
  
  for (j in 1:length(weight_h_exp[1,]))
  {
    for (k in 0:(K))
    {
      trffkt[k+1,j]<-(b[,j]%*%exp(1.i*k*lag_mat[,j+1]*pi/(K)))
    }
  }
  
  
  trt<-apply(((trffkt)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
  # DFA criterion which accounts for customization but not for regularization term
  # MDFA-Legacy : new normalization for all error terms below
  rever<-sum(abs(Gamma*weight_target-Re(trt)-1.i*sqrt(1+lambda*Gamma)*Im(trt))^2)*pi/(K+1)
  # MS-filter error : DFA-criterion without effects by lambda or eta (one must divide spectrum by eta_vec)
  MS_error<-sum((abs(Gamma*weight_target-trt)/eta_vec)^2)*pi/(K+1)
  # Definition of Accuracy, time-shift and noise suppression terms
  Gamma_cp<-Gamma[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gamma[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-(trt/eta_vec)[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-(trt/eta_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-(weight_target/eta_vec)[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-(weight_target/eta_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  # define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2)*2*pi/(K+1)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp)*2*pi/(K+1)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2)*2*pi/(K+1)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)*2*pi/(K+1)
  
  # The following computations are time-consuming if K is `large'
  # By default they are skipped i.e. troikaner=F
  # Additional output: degrees of freedom, AIC
  if (troikaner)
  {
    # The following derivations of the DFA-criterion are equivalent
    # They are identical with rever (up to normalization by (2*(K+1)^2))
    trth<-((X_new)%*%(X_inv%*%t(Re(X_new))))%*%(weight_target*Gamma)
    # The projection matrix
    Proj_mat<-((X_new)%*%(X_inv%*%t(Re(X_new))))
    # The residual projection matrix
    res_mat<-diag(rep(1,dim(Proj_mat)[1]))-Proj_mat
    # DFA criterion: first possibility (all three variants are identical)
    sum(abs(res_mat%*%(weight_target*Gamma))^2)*pi/(K+1)
    # Residuals (DFT of filter error)
    resi<-res_mat%*%(weight_target*Gamma)
    # DFA criterion: second possibility
    t(Conj(resi))%*%resi*pi/(K+1)
    t((weight_target*Gamma))%*%(t(Conj(res_mat))%*%(res_mat))%*%(weight_target*Gamma)*pi/(K+1)
    freezed_degrees_freedom<-2*Re(sum(diag(t(Conj(res_mat))%*%(res_mat))))
    #  Re(t(Conj(res_mat))%*%(res_mat))-Re(res_mat)
    degrees_freedom<-2*K+1-freezed_degrees_freedom
    # The following expression equates to zero i.e. projection matrix is a projection for real part
    #   Re(t(Conj(res_mat))%*%(res_mat))-Re(res_mat)
    
    M<-((X_new)%*%(X_inv%*%t(Conj(X_new))))
    res_M<-diag(rep(1,dim(M)[1]))-M
    # The following is a real number but due to numerical rounding errors we prefer to take the real part of the result
    edof<-Re(K+1-sum(eigen(res_M)$values))
    # Troikaner
    tic<-ifelse(freezed_degrees_freedom<2*K+1&freezed_degrees_freedom>1,log(rever)+2*(K-freezed_degrees_freedom+1)/(freezed_degrees_freedom-2),NA)
    tic<-log(rever)+2*edof/(2*K)
    return(list(b=b,trffkt=trffkt,rever=rever,freezed_degrees_freedom=freezed_degrees_freedom,tic=tic,degrees_freedom=degrees_freedom,Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error,edof=edof))
  } else
  {
    # Simplified (shorter) return
    return(list(b=b,trffkt=trffkt,rever=rever,Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error))
  }
}

# Functions for MDFA  
white_noise_synchronicity<-function(weight_func,white_noise,synchronicity)
{
  
  return(list(weight_func=weight_func))
}  

spec_mat_comp<-function(weight_func,L,Lag,c_eta,lag_mat)
{
  K<-length(weight_func[,1])-1
  weight_h<-weight_func
  # Frequency zero receives half weight
  # Chris mod: avoid halving
  # MDFA-Legacy: the DFT in frequency zero is weighted by 1/sqrt(2) (the weight of frequency zero is 0.5 in mean-square)
  if (!c_eta)
    weight_h[1,]<-weight_h[1,]/sqrt(2)
  # Extract DFT target variable (first column)
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))#Im(exp(-1.i*Arg(weight_target)))
  weight_target<-weight_target*exp(-1.i*Arg(weight_target))
  # DFT's explaining variables (target variable can be an explaining variable too)
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])
  spec_mat<-as.vector(t(as.matrix(weight_h_exp[1,])%*%t(as.matrix(rep(1,L)))))
  
  for (j in 1:(K))#j<-1  h<-2  lag_mat<-matrix(rep(0:1,3),ncol=3)   Lag<-2
  {
    omegak<-j*pi/K
    # We feed the daily structure of the mixed_frequency data lag_mat
    exp_mat<-matrix(nrow=L,ncol=ncol(weight_h_exp))
    for (h in 1:ncol(exp_mat))
      exp_mat[,h]<-exp(1.i*omegak*(lag_mat[,h+1]-Lag))
    # Inclusion of the time-varying lag-structure of the mixed-frequency data
    spec_mat<-cbind(spec_mat,as.vector(t(weight_h_exp[j+1,]*t(exp_mat))))
  }
  dim(spec_mat)
  return(list(spec_mat=spec_mat))#as.matrix(Re(spec_mat[1,]))
}

mat_func<-function(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean,b0_H0,c_eta,lag_mat)
{
  
  
  
  # MDFA-Legacy: new functions
  # Regularization
  reg_mat_obj<-reg_mat_func(weight_h_exp,L,c_eta,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth,lag_mat)
  
  Q_smooth<-reg_mat_obj$Q_smooth
  Q_cross<-reg_mat_obj$Q_cross
  Q_decay<-reg_mat_obj$Q_decay
  lambda_decay<-reg_mat_obj$lambda_decay
  lambda_smooth<-reg_mat_obj$lambda_smooth
  lambda_cross<-reg_mat_obj$lambda_cross
  
  # MDFA-Legacy: new functions
  # weight vector for constraints
  w_eight<-w_eight_func(i1,i2,Lag,weight_constraint,shift_constraint,L,weight_h_exp)$w_eight
  
  # MDFA-Legacy: new functions
  # Grand-mean and Constraints
  
  des_mat<-des_mat_func(i2,i1,L,weight_h_exp,weight_constraint,shift_constraint,Lag)$des_mat
  
  # Transforming back from grand-mean if necessary (des_mat assumes grand-mean)
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
    Q_centraldev_original<-centraldev_original_func(L,weight_h_exp)$Q_centraldev_original
  } else
  {
    Q_centraldev_original<-NULL
  }
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
    # apply A^{-1} to A*R giving R where A and R are defined in MDFA_Legacy book
    #and des_mat=A*R i.e.  t(Q_centraldev_original%*%t(des_mat)) is R (called des_mat in my code)
    des_mat<-t(Q_centraldev_original%*%t(des_mat))
  }
  
  # Here we fold all three regularizations (cross, smooth and decay) into a single reg-matrix
  if ((length(weight_h_exp[1,])>1))
  {
    if (grand_mean)
    {
      reg_t<-(Q_smooth+Q_decay+t(Q_centraldev_original)%*%Q_cross%*%Q_centraldev_original)
    } else
    {
      reg_t<-(Q_smooth+Q_decay+Q_cross)
    }
  } else
  {
    reg_t<-(Q_smooth+Q_decay)
  }
  
  
  
  # Normalize regularization terms (which are multiplied by des_mat) in order to disentangle i1/i2 effects
  reg_mat<-(des_mat)%*%reg_t%*%t(des_mat)#sum(diag(reg_mat))
  if (is.null(b0_H0))
    b0_H0<-rep(0,length(w_eight))
  b0_H0<-as.vector(b0_H0)
  if (lambda_smooth+lambda_decay[2]+lambda_cross>0)
  {
    disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))
    if (abs(sum(diag(reg_mat)))>0)
    {
      disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))
    } else
    {
      disentangle_des_mat_effect<-1
    }
    
    reg_mat<-reg_mat*disentangle_des_mat_effect
    reg_xtxy<-des_mat%*%reg_t%*%(w_eight-b0_H0)*(disentangle_des_mat_effect)#+t(w_eight)%*%reg_t%*%t(des_mat)
  } else
  {
    reg_xtxy<-des_mat%*%reg_t%*%(w_eight-b0_H0)#+t(w_eight)%*%reg_t%*%t(des_mat)
    dim(des_mat)
    dim(reg_t)
  }
  return(list(des_mat=des_mat,reg_mat=reg_mat,reg_xtxy=reg_xtxy,w_eight=w_eight))
}

structure_func<-function(weight_func,spec_mat)
{
  K<-nrow(weight_func)-1
  
  # We have to initialize spect_mat_structure and spec_mat_structure_diff by an arbitrary matrix of right dimensions
  spec_mat_structure<-spec_mat_structure_diff<-spec_mat
  Gamma_structure<-Gamma_structure_diff<-rep(0,K+1)
  
  return(list(Gamma_structure=Gamma_structure,spec_mat_structure=spec_mat_structure,
              Gamma_structure_diff=Gamma_structure_diff,spec_mat_structure_diff=spec_mat_structure_diff        ))
}

mat_func<-function(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean,b0_H0,c_eta,lag_mat)
{
  
  
  
  # MDFA-Legacy: new functions
  # Regularization
  reg_mat_obj<-reg_mat_func(weight_h_exp,L,c_eta,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth,lag_mat)
  
  Q_smooth<-reg_mat_obj$Q_smooth
  Q_cross<-reg_mat_obj$Q_cross
  Q_decay<-reg_mat_obj$Q_decay
  lambda_decay<-reg_mat_obj$lambda_decay
  lambda_smooth<-reg_mat_obj$lambda_smooth
  lambda_cross<-reg_mat_obj$lambda_cross
  
  # MDFA-Legacy: new functions
  # weight vector for constraints
  w_eight<-w_eight_func(i1,i2,Lag,weight_constraint,shift_constraint,L,weight_h_exp)$w_eight
  
  # MDFA-Legacy: new functions
  # Grand-mean and Constraints
  
  des_mat<-des_mat_func(i2,i1,L,weight_h_exp,weight_constraint,shift_constraint,Lag)$des_mat
  
  # Transforming back from grand-mean if necessary (des_mat assumes grand-mean)
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
    Q_centraldev_original<-centraldev_original_func(L,weight_h_exp)$Q_centraldev_original
  } else
  {
    Q_centraldev_original<-NULL
  }
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
    # apply A^{-1} to A*R giving R where A and R are defined in MDFA_Legacy book
    #and des_mat=A*R i.e.  t(Q_centraldev_original%*%t(des_mat)) is R (called des_mat in my code)
    des_mat<-t(Q_centraldev_original%*%t(des_mat))
  }
  
  # Here we fold all three regularizations (cross, smooth and decay) into a single reg-matrix
  if ((length(weight_h_exp[1,])>1))
  {
    if (grand_mean)
    {
      reg_t<-(Q_smooth+Q_decay+t(Q_centraldev_original)%*%Q_cross%*%Q_centraldev_original)
    } else
    {
      reg_t<-(Q_smooth+Q_decay+Q_cross)
    }
  } else
  {
    reg_t<-(Q_smooth+Q_decay)
  }
  
  
  
  # Normalize regularization terms (which are multiplied by des_mat) in order to disentangle i1/i2 effects
  reg_mat<-(des_mat)%*%reg_t%*%t(des_mat)#sum(diag(reg_mat))
  if (is.null(b0_H0))
    b0_H0<-rep(0,length(w_eight))
  b0_H0<-as.vector(b0_H0)
  if (lambda_smooth+lambda_decay[2]+lambda_cross>0)
  {
    disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))
    if (abs(sum(diag(reg_mat)))>0)
    {
      disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))
    } else
    {
      disentangle_des_mat_effect<-1
    }
    
    reg_mat<-reg_mat*disentangle_des_mat_effect
    reg_xtxy<-des_mat%*%reg_t%*%(w_eight-b0_H0)*(disentangle_des_mat_effect)#+t(w_eight)%*%reg_t%*%t(des_mat)
  } else
  {
    reg_xtxy<-des_mat%*%reg_t%*%(w_eight-b0_H0)#+t(w_eight)%*%reg_t%*%t(des_mat)
    dim(des_mat)
    dim(reg_t)
  }
  return(list(des_mat=des_mat,reg_mat=reg_mat,reg_xtxy=reg_xtxy,w_eight=w_eight))
}

reg_mat_func<-function(weight_h_exp,L,c_eta,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth,lag_mat)
{
  lambda_smooth<-100*tan(min(abs(lambda_smooth),0.999999999)*pi/2)
  lambda_cross<-100*tan(min(abs(lambda_cross),0.9999999999)*pi/2)
  if (c_eta)
  {
    # Chris replication: imposing non-linear transform to lambda_decay[1] too
    lambda_decay<-c(min(abs(tan(min(abs(lambda_decay[1]),0.999999999)*pi/2)),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  } else
  {
    lambda_decay<-c(min(abs(lambda_decay[1]),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  }
  # New regularization function
  
  Q_obj<-Q_reg_func(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat,grand_mean,c_eta)
  
  Q_smooth<-Q_obj$Q_smooth
  Q_decay<-Q_obj$Q_decay
  Q_cross<-Q_obj$Q_cross
  # The Q_smooth and Q_decay matrices address regularizations for original unconstrained parameters (Therefore dimension L^2)
  # At the end, the matrix des_mat is used to map these regularizations to central-deviance parameters
  # accounting for first order constraints!
  
  return(list(Q_smooth=Q_smooth,Q_cross=Q_cross,Q_decay=Q_decay,lambda_decay=lambda_decay,
              lambda_smooth=lambda_smooth,lambda_cross=lambda_cross))
}


Q_reg_func<-function(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat,grand_mean,c_eta)
{
  Q_smooth<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  Q_decay<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  # Cross-sectional regularization if dimension>1
  if ((length(weight_h_exp[1,])>1))
  {
    # The cross-sectional regularization is conveniently implemented on central-deviance parameters. The regularization is expressed on the
    # unconstrained central-deviance parameters (dimension L), then mapped to the original (unconstrained) parameters (dimension L) with Q_centraldev_original
    # and then maped back to central-deviance with constraint (dim L-1) with des_mat (mathematically unnecessarily complicate but more convenient to implement in code).
    Q_cross<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
    Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  } else
  {
    # 16.08.2012
    Q_cross<-NULL
  }
  for (i in 1:L)
  {
    # For symmetric filters or any historical filter with Lag>0 the decay must be symmetric about b_max(0,Lag)
    # lambda_decay is a 2-dim vector: the first component controls for the exponential decay and the second accounts for the strength of the regularization
    # The maximal weight is limited to 1e+4
    Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(i-1-max(0,Lag))),1e+4)
    # Original idea: with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale: this is not a clever idea...
    #   Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+1]-max(0,Lag))),1e+4)
    
    if (L>4)
    {
      if(i==1)
      {
        Q_smooth[i,i:(i+2)]<-c(1,-2,1)
      } else
      {
        if(i==2)
        {
          Q_smooth[i,(i-1):(i+2)]<-c(-2,5,-4,1)
        } else
        {
          if(i==L)
          {
            Q_smooth[i,(i-2):i]<-c(1,-2,1)
          } else
          {
            if(i==L-1)
            {
              Q_smooth[i,(i-2):(i+1)]<-c(1,-4,5,-2)
            } else
            {
              Q_smooth[i,(i-2):(i+2)]<-c(1,-4,6,-4,1)
            }
          }
        }
      }
    } else
    {
      print("L<=4: no smoothness regularization imposed!!!!!!!!!")
    }
  }
  
  if (length(weight_h_exp[1,])>1)
  {
    
    for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
    {
      if (L>4)
        Q_smooth[j*L+1:L,j*L+1:L]<-Q_smooth[1:L,1:L]
      for (i in 1:L)
        Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(i-1-max(0,Lag))),1e+4)
      # Original idea: with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale: this is not a clever idea...
      #      Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+j+1]-max(0,Lag))),1e+4)
      
    }
    Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
    if (L>1)
    {
      diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
      for (i in 2:length(weight_h_exp[1,]))   #i<-2
      {
        diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
        diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
        diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
      }
    }
    Q_centraldev_original<-solve(Q_centraldev_original,tol = 1e-25)
    
    # 06.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
    #   If grand_mean==T then the code replicates I-MDFA as released prior 30.07.2012
    #   If grand_mean==F then the new parametrization is used.
    #   Differences between both approaches: see section 7.2 of my elements paper posted on SEFBlog (both codes are identical when no regularization is imposed. Otherwise the later version (grand_mean==F) is logically more consistent becuase it treats all series identically (no asymmetry)).
    if (grand_mean)
    {
      diag(Q_cross[L+1:((length(weight_h_exp[1,])-1)*L),L+1:((length(weight_h_exp[1,])-1)*L)])<-
        rep(1,((length(weight_h_exp[1,])-1)*L))
    } else
    {
      #30.07.2012:new definition (parametrization) of Q_cross (Lambda_{cross} in the elements-paper)
      diag(Q_cross)<-1
      for (i in 1:length(weight_h_exp[1,]))
      {
        for (j in 1:L)
        {
          Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]<-Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]-1/length(weight_h_exp[1,])
        }
      }
    }
  } else
  {
    # define matrix for univariate case
    Q_centraldev_original<-NULL
  }
  # Normalizing the troika are new: disentangle the effect by L
  Q_decay<-Q_decay*lambda_decay[2]
  Q_cross<-Q_cross*lambda_cross                   #Qh<-Q_cross
  Q_smooth<-Q_smooth*lambda_smooth
  if (lambda_decay[2]>0)
  {
    # The second parameter in lambda_decay accounts for the strength of the regularization
    Q_decay<-lambda_decay[2]*(Q_decay/(sum(diag(Q_decay))))
  }
  if (lambda_cross>0)
  {
    if (c_eta)
    {
      Q_cross<-lambda_cross^2*(Q_cross/(sum(diag(Q_cross))))
    } else
    {
      Q_cross<-lambda_cross*(Q_cross/(sum(diag(Q_cross))))
    }
  }
  if (lambda_smooth>0&L>4)
  {
    Q_smooth<-lambda_smooth*(Q_smooth/(sum(diag(Q_smooth))))
  }
  return(list(Q_smooth=Q_smooth,Q_decay=Q_decay,Q_cross=Q_cross))
}

w_eight_func<-function(i1,i2,Lag,weight_constraint,shift_constraint,L,weight_h_exp)
{
  if (i1)
  {
    if (i2)
    {
      #               impose constraints to b_Lag, b_{Lag+1} instead of b_{L-1} and b_L
      #               Therefore the decay regularization does not potentially conflict with filter constraints
      if (Lag<1)
      {
        # corrected time-shift expression if A(0) different from 1
        w_eight<-c(-(Lag-1)*weight_constraint[1]-shift_constraint[1]*weight_constraint[1],
                   Lag*weight_constraint[1]+shift_constraint[1]*weight_constraint[1],rep(0,L-2))
      } else
      {
        # corrected time-shift expression if A(0) different from 1
        w_eight<-c(rep(0,Lag),weight_constraint[1]-shift_constraint[1]*weight_constraint[1],
                   shift_constraint[1]*weight_constraint[1],rep(0,L-Lag-2))
      }
      
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          
          if (Lag<1)
          {
            # corrected time-shift expression if A(0) different from 1
            w_eight<-c(w_eight,-(Lag-1)*weight_constraint[j]-shift_constraint[j]*weight_constraint[j],
                       Lag*weight_constraint[j]+shift_constraint[j]*weight_constraint[j],rep(0,L-2))
          } else
          {
            # corrected time-shift expression if A(0) different from 1
            w_eight<-c(w_eight,c(rep(0,Lag),weight_constraint[j]-shift_constraint[j]*weight_constraint[j],
                                 shift_constraint[j]*weight_constraint[j],rep(0,L-Lag-2)))
          }
        }
      }
    } else
    {
      if (Lag<1)
      {
        w_eight<-c(weight_constraint[1],rep(0,L-1))
      } else
      {
        w_eight<-c(rep(0,Lag),weight_constraint[1],rep(0,L-Lag-1))
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          if (Lag<1)
          {
            w_eight<-c(w_eight,weight_constraint[j],rep(0,L-1))
          } else
          {
            w_eight<-c(w_eight,rep(0,Lag),weight_constraint[j],rep(0,L-Lag-1))
          }
        }
      }
    }
  } else
  {
    # MDFA-Legacy : the case i2==T i1==F is fixed, at last.
    if (i2)
    {
      if (Lag<1)
      {
        w_eight<-rep(0,L)
      } else
      {
        w_eight<-rep(0,L)
      }
      
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          if (Lag<1)
          {
            w_eight<-c(w_eight,rep(0,L))
          } else
          {
            w_eight<-c(w_eight,rep(0,L))
          }
        }
      }
    } else
    {
      w_eight<-rep(0,L*length(weight_h_exp[1,]))
    }
  }
  return(list(w_eight= w_eight))
  
}

des_mat_func<-function(i2,i1,L,weight_h_exp,weight_constraint,shift_constraint,Lag)
{
  # Here we implement the matrix which links freely determined central-deviance parameters and constrained original parameters
  # In the MDFA_Legacy book t(des_mat) corresponds to A%*%R (R links constrained and unconstrained parameters and A maps central-deviance to original parameters)
  #   Please note that:
  #   1. Here I'm working with central-deviance parameters
  #   2. The same matrix R applies to either parameter set
  #   3. If I work with central-deviance parameters then R maps the freely determined set to the constrained (central-deviance)
  #       and A then maps the constrained (central-deviance) set to original constrained parameters.
  
  if (i2)
  {
    if (i1)
    {
      # First and second order restrictions
      des_mat<-matrix(data=rep(0,(L-2)*L*(length(weight_h_exp[1,]))^2),nrow=(L-2)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      
      for (i in 1:(L-2))
      {
        if (Lag<1)
        {
          des_mat[i,i+2+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<-i
          des_mat[i,2+(0:(length(weight_h_exp[1,])-1))*L]<--(i+1)
          
        } else
        {
          des_mat[i,ifelse(i<Lag+1,i,i+2)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
          des_mat[i,Lag+2+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
        }
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
        {
          for (i in 1:(L-2))
          {
            
            if (Lag<1)
            {
              des_mat[i+j*(L-2),i+2]<--1
              des_mat[i+j*(L-2),1]<--i
              des_mat[i+j*(L-2),2]<-(i+1)
            } else
            {
              des_mat[i+j*(L-2),ifelse(i<Lag+1,i,i+2)]<--1
              des_mat[i+j*(L-2),Lag+1]<--ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
              des_mat[i+j*(L-2),Lag+2]<--ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
            }
            
            if (Lag<1)
            {
              des_mat[i+j*(L-2),i+2+j*L]<-1
              des_mat[i+j*(L-2),1+j*L]<-i
              des_mat[i+j*(L-2),2+j*L]<--(i+1)
            } else
            {
              des_mat[i+j*(L-2),ifelse(i<Lag+1,i+j*L,i+2+j*L)]<-1
              des_mat[i+j*(L-2),Lag+1+j*L]<-ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
              des_mat[i+j*(L-2),Lag+2+j*L]<-ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
            }
          }
        }
      }
    } else
    {
      
      des_mat<-matrix(data=rep(0,(L-1)*L*(length(weight_h_exp[1,]))^2),nrow=(L-1)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      # MDFA-Legacy : the case i2==T i1==F is fixed, at last
      # Avoid singularities in the time-shift specification
      epsi<-1.e-02
      if (Lag>=1)
        shift_constraint[which(abs(shift_constraint)<epsi)]<-shift_constraint[which(abs(shift_constraint)<epsi)]+epsi
      if (Lag<1)
        shift_constraint[abs(-Lag-shift_constraint)<epsi]<-shift_constraint[which(abs(-Lag-shift_constraint)<epsi)]+epsi
      for (i in 1:(L-1))
      {
        if (Lag<1)
        {
          # Forecast and nowcast
          # we have to differentiate the case with vanishing shift and non-vansihing shift
          #    For a vanishing shift we select b2 for the constraint because then we are sure that the equations can be solved
          #    For a non-vanishing shift we select b1 because then we are sure that the equations can be solved (irrespective of the shift)
          # See MDFA-legacy book for reference
          
          # shift is different from zero: b1 is isolated (new formulas)
          des_mat[i,ifelse(i<1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<--(-Lag+i-shift_constraint[1])/
            (-Lag-shift_constraint[1])
        } else
        {
          # Backcast: we have to differentiate the case with vanishing shift and non-vansihing shift
          #    For a vanishing shift we select b2 for the constraint because then we are sure that the equations can be solved
          #    For a non-vanishing shift we select b1 because then we are sure that the equations can be solved (irrespective of the shift)
          # See MDFA-legacy book for reference
          # shift is different from zero: b1 is isolated (new formulas: the signs change because of multiplication with -s)
          des_mat[i,ifelse(i<Lag+1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<-
            ifelse(i<Lag+1,(-(Lag+1-i)-shift_constraint[1])/shift_constraint[1],(-(Lag-i)-shift_constraint[1])/shift_constraint[1])
          
        }
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))#j<-1
        {
          for (i in 1:(L-1))
          {
            if (Lag<1)
            {
              # shift is different from zero: b1 is isolated (new formulas)
              des_mat[i+j*(L-1),ifelse(i<1,i,i+1)]<--1
              des_mat[i+j*(L-1),1]<-(-Lag+i-shift_constraint[j+1])/(-Lag-shift_constraint[j+1])
              
            } else
            {
              # shift is different from zero: b1 is isolated (new formulas: the signs change because of multiplication with -s)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)]<--1
              des_mat[i+j*(L-1),Lag+1]<-
                -ifelse(i<Lag+1,(-(Lag+1-i)-shift_constraint[j+1])/shift_constraint[j+1],(-(Lag-i)-shift_constraint[j+1])/shift_constraint[j+1])
              
            }
            
            if (Lag<1)
            {
              # shift is different from zero: b1 is isolated (new formulas)
              des_mat[i+j*(L-1),ifelse(i<1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),1+j*L]<--(-Lag+i-shift_constraint[j+1])/(-Lag-shift_constraint[j+1])
              
            } else
            {
              # shift is different from zero: b1 is isolated (new formulas: the signs change because of multiplication with -s)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),Lag+1+j*L]<-
                ifelse(i<Lag+1,(-(Lag+1-i)-shift_constraint[j+1])/shift_constraint[j+1],(-(Lag-i)-shift_constraint[j+1])/shift_constraint[j+1])
              
            }
            
          }
        }
      }
      
    }
  } else
  {
    if (i1)
    {
      des_mat<-matrix(data=rep(0,(L-1)*L*(length(weight_h_exp[1,]))^2),nrow=(L-1)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      for (i in 1:(L-1))
      {
        # The i1-constraint is imposed on b_max(0,Lag) (instead of b_L ) in order to avoid a conflict with the exponential decay requirement
        if (Lag<1)
        {
          des_mat[i,i+1+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<--1
        } else
        {
          # Lag cannot be larger than (L-1)/2 (symmetric filter)
          des_mat[i,ifelse(i<Lag+1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<--1
        }
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
        {
          for (i in 1:(L-1))
          {
            # The i1-constraint is imposed on b_max(0,Lag) (instead of b_L ) in order to avoid a conflict with the exponential decay requirement
            if (Lag<1)
            {
              des_mat[i+j*(L-1),i+1]<--1
              des_mat[i+j*(L-1),1]<-1
              des_mat[i+j*(L-1),i+1+j*L]<-1
              des_mat[i+j*(L-1),1+j*L]<--1
            } else
            {
              # Lag cannot be larger than (L-1)/2 (symmetric filter)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)]<--1
              des_mat[i+j*(L-1),Lag+1]<-1
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),Lag+1+j*L]<--1
            }
          }
        }
      }
    } else
    {
      des_mat<-matrix(data=rep(0,(L)*L*(length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      for (i in 1:(L))
      {
        des_mat[i,i+(0:(length(weight_h_exp[1,])-1))*L]<-1
      }
      
      if (length(weight_h_exp[1,])>1)
      {
        
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
        {
          for (i in 1:(L))
          {
            des_mat[i+(j)*(L),i]<--1
            des_mat[i+(j)*(L),i+j*L]<-1
          }
        }
      }
    }
  }
  return(list(des_mat=des_mat))
  
}


Q_reg_func<-function(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat,grand_mean,c_eta)
{
  Q_smooth<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  Q_decay<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  # Cross-sectional regularization if dimension>1
  if ((length(weight_h_exp[1,])>1))
  {
    # The cross-sectional regularization is conveniently implemented on central-deviance parameters. The regularization is expressed on the
    # unconstrained central-deviance parameters (dimension L), then mapped to the original (unconstrained) parameters (dimension L) with Q_centraldev_original
    # and then maped back to central-deviance with constraint (dim L-1) with des_mat (mathematically unnecessarily complicate but more convenient to implement in code).
    Q_cross<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
    Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  } else
  {
    # 16.08.2012
    Q_cross<-NULL
  }
  for (i in 1:L)
  {
    # For symmetric filters or any historical filter with Lag>0 the decay must be symmetric about b_max(0,Lag)
    # lambda_decay is a 2-dim vector: the first component controls for the exponential decay and the second accounts for the strength of the regularization
    # The maximal weight is limited to 1e+4
    Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(i-1-max(0,Lag))),1e+4)
    # Original idea: with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale: this is not a clever idea...
    #   Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+1]-max(0,Lag))),1e+4)
    
    if (L>4)
    {
      if(i==1)
      {
        Q_smooth[i,i:(i+2)]<-c(1,-2,1)
      } else
      {
        if(i==2)
        {
          Q_smooth[i,(i-1):(i+2)]<-c(-2,5,-4,1)
        } else
        {
          if(i==L)
          {
            Q_smooth[i,(i-2):i]<-c(1,-2,1)
          } else
          {
            if(i==L-1)
            {
              Q_smooth[i,(i-2):(i+1)]<-c(1,-4,5,-2)
            } else
            {
              Q_smooth[i,(i-2):(i+2)]<-c(1,-4,6,-4,1)
            }
          }
        }
      }
    } else
    {
      print("L<=4: no smoothness regularization imposed!!!!!!!!!")
    }
  }
  
  if (length(weight_h_exp[1,])>1)
  {
    
    for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
    {
      if (L>4)
        Q_smooth[j*L+1:L,j*L+1:L]<-Q_smooth[1:L,1:L]
      for (i in 1:L)
        Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(i-1-max(0,Lag))),1e+4)
      # Original idea: with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale: this is not a clever idea...
      #      Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+j+1]-max(0,Lag))),1e+4)
      
    }
    Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
    if (L>1)
    {
      diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
      for (i in 2:length(weight_h_exp[1,]))   #i<-2
      {
        diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
        diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
        diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
      }
    }
    Q_centraldev_original<-solve(Q_centraldev_original)
    
    # 06.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
    #   If grand_mean==T then the code replicates I-MDFA as released prior 30.07.2012
    #   If grand_mean==F then the new parametrization is used.
    #   Differences between both approaches: see section 7.2 of my elements paper posted on SEFBlog (both codes are identical when no regularization is imposed. Otherwise the later version (grand_mean==F) is logically more consistent becuase it treats all series identically (no asymmetry)).
    if (grand_mean)
    {
      diag(Q_cross[L+1:((length(weight_h_exp[1,])-1)*L),L+1:((length(weight_h_exp[1,])-1)*L)])<-
        rep(1,((length(weight_h_exp[1,])-1)*L))
    } else
    {
      #30.07.2012:new definition (parametrization) of Q_cross (Lambda_{cross} in the elements-paper)
      diag(Q_cross)<-1
      for (i in 1:length(weight_h_exp[1,]))
      {
        for (j in 1:L)
        {
          Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]<-Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]-1/length(weight_h_exp[1,])
        }
      }
    }
  } else
  {
    # define matrix for univariate case
    Q_centraldev_original<-NULL
  }
  # Normalizing the troika are new: disentangle the effect by L
  Q_decay<-Q_decay*lambda_decay[2]
  Q_cross<-Q_cross*lambda_cross                   #Qh<-Q_cross
  Q_smooth<-Q_smooth*lambda_smooth
  if (lambda_decay[2]>0)
  {
    # The second parameter in lambda_decay accounts for the strength of the regularization
    Q_decay<-lambda_decay[2]*(Q_decay/(sum(diag(Q_decay))))
  }
  if (lambda_cross>0)
  {
    if (c_eta)
    {
      Q_cross<-lambda_cross^2*(Q_cross/(sum(diag(Q_cross))))
    } else
    {
      Q_cross<-lambda_cross*(Q_cross/(sum(diag(Q_cross))))
    }
  }
  if (lambda_smooth>0&L>4)
  {
    Q_smooth<-lambda_smooth*(Q_smooth/(sum(diag(Q_smooth))))
  }
  return(list(Q_smooth=Q_smooth,Q_decay=Q_decay,Q_cross=Q_cross))
}

centraldev_original_func<-function(L,weight_h_exp)
{
  # Grand-mean parametrization: des_mat is implemented in terms of grand-mean and Q_centraldev_original can
  # be used to invert back to original coefficients
  Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  
  
  Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
  if (L>1)
  {
    diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
    for (i in 2:length(weight_h_exp[1,]))   #i<-2
    {
      diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
      diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
      diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
    }
  }
  Q_centraldev_original<-solve(Q_centraldev_original)
  return(list(Q_centraldev_original=Q_centraldev_original))
}

###################################
### Daten Download und Manipulation
###################################

### Daten
## Load
Load.Daten <- function(ticker,date1,date2){
  
  Stock <- load.Data(ticker,date1,date2)
  return(Stock)
}

## Manipulate
Choice <- function(data,Frequence,OOS = FALSE,date2 = NULL){
  
  data <- data[c("date","adjusted")]
  
  data <- data[complete.cases(data),]
  
  # Select Data by Frequence
  if(Frequence == "4-Daily"){
    # We take always the last available date of each month
    data <- data[rev(seq(nrow(data),1,-3)),]
  }
  
  if(Frequence == "3-Daily"){
    data <- data[rev(seq(nrow(data),1,-2)),]
  }
  
  if(Frequence == "2-Daily"){
    data <- data[rev(seq(nrow(data),1,-1)),]
  }
  
  if(OOS){
    data <- subset(data,date<=date2)
  }
  
  # Select Data by Duration & Explanatory Variables
  # Prepare Data for calculations
  x_md <- as.data.frame(data[,2])
  
  # Shifted and scaled Log-Returns
  x_md <- na.locf(x_md)
  n <- nrow(x_md)
  logReturn <- log(t(t(data.matrix(x_md[-1,])))/t(t(data.matrix(x_md[-n,]))))
  
  
  list(data,logReturn)
}

#####################################
### Univariate Direct Filter Approach
#####################################

### Output of Univariate Direct Filter
Univariate.Filter.cust <- function(data,Lag,L,k,lambda,eta,cutoff,Boolean){
  
  if("i1" %in% Boolean){
    i1 <- TRUE
  } else {
    i1 <- FALSE
  }
  
  if("i2" %in% Boolean){
    i2 <- TRUE
  } else {
    i2 <- FALSE
  }
  
  leng1 <- length(data)
  
  periodogram<-matrix(ncol=1,nrow=ifelse(leng1 %% 2 == 0,leng1/2+1,(leng1-1)/2+1))
  yhat<-rep(NA,leng1)
  xi <- data
  
  Gamma <- c(1,(1:ifelse(leng1 %% 2 == 0,(leng1)/2,(leng1-1)/2))<(leng1+1)/2*k)
  b <- matrix(nrow=L,ncol=1)
  # Compute real-time filters
  # Compute the periodogram based on the data (length 120)
  DFT_L <- per(xi,F)
  periodogram[,1]<-DFT_L$per
  DFT_M <- DFT_L$DFT
  # Optimize filters
  filt<-dfa_analytic(L,lambda,periodogram[,1],Lag,Gamma,eta,cutoff,i1,i2)
  for (j in (L+1):(leng1)){
    yhat[j]<-filt$b%*%xi[(j+Lag):(j-L+1+Lag)]}
  trffkt <- filt$trffkt
  ABS <- abs(trffkt)
  
  omega_k<-pi*0:(leng1/2)/(leng1/2)
  Shift <- Arg(trffkt)/omega_k
  
  DFT_M <- cbind(DFT_M,DFT_M)
  trffkt <- as.matrix(trffkt,ncol = 1)
  
  list(yhat,ABS,Shift,trffkt,Gamma,DFT_M,cutoff,filt$b)
  
}

# Amplitude
Amplitude <- function(Amp,k){
  leng1 <- length(Amp)*2
  p <- plot_ly() %>%
    layout(yaxis = list(
      title = ""),
      xaxis = list(
        ticktext = list("0","pi/6","2pi/6","3pi/6",
                        "4pi/6","5pi/6","pi"), 
        tickvals = as.list((1+0:6*leng1/12)/length(Amp)*pi),
        tickmode = "array"
      )) %>%
    add_lines(x = seq(1,length(Amp))/length(Amp)*pi,y = Amp,line = list(color = "black"), name = "Amplitude") %>% 
    add_lines(x = seq(1,length(c(rep(1,round(leng1/2*k)),rep(0,round(leng1/2*(1-k))))))/length(c(rep(1,round(leng1/2*k)),rep(0,round(leng1/2*(1-k)))))*pi,y = c(rep(1,round(leng1/2*k)),rep(0,round(leng1/2*(1-k)))),line = list(color = "red"), name = "Target") 
}

# Shift
Shift <- function(shf){
  leng1 <- length(shf)*2
  p <- plot_ly() %>%
    layout(yaxis = list(
      title = ""),
      xaxis = list(
        ticktext = list("0","pi/6","2pi/6","3pi/6",
                        "4pi/6","5pi/6","pi"), 
        tickvals = as.list((1+0:6*leng1/12)/length(shf)*pi),
        tickmode = "array"
      )) %>%
    add_lines(x = seq(1,length(shf))/length(shf)*pi,y = shf,line = list(color = "black"), name = "Shift") %>% 
    add_lines(x = seq(1,length(shf))/length(shf)*pi,y = rep(0,round(leng1/2)),line = list(color = "red"), name = "Target") 
}

### Double Filter Univariate Direct Filter
Univariate.Filter.cust.conv <- function(data,Lag,L,k,lambda,eta,cutoff,Boolean,
                                        Lag1,L1,k1,lambda1,eta1,cutoff1,Boolean1){
  
  if("i1" %in% Boolean){
    i1 <- TRUE
  } else {
    i1 <- FALSE
  }
  
  if("i2" %in% Boolean){
    i2 <- TRUE
  } else {
    i2 <- FALSE
  }
  
  if("i1" %in% Boolean1){
    i11 <- TRUE
  } else {
    i11 <- FALSE
  }
  
  if("i2" %in% Boolean1){
    i21 <- TRUE
  } else {
    i21 <- FALSE
  }
  
  leng1 <- length(data)
  
  periodogram<-matrix(ncol=1,nrow=ifelse(leng1 %% 2 == 0,leng1/2+1,(leng1-1)/2+1))
  yhat<-rep(NA,leng1)
  yhatc<-rep(NA,leng1)
  xi <- data
  
  Gamma <- c(1,(1:ifelse(leng1 %% 2 == 0,(leng1)/2,(leng1-1)/2))<(leng1+1)/2*k)
  Gamma1 <- c(1,(1:ifelse(leng1 %% 2 == 0,(leng1)/2,(leng1-1)/2))<(leng1+1)/2*k1)
  # Compute real-time filters
  # Compute the periodogram
  DFT_L <- per(xi,F)
  periodogram[,1]<-DFT_L$per
  DFT_M <- DFT_L$DFT
  # Optimize filters
  filt<-dfa_analytic(L,lambda,periodogram[,1],Lag,Gamma,eta,cutoff,i1,i2)
  filt1<-dfa_analytic(L1,lambda1,periodogram[,1],Lag1,Gamma1,eta1,cutoff1,i11,i21)
  n <- max(length(filt1$b),length(filt$b))
  b1 <- c(filt$b,rep(0,n-length(filt$b)))
  b2 <- c(filt1$b,rep(0,n-length(filt1$b)))
  #b <- as.numeric(fft(fft(b1)*fft(b2),inverse = T))
  b <- convolve(filt$b,rev(filt1$b),type = "o")
  #yhat <- convolve(convolve(xi,rev(filt$b),type = "f"),rev(filt1$b),type = "f")
  #for (j in (L+1):(leng1)){
  #  yhatc[j]<-filt$b%*%xi[(j+Lag):(j-L+1+Lag)]}
  #for (j in (L1+1):(leng1)){
  #  yhat[j]<-filt1$b%*%yhatc[(j+Lag1):(j-L1+1+Lag1)]}
  
  for (j in (L+L1):(leng1)){
    yhat[j]<-b%*%xi[(j+Lag):(j-(L+L1-1)+1+Lag)]}
  #yhatc <- stats::filter(xi,filt$b,"convolution",1)
  #yhat <- stats::filter(yhatc,filt1$b,"convolution",1)
  
  trffkt <- filt$trffkt*filt1$trffkt
  ABS <- abs(trffkt)
  
  omega_k<-pi*0:(leng1/2)/(leng1/2)
  Shift <- Arg(trffkt)/omega_k
  
  DFT_M <- cbind(DFT_M,DFT_M)
  trffkt <- as.matrix(trffkt,ncol = 1)
  
  list(yhat,ABS,Shift,trffkt,Gamma,DFT_M,cutoff,b)
  
}

### Chart with Out of Sample Periode DFA
Out_of_Sample_Chart <- function(x1,y1,date,filter.data,n){
  
  leng <- length(x1)
  leng1 <- length(filter.data)
  data <- rep(NA,leng)
  for (j in (leng1+1):(leng)){
    data[j]<-filter.data%*%y1[(j):(j-(leng1)+1)]}
  
  p <- plot_ly(x = x1, y = y1,
               type = "scatter", mode = "lines", name = "Underlying Data",line = list(color = "black")) %>%
    layout(yaxis = list(
      title = "",
      showticklabels = FALSE
    ),legend = list(orientation = 'h')) %>%
    add_lines(y = Skalieren(data,n), name = "Direct Filter RT", line = list(color = "#ff7f0e")) %>%
    add_segments(x = as.Date(as.Date(date)+1), 
                 xend = as.Date(as.Date(date)+1), 
                 y = min(y1), yend = max(y1), 
                 line = list(color = "blue"), opacity = 0.5, name = "Start of Out of Sample")
  
  return(p)
}

#########################
### Explanatory Variables
#########################

### Explanatory Plot
Exp.Plot <- function(data, Cur){
  
  Currency.Vect <- c("CHF","USD","EUR")
  Currency <- Currency.Vect[which(Currency.Vect != Cur)]
  
  time <- data[[1]][,1]
  data <- data[[1]][,-1]
  
  ay <- list(
    tickfont = list(color = "red"),
    overlaying = "y",
    side = "right",
    title = "PMI"
  )
  
  
  p <- plot_ly(x = time, y = data[,1],
               type = 'scatter', mode = 'lines', name = "Target") %>% 
    layout(xaxis = list(title = "Time"),yaxis2 = ay) %>% 
    add_lines(y = data[,2], name = "PMI", yaxis = "y2", opacity = 0.7) %>% 
    add_lines(y = data[,3], name = "Inflation", opacity = 0.7) %>% 
    add_lines(y = data[,4], name = "Short-Duration Yield", opacity = 0.7) %>% 
    add_lines(y = data[,5], name = Currency[1], opacity = 0.7) %>% 
    add_lines(y = data[,6], name = Currency[2], opacity = 0.7)
  
  p
  
}

### Cross-Correlation Explanatory Variables
Cross.Corr.Exp <- function(data,Cur,Freq){
  
  Currency.Vect <- c("CHF","USD","EUR")
  Currency <- Currency.Vect[which(Currency.Vect != Cur)]
  
  ccf.data.1 <- ccf(data[,3],data[,2],lag.max = 100,plot = FALSE)
  ccf.data.2 <- ccf(data[,4],data[,2],lag.max = 100,plot = FALSE)
  ccf.data.3 <- ccf(data[,5],data[,2],lag.max = 100,plot = FALSE)
  ccf.data.4 <- ccf(data[,6],data[,2],lag.max = 100,plot = FALSE)
  ccf.data.5 <- ccf(data[,7],data[,2],lag.max = 100,plot = FALSE)
  
  Best.1 <- which(as.numeric(ccf.data.1$acf) == max(as.numeric(ccf.data.1$acf),na.rm = T))
  Best.2 <- which(as.numeric(ccf.data.2$acf) == max(as.numeric(ccf.data.2$acf),na.rm = T))
  Best.3 <- which(as.numeric(ccf.data.3$acf) == max(as.numeric(ccf.data.3$acf),na.rm = T))
  Best.4 <- which(as.numeric(ccf.data.4$acf) == max(as.numeric(ccf.data.4$acf),na.rm = T))
  Best.5 <- which(as.numeric(ccf.data.5$acf) == max(as.numeric(ccf.data.5$acf),na.rm = T))
  
  mat <- rbind(c(as.numeric(ccf.data.1$lag)[Best.1],as.numeric(ccf.data.2$lag)[Best.2],as.numeric(ccf.data.3$lag)[Best.3],
                 as.numeric(ccf.data.4$lag)[Best.4],as.numeric(ccf.data.5$lag)[Best.5]),
               c(as.numeric(ccf.data.1$acf)[Best.1],as.numeric(ccf.data.2$acf)[Best.2],as.numeric(ccf.data.3$acf)[Best.3],
                 as.numeric(ccf.data.4$acf)[Best.4],as.numeric(ccf.data.5$acf)[Best.5]))
  
  colnames(mat) <- c("PMI","Inflation","Short-Duration Yield",Currency[1],Currency[2])
  rownames(mat) <- c("Lag","Correlation")
  
  p <- plot_ly(x = as.numeric(ccf.data.1$lag), y = as.numeric(ccf.data.1$acf),
               type = "bar", name = "PMI") %>%
    layout(legend = list(),  xaxis = list(title = 'Lag'), yaxis = list(title = 'Correlation')) %>%
    add_bars(y = as.numeric(ccf.data.2$acf),name = 'Inflation') %>%
    add_bars(y = as.numeric(ccf.data.3$acf),name = 'Short-Duration Yield') %>%
    add_bars(y = as.numeric(ccf.data.4$acf),name = Currency[1]) %>%
    add_bars(y = as.numeric(ccf.data.5$acf),name = Currency[2]) 
  
  Freq1 <- ifelse(Freq == "T?glich","Tage",ifelse(Freq == "Monatlich","Monate","Wochen"))
  
  if(Freq == "2-W?chentlich"){
    mat[1,] <- 2*mat[1,]
    mat[1,] <- paste(mat[1,],Freq1)
  } else {
    mat[1,] <- paste(mat[1,],Freq1)
  }
  
  mat[2,] <- round(as.numeric(mat[2,]),2)
  
  list(p,mat)
}

#######################################
### Multivariate Direct Filter Approach
#######################################

### Output of Multivariate Direct Filter
Multivariate.Filter.cust <- function(data,Lag,L,k,lambda,eta,cutoff,Boolean){
  
  if("i1" %in% Boolean){
    i1 <- TRUE
  } else {
    i1 <- FALSE
  }
  
  if("i2" %in% Boolean){
    i2 <- TRUE
  } else {
    i2 <- FALSE
  }
  
  Spect <- spec_comp(insamp = nrow(data),x=data,d = 0)
  weight_func <- Spect$weight_func 
  grand_mean <- F
  c_eta <- F
  weight_structure<-c(0,0)
  troikaner <- F
  weight_constraint <- 0
  lin_eta<-F
  weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
  #weight_constraint<-seq(1,1/(ncol(weight_func)-1),-1/ncol(weight_func)-1)
  lambda_cross<-0
  lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  shift_constraint<-rep(0,nrow(weight_func)-1)
  grand_mean<-F
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  b0_H0<-NULL
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  troikaner<-F
  leng1 <- nrow(data)
  
  yhat<-rep(NA,leng1)
  xi <- data
  
  Gamma <- c(1,(1:ifelse(leng1 %% 2 == 0,(leng1)/2,(leng1-1)/2))<(leng1+1)/2*k)
           
  b <- matrix(nrow=L,ncol=1)
  # Compute real-time filters
  # Compute the periodogram based on the data (length 120) 
  
  # Optimize filters
  filt<-mdfa_analytic(L,lambda,Spect$weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)
  
  for (j in (L+1):(leng1)){
    yhat[j]<-sum(apply(filt$b*xi[(j+Lag):(j-L+1+Lag),2:ncol(xi)],1,sum))}
  trffkt <- filt$trffkt
  ABS <- abs(trffkt)
  
  Shift <- Arg(t(sign(apply(filt$b,2,sum))*t(trffkt)))/((0:(nrow(trffkt)-1))*pi/(nrow(trffkt)-1))
  Shift[1,] <- apply(filt$b*((0:(L-1))),2,sum)/apply(filt$b,2,sum)
  list(yhat,ABS,Shift,filt$trffkt,Gamma, weight_func,cutoff,filt$b)
}

# Amplitude MDFA
Amplitude.MDF <- function(Amp,k){
  leng1 <- nrow(Amp)*2
  p <- plot_ly() %>%
    layout(yaxis = list(
      title = ""),
      xaxis = list(
        ticktext = list("0","pi/6","2pi/6","3pi/6",
                        "4pi/6","5pi/6","pi"), 
        tickvals = as.list((1+0:6*leng1/12)/nrow(Amp)*pi),
        tickmode = "array"
      )) %>%
    add_lines(x = seq(1,nrow(Amp))/nrow(Amp)*pi,y = Amp[,1],line = list(color = "black"), name = "Yield") %>% 
    add_lines(x = seq(1,nrow(Amp))/nrow(Amp)*pi,y = c(rep(1,round(leng1/2*k)),rep(0,round(leng1/2*(1-k)))),line = list(color = "red"), name = "Target")
  if(ncol(Amp)>1){
    for (i in 2:ncol(Amp)) {
      p <- p %>% 
        add_lines(x = seq(1,nrow(Amp))/nrow(Amp)*pi,y = Amp[,i],name = paste("Explanatory",i-1)) 
    }
  }
  p
}

# Shift MDFA
Shift.MDF <- function(shf){
  leng1 <- nrow(shf)*2
  p <- plot_ly() %>%
    layout(yaxis = list(
      title = ""),
      xaxis = list(
        ticktext = list("0","pi/6","2pi/6","3pi/6",
                        "4pi/6","5pi/6","pi"), 
        tickvals = as.list((1+0:6*leng1/12)/nrow(shf)*pi),
        tickmode = "array"
      )) %>%
    add_lines(x = seq(1,nrow(shf))/nrow(shf)*pi,y = shf[,1],line = list(color = "black"), name = "Yield") %>% 
    add_lines(x = seq(1,nrow(shf))/nrow(shf)*pi,y = rep(0,round(leng1/2)),line = list(color = "red"), name = "Target") 
  
  if(ncol(shf)>1){
    for (i in 2:ncol(shf)) {
      p <- p %>% 
        add_lines(x = seq(1,nrow(shf))/nrow(shf)*pi,y = shf[,i],name = paste("Explanatory",i-1)) 
    }
  }
  p
}

### Double Filter Multivariate Direct Filter
Multivariate.Filter.cust.conv <- function(data,Lag,L,k,lambda,eta,cutoff,Boolean,
                                          Lag1,L1,k1,lambda1,eta1,cutoff1,Boolean1){
  
  if("i1" %in% Boolean){
    i1 <- TRUE
  } else {
    i1 <- FALSE
  }
  
  if("i2" %in% Boolean){
    i2 <- TRUE
  } else {
    i2 <- FALSE
  }
  
  if("i1" %in% Boolean1){
    i11 <- TRUE
  } else {
    i11 <- FALSE
  }
  
  if("i2" %in% Boolean1){
    i21 <- TRUE
  } else {
    i21 <- FALSE
  }
  
  leng2 <- nrow(data)
  
  Spect <- spec_comp(insamp = nrow(data),x=data,d = 0)
  weight_func <- Spect$weight_func 
  
  grand_mean <- F
  c_eta <- F
  weight_structure<-c(0,0)
  troikaner <- F
  weight_constraint <- 0
  lin_eta<-F
  weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
  #weight_constraint<-seq(1,1/(ncol(weight_func)-1),-1/ncol(weight_func)-1)
  lambda_cross<-0
  lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  shift_constraint<-rep(0,nrow(weight_func)-1)
  grand_mean<-F
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  b0_H0<-NULL
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  lag_mat1<-matrix(rep(0:(L1-1),ncol(weight_func)),nrow=L1)
  troikaner<-F
  
  yhat<-rep(NA,leng2)
  
  xi <- data
  yhat2<-matrix(NA,nrow = leng2,ncol = ncol(xi)-1)
  Gamma <- c(1,(1:ifelse(leng2 %% 2 == 0,(leng2)/2,(leng2-1)/2))<(leng2+1)/2*k)
  Gamma1 <- c(1,(1:ifelse(leng2 %% 2 == 0,(leng2)/2,(leng2-1)/2))<(leng2+1)/2*k1)
  
  # Optimize filters
  filt<-mdfa_analytic(L,lambda,Spect$weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                      lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                      b0_H0,c_eta,weight_structure,white_noise,
                      synchronicity,lag_mat,troikaner)
  filt1<-mdfa_analytic(L1,lambda1,Spect$weight_func,Lag1,Gamma1,eta1,cutoff1,i11,i21,weight_constraint,
                       lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                       b0_H0,c_eta,weight_structure,white_noise,
                       synchronicity,lag_mat1,troikaner)
  
  n <- max(nrow(filt1$b),nrow(filt$b))
  b1 <- rbind(filt$b,matrix(0,nrow = n-nrow(filt$b),ncol = ncol(filt$b)))
  b2 <- rbind(filt1$b,matrix(0,nrow = n-nrow(filt1$b),ncol = ncol(filt1$b)))
  b <- matrix(0,nrow = (nrow(filt$b)+nrow(filt1$b)-1),ncol = ncol(filt$b))
  for(i in 1:ncol(filt$b)){
    b[,i] <- convolve(filt$b[,i],rev(filt1$b[,i]),type = "o")
  }
  
  for (j in (L1+L):(leng2)){
    yhat[j]<-sum(apply(b*xi[(j+Lag):(j-(L1+L-1)+1+Lag),2:ncol(xi)],1,sum))}
  
  trffkt <- filt$trffkt*filt1$trffkt
  ABS <- abs(trffkt)
  
  Shift <- Arg(t(sign(apply(b,2,sum))*t(trffkt)))/((0:(nrow(trffkt)-1))*pi/(nrow(trffkt)-1))
  Shift[1,] <- apply(b*((0:(L+L1-1-1))),2,sum)/apply(b,2,sum)
  trffkt <- as.matrix(trffkt,ncol = ncol(b))
  
  list(yhat,ABS,Shift,trffkt,Gamma, weight_func,cutoff,b)
}

### Chart with Out of Sample Periode MDFA
Out_of_Sample_Chart_M <- function(x1,y1,date,filter.data,n){
  
  leng <- length(x1)
  leng1 <- nrow(filter.data)
  data <- rep(NA,leng)
  for (j in (leng1+1):(leng)){
    data[j]<-sum(apply(filter.data*y1[(j):(j-(leng1)+1),2:ncol(y1)],1,sum))}
  
  p <- plot_ly(x = x1, y = y1[,1],
               type = "scatter", mode = "lines", name = "Underlying Data",line = list(color = "black")) %>%
    layout(showlegend = FALSE,yaxis = list(
      title = "",
      showticklabels = FALSE
    )) %>%
    add_lines(y = Skalieren(data,n), name = "Direct Filter RT", line = list(color = "#ff7f0e")) %>%
    add_segments(x = as.Date(as.Date(date)+1), 
                 xend = as.Date(as.Date(date)+1), 
                 y = min(y1[,1]), yend = max(y1[,1]), 
                 line = list(color = "blue"), opacity = 0.5, name = "Start of Out of Sample")
  
  return(p)
}

#####################
### Error Computation
#####################

# Computing Errors (ATS)
MS_decomp_total<-function(Gamma,trffkt,weight_func,cutoff)
{
  Lag <- 0
  if (!(length(trffkt[,1])==length(weight_func[,1])))
  {
    len_w<-min(length(trffkt[,1]),length(weight_func[,1]))
    if (length(trffkt[,1])<length(weight_func[,1]))
    {
      len_r<-(length(weight_func[,1])-1)/(length(trffkt[,1])-1)
      weight_funch<-weight_func[c(1,(1:(len_w-1))*len_r),]
      trffkth<-trffkt
    } else
    {
      len_r<-1/((length(weight_func[,1])-1)/(length(trffkt[,1])-1))
      trffkth<-trffkt[c(1,(1:(len_w-1))*len_r),]
      weight_funch<-weight_func
    }
  } else
  {
    len_w<-length(trffkt[,1])
    weight_funch<-weight_func
    trffkth<-trffkt
    Gammah<-Gamma
  }
  if (length(Gamma)>len_w)
  {
    len_r<-(length(Gamma)-1)/(len_w-1)
    Gammah<-Gamma[c(1,(1:(len_w-1))*len_r)]
  }
  
  
  
  weight_h<-weight_funch
  K<-length(weight_funch[,1])-1
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))
  # DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])
  
  trt<-apply(((trffkth)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
  # MS-filter error : DFA-criterion without effects by lambda or eta (one must divide spectrum by eta_vec)
  MS_error<-sum((abs(Gammah*weight_target-trt))^2)/(2*(K+1)^2)
  Gamma_cp<-Gammah[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gammah[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-trt[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-trt[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-weight_target[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-weight_target[(2+as.integer(K*(cutoff/pi))):(K+1)]
  # define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2,na.rm = T)/(2*(K+1)^2)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp,na.rm = T)/(2*(K+1)^2)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2,na.rm = T)/(2*(K+1)^2)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)/(2*(K+1)^2)
  
  Mat <- as.data.frame(matrix(c(Accuracy,Smoothness, Timeliness, MS_error),ncol = 4,nrow = 1))
  colnames(Mat) <- c("Accuracy","Smoothness","Timeliness","MS-Error")
  return(Mat)
}

####################
### Hilfs-Funktionen
####################

# Scale
Skalieren <- function(y,n=1){
  y <- y/max(abs(y),na.rm = T)*n
}

# Mean of Filters
Mean.of.Vectors <- function(d.x,w.x,w2.x,m.x,d.y,w.y,w2.y,m.y){
  m <- rep(NA,length(as.numeric(d.y)))
  m[which(as.character(d.x) %in% as.character(m.x))] <- m.y
  m <- na.approx(m)
  m <- c(rep(NA,length(as.numeric(d.y))-length(m)),m)
  
  w <- rep(NA,length(as.numeric(d.y)))
  w[which(as.character(d.x) %in% as.character(w.x))] <- w.y
  w <- na.approx(w)
  w <- c(rep(NA,length(as.numeric(d.y))-length(w)),w)
  
  w2 <- rep(NA,length(as.numeric(d.y)))
  w2[which(as.character(d.x) %in% as.character(w2.x))] <- w2.y
  w2 <- na.approx(w2)
  w2 <- c(rep(NA,length(as.numeric(d.y))-length(w2)),w2)
  
  mv <- rowMeans(cbind(m,w,w2,as.numeric(d.y)), na.rm=FALSE)
  list(mv,m,w2,w,as.numeric(d.y))
}

#################################
### Filter Performance Funktionen
#################################

### Trend Visualisierung
# All Out of Sample Performance Plots
Segmente.Trend <- function(x,y1,y,date,filter.data,OS,Variance = 0,Fee = FALSE,Fee.num = 1){
  
  Fee.num = Fee.num/100
  
  if(!is.matrix(filter.data)){
    leng <- length(x)-1
    leng1 <- length(filter.data)
    data <- rep(NA,leng)
    for (j in (leng1+1):(leng)){
      data[j]<-filter.data%*%y1[(j):(j-(leng1)+1)]}
  } else {
    leng <- length(x)-1
    leng1 <- nrow(filter.data)
    data <- rep(NA,leng)
    for (j in (leng1+1):(leng)){
      data[j]<-sum(apply(filter.data*y1[(j):(j-leng1+1),2:ncol(y1)],1,sum))}
  }
  
  
  res <- c(NA,data)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  
  trend <- c(NA,data)
  trend[trend == 0] <- NA
  trend_min <- c(NA,trend)[1:(length(trend))]
  trend1 <- ifelse(trend>0 | (trend < 0 & trend_min > 0),1,0)
  trend2 <- ifelse(trend<0 | (trend > 0 & trend_min < 0),-1,0)
  
  ind <- which(trend1==1)
  ind1 <- which(trend2==-1)
  
  indx <- which(res1==1)
  indx1 <- which(res1==-1)
  
  x.trend1 <- as.Date(x)[indx]
  x.trend2 <- as.Date(x)[indx1]
  
  if(x.trend1[1]<x.trend2[1]){
    x.trend1.new <- x.trend1[-1]
  } else {
    x.trend1.new <- x.trend1
  }
  
  if(x.trend1[length(x.trend1)]<x.trend2[length(x.trend2)]){
    x.trend1.new <- c(x.trend1.new,as.Date(x)[length(as.character(x))])
  } else {
    x.trend1.new <- x.trend1.new
  }
  
  if(x.trend1[1]>x.trend2[1]){
    x.trend2.new <- x.trend2[-1]
  } else {
    x.trend2.new <- x.trend2
  }
  
  if(x.trend1[length(x.trend1)]>x.trend2[length(x.trend2)]){
    x.trend2.new <- c(x.trend2.new,as.Date(x)[length(as.character(x))])
  } else {
    x.trend2.new <- x.trend2.new
  }
  
  
  line <- list(
    type = "rect",
    fillcolor = "green",
    line = list(color = "green", width = 0.1),
    opacity = 0.3,
    xref = "x",
    yref = "y",
    name = "Positive Trend"
  )
  
  lines <- list()
  for (i in 1:length(x.trend1)) {
    line[["x0"]] <- x.trend1[i]
    line[["x1"]] <- x.trend2.new[i]
    line[c("y0", "y1")] <- c(min(y), max(y))
    lines <- c(lines, list(line))
  }
  
  line1 <- list(
    type = "rect",
    fillcolor = "red",
    line = list(color = "red", width = 0.1),
    opacity = 0.3,
    xref = "x",
    yref = "y",
    name = "Negative Trend"
  )
  
  lines1 <- list()
  for (i in 1:length(x.trend2)) {
    line1[["x0"]] <- x.trend2[i]
    line1[["x1"]] <- x.trend1.new[i]
    line1[c("y0", "y1")] <- c(min(y), max(y))
    lines1 <- c(lines1, list(line1))
  }
  
  Ret.per <- CalculateReturns(ts(y), method = "discrete")
  Ret.per <- Ret.per[-1]
  ynew.perf <- cumprod(c(100,1+Ret.per))
  
  Ret.per2 <- as.data.frame(cbind(x[-1],Ret.per))
  colnames(Ret.per2) <- c("date","Return")
  Ret.per2$date <- as.Date(Ret.per2$date)
  Ret.per3 <- Ret.per2
  
  for (i in 1:length(x.trend2)) {
    
    Ret.per2$Return[which(Ret.per2$date > x.trend2[i] & Ret.per2$date <= x.trend1.new[i])] <- 0
    
    
    
  }
  
  ynew.perf2 <- cumprod(c(100,1+Ret.per2$Return))
  
  for (i in 1:length(x.trend2)) {
    if(as.Date(x.trend1.new[i])>as.Date(date)){
      dat.interim <- max(c(x.trend2[i],as.Date(date)))
      Ret.per3$Return[which(Ret.per3$date > as.Date(dat.interim) & Ret.per3$date <= x.trend1.new[i])] <- 0
    }
  }
  
  
  ynew.perf3 <- cumprod(c(100,1+Ret.per3$Return))
  
  indx <- which(Ret.per3$date %in% c(x.trend2,x.trend1) & Ret.per3$date >= as.Date(date))
  Ret.per4 <- Ret.per3$Return
  Ret.per4[indx] <- Ret.per4[indx]-Fee.num
  ynew.perf4 <- cumprod(c(100,1+Ret.per4))
  
  if(Fee){
    ynew.perf3 <- ynew.perf4
  }
  
  ynew.perf3 <- (ynew.perf3-100)/100
  ynew.perf2 <- (ynew.perf2-100)/100
  ynew.perf <- (ynew.perf-100)/100
  p.perf <- plot_ly(x = x, y = ynew.perf,
               type = "scatter", mode = "lines", name = "Hold",line = list(color = "black"))
  
  p.perf <- p.perf %>% add_lines(y = ynew.perf2, name = "DFA Strategy",line = list(color = "green"))
  
  
  outperf1 <- plot_ly(x = x, y = ynew.perf2-ynew.perf,
                      type = "bar", name = "Hold",line = list(color = "black"))
  
  if(date != 0){
    p.perf <- p.perf %>%
      add_segments(x = as.Date(as.Date(date)+1), 
                   xend = as.Date(as.Date(date)+1), 
                   y = min(c(ynew.perf,ynew.perf2)), yend = max(c(ynew.perf,ynew.perf2)), 
                   line = list(color = "blue"), opacity = 1, name = "Start of Out of Sample")
  }
  
  p.perf2 <- plot_ly(x = x, y = ynew.perf,
                    type = "scatter", mode = "lines", name = "Hold",line = list(color = "black"))
  
  p.perf2 <- p.perf2 %>% add_lines(y = ynew.perf3, name = "DFA Strategy",line = list(color = "green"))
  
  outperf2 <- plot_ly(x = x, y = ynew.perf3-ynew.perf,
                      type = "bar", name = "Hold",line = list(color = "black"))
  
  if(date != 0){
    p.perf2 <- p.perf2 %>%
      add_segments(x = as.Date(as.Date(date)+1), 
                   xend = as.Date(as.Date(date)+1), 
                   y = min(c(ynew.perf,ynew.perf3)), yend = max(c(ynew.perf,ynew.perf3)), 
                   line = list(color = "blue"), opacity = 1, name = "Start of Out of Sample")
  }
  
  lines4 <- list()
  for (i in 1:length(x.trend1)) {
    line[["x0"]] <- x.trend1[i]
    line[["x1"]] <- x.trend2.new[i]
    line[c("y0", "y1")] <- c(min(c(ynew.perf,ynew.perf2)), max(c(ynew.perf,ynew.perf2)))
    lines4 <- c(lines4, list(line))
  }
  
  lines5 <- list()
  for (i in 1:length(x.trend2)) {
    line1[["x0"]] <- x.trend2[i]
    line1[["x1"]] <- x.trend1.new[i]
    line1[c("y0", "y1")] <- c(min(c(ynew.perf,ynew.perf2)), max(c(ynew.perf,ynew.perf2)))
    lines5 <- c(lines5, list(line1))
  }
  
  lines6 <- list()
  for (i in 1:length(x.trend1)) {
    line[["x0"]] <- x.trend1[i]
    line[["x1"]] <- x.trend2.new[i]
    line[c("y0", "y1")] <- c(min(c(ynew.perf,ynew.perf3)), max(c(ynew.perf,ynew.perf3)))
    lines6 <- c(lines6, list(line))
  }
  
  lines7 <- list()
  for (i in 1:length(x.trend2)) {
    line1[["x0"]] <- x.trend2[i]
    line1[["x1"]] <- x.trend1.new[i]
    line1[c("y0", "y1")] <- c(min(c(ynew.perf,ynew.perf3)), max(c(ynew.perf,ynew.perf3)))
    lines7 <- c(lines7, list(line1))
  }
  
  p.perf <- p.perf %>%
    layout(xaxis = list(title = ''), yaxis = list(title = '',tickformat = ".3%"),
           shapes = c(lines5,lines4),legend = list(orientation = 'h'))
  
  p.perf2 <- p.perf2 %>%
    layout(xaxis = list(title = ''), yaxis = list(title = '',tickformat = ".3%"),
           shapes = c(lines6,lines7),legend = list(orientation = 'h'))
  
  
  p <- plot_ly(x = x, y = y,
            type = "scatter", mode = "lines", name = "Underlying Data",line = list(color = "black"))
  
  if(date != 0){
    p <- p %>%
      add_segments(x = as.Date(as.Date(date)+1), 
                   xend = as.Date(as.Date(date)+1), 
                   y = min(y)-0.1, yend = max(y)+0.1, 
                   line = list(color = "blue"), opacity = 1, name = "Start of Out of Sample")
  } else {
    p <- p
  }
  
  y.int <- y
  
  differences <- c(NA,diff(y))
  
  for (i in 1:length(x.trend2)) {
    differences[as.character(x)> x.trend2[i] & as.character(x)<= x.trend1.new[i]] <- -differences[as.character(x)> x.trend2[i] & as.character(x)<= x.trend1.new[i]]
  }
  
  y <- cumsum(c(y[1],differences[-1]))
  
  lines2 <- list()
  for (i in 1:length(x.trend1)) {
    line[["x0"]] <- x.trend1[i]
    line[["x1"]] <- x.trend2.new[i]
    line[c("y0", "y1")] <- c(min(y), max(y))
    lines2 <- c(lines2, list(line))
  }
  
  lines3 <- list()
  for (i in 1:length(x.trend2)) {
    line1[["x0"]] <- x.trend2[i]
    line1[["x1"]] <- x.trend1.new[i]
    line1[c("y0", "y1")] <- c(min(y), max(y))
    lines3 <- c(lines3, list(line1))
  }
  
  
  p1 <- plot_ly(x = x, y = y,
                type = "scatter", mode = "lines", name = "Performance",line = list(color = "black")) %>%
    layout(xaxis = list(title = 'Date'), yaxis = list(title = ''),
           shapes = c(lines2,lines3)) 
  
  if(date != 0){
    p1 <- p1 %>%
      add_segments(x = as.Date(as.Date(date)+1), 
                   xend = as.Date(as.Date(date)+1), 
                   y = min(y)-0.1, yend = max(y)+0.1, 
                   line = list(color = "blue"), opacity = 1, name = "Start of Out of Sample")
  }
  
  
  dates1 <- as.data.frame(as.Date(x.trend1))
  colnames(dates1) <- c("Date")
  dates1$Trend <- rep("Positive",length(x.trend1))
  dates2 <- as.data.frame(as.Date(x.trend2))
  colnames(dates2) <- c("Date")
  dates2$Trend <- rep("Negative",length(x.trend2))
  
  dates <- rbind(dates1,dates2)
  dates <- dates[order(dates$Date),]
  
  y <- as.data.frame(y)
  y$dates <- as.Date(x)
  
  Mat <- 0
  Mat1 <- 0
  
  p2 <- p %>%
    layout(xaxis = list(title = 'Date'), yaxis = list(title = 'Yield'),
           shapes = c(lines1,lines))
  
  p <- p %>%
    layout(showlegend = FALSE,xaxis = list(title = 'Date'), shapes = c(lines1,lines))
  
  
  
  list(p,p1,Mat,Mat1,p2,p.perf,p.perf2,outperf1,outperf2)
}

# Segmente Plot
Only.Segmente <- function(x,data,y,Freq){
  
  if(Freq == "Average"){
    data <- data[[1]]
  }
  if(Freq == "4-Daily"){
    data <- data[[2]]
  }
  if(Freq == "3-Daily"){
    data <- data[[3]]
  }
  if(Freq == "2-Daily"){
    data <- data[[4]]
  }
  if(Freq == "Daily"){
    data <- data[[5]]
  }
  res <- c(NA,data)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  
  trend <- c(NA,data)
  trend[trend == 0] <- NA
  trend_min <- c(NA,trend)[1:(length(trend))]
  trend1 <- ifelse(trend>0 | (trend < 0 & trend_min > 0),1,0)
  trend2 <- ifelse(trend<0 | (trend > 0 & trend_min < 0),-1,0)
  
  ind <- which(trend1==1)
  ind1 <- which(trend2==-1)
  
  indx <- which(res1==1)
  indx1 <- which(res1==-1)
  
  x.trend1 <- as.Date(x)[indx]
  x.trend2 <- as.Date(x)[indx1]
  
  if(x.trend1[1]<x.trend2[1]){
    x.trend1.new <- x.trend1[-1]
  } else {
    x.trend1.new <- x.trend1
  }
  
  if(x.trend1[length(x.trend1)]<x.trend2[length(x.trend2)]){
    x.trend1.new <- c(x.trend1.new,as.Date(x)[length(as.character(x))])
  } else {
    x.trend1.new <- x.trend1.new
  }
  
  if(x.trend1[1]>x.trend2[1]){
    x.trend2.new <- x.trend2[-1]
  } else {
    x.trend2.new <- x.trend2
  }
  
  if(x.trend1[length(x.trend1)]>x.trend2[length(x.trend2)]){
    x.trend2.new <- c(x.trend2.new,as.Date(x)[length(as.character(x))])
  } else {
    x.trend2.new <- x.trend2.new
  }
  
  
  line <- list(
    type = "rect",
    fillcolor = "green",
    line = list(color = "green", width = 0.1),
    opacity = 0.3,
    xref = "x",
    yref = "y",
    name = "Positive Trend"
  )
  
  lines <- list()
  for (i in 1:length(x.trend1)) {
    line[["x0"]] <- x.trend1[i]
    line[["x1"]] <- x.trend2.new[i]
    line[c("y0", "y1")] <- c(min(y), max(y))
    lines <- c(lines, list(line))
  }
  
  line1 <- list(
    type = "rect",
    fillcolor = "red",
    line = list(color = "red", width = 0.1),
    opacity = 0.3,
    xref = "x",
    yref = "y",
    name = "Negative Trend"
  )
  
  lines1 <- list()
  for (i in 1:length(x.trend2)) {
    line1[["x0"]] <- x.trend2[i]
    line1[["x1"]] <- x.trend1.new[i]
    line1[c("y0", "y1")] <- c(min(y), max(y))
    lines1 <- c(lines1, list(line1))
  }
  
  p <- plot_ly(x = x, y = y,
               type = "scatter", mode = "lines", name = "Underlying Data",line = list(color = "black"))
  
  p <- p %>%
    layout(showlegend = FALSE,xaxis = list(title = 'Date'), shapes = c(lines1,lines))
  
  p
}

### Cross-Correlation
Cross.Corr <- function(x1,y,filter.data,date,OS,Frequ){
  
  x1 <- as.character(x1)
  
  if(str_length(x1[1])==7){
  x1 <- paste0(x1,"-01")}
  
  if(OS){
    if(is.matrix(filter.data)){
      y = y[which(x1 >= date),]
    } else {
      y = y[which(x1 >= date)]
    }
  }
  
  if(!is.matrix(filter.data)){
    leng <- length(y)-1
    leng1 <- length(filter.data)
    data <- rep(NA,leng)
    for (j in (leng1+1):(leng)){
      data[j]<-filter.data%*%y[(j):(j-(leng1)+1)]}
  } else {
    leng <- nrow(y)-1
    leng1 <- nrow(filter.data)
    data <- rep(NA,leng)
    for (j in (leng1+1):(leng)){
      data[j]<-sum(apply(filter.data*y[(j):(j-leng1+1),2:ncol(y)],1,sum))}
  }
  
  non_na_vec <- which(!is.na(data))
  first_non_na.1 <- min(non_na_vec)-1
  
  x <- data[-(1:first_non_na.1)]
  y <- y[-(1:first_non_na.1)]
  
  ccf.data <- ccf(Skalieren(x),Skalieren(y),plot = FALSE)
  n <- (length(ccf.data$lag)-1)/2
  p <- plot_ly(x = as.numeric(ccf.data$lag), y = as.numeric(ccf.data$acf),
               type = "bar",marker = list(color = c(rep("black",n),"green",rep("black",n)), name = "Underlying Data",line = list(color = "black"))) %>%
    layout(legend = list(),  xaxis = list(title = 'Lag'), yaxis = list(title = 'Correlation'))
  
  # Best Correlation
  Best <- which(as.numeric(ccf.data$acf) == max(as.numeric(ccf.data$acf),na.rm = T))
  Best_Lagg <- as.numeric(ccf.data$lag)[Best]
  # Text
  text_part <- ifelse(Best_Lagg>=0,"The Filter has a positive Lag of",
                      "The Filter has a negative Lag of")
  text1 <- paste("The Correlation is highest at a Lag of",Best_Lagg, ".",text_part,
                 ifelse(Frequ == "3-Daily",abs(Best_Lagg)*3,ifelse(Frequ == "2-Daily",abs(Best_Lagg)*2,
                                                                   ifelse(Frequ == "4-Daily",abs(Best_Lagg)*4,abs(Best_Lagg)))),
                 "Days",".")
  list(p,text1)
}

### Cross-Correlation aktuelle Resultate
Cross.Corr.akt.R <- function(data,data.1,data.2,data.3,data.4,data.5){
  
  non_na_vec <- which(!is.na(data.1))
  first_non_na.1 <- min(non_na_vec)-1
  non_na_vec <- which(!is.na(data.2))
  first_non_na.2 <- min(non_na_vec)-1
  non_na_vec <- which(!is.na(data.3))
  first_non_na.3 <- min(non_na_vec)-1
  non_na_vec <- which(!is.na(data.4))
  first_non_na.4 <- min(non_na_vec)-1
  non_na_vec <- which(!is.na(data.5))
  first_non_na.5 <- min(non_na_vec)-1
  
  ccf.data.1 <- ccf(Skalieren(data.1[-(1:first_non_na.1)]),Skalieren(data[-(1:first_non_na.1)]),plot = FALSE,na.action = na.pass,lag.max = 260)
  ccf.data.2 <- ccf(Skalieren(data.2[-(1:first_non_na.2)]),Skalieren(data[-(1:first_non_na.2)]),plot = FALSE,na.action = na.pass,lag.max = 260)
  ccf.data.3 <- ccf(Skalieren(data.3[-(1:first_non_na.3)]),Skalieren(data[-(1:first_non_na.3)]),plot = FALSE,na.action = na.pass,lag.max = 260)
  ccf.data.4 <- ccf(Skalieren(data.4[-(1:first_non_na.4)]),Skalieren(data[-(1:first_non_na.4)]),plot = FALSE,na.action = na.pass,lag.max = 260)
  ccf.data.5 <- ccf(Skalieren(data.5[-(1:first_non_na.5)]),Skalieren(data[-(1:first_non_na.5)]),plot = FALSE,na.action = na.pass,lag.max = 260)
  
  Best.1 <- which(as.numeric(ccf.data.1$acf) == max(as.numeric(ccf.data.1$acf),na.rm = T))
  Best.2 <- which(as.numeric(ccf.data.2$acf) == max(as.numeric(ccf.data.2$acf),na.rm = T))
  Best.3 <- which(as.numeric(ccf.data.3$acf) == max(as.numeric(ccf.data.3$acf),na.rm = T))
  Best.4 <- which(as.numeric(ccf.data.4$acf) == max(as.numeric(ccf.data.4$acf),na.rm = T))
  Best.5 <- which(as.numeric(ccf.data.5$acf) == max(as.numeric(ccf.data.5$acf),na.rm = T))
  
  mat <- rbind(c(as.numeric(ccf.data.1$lag)[Best.1],as.numeric(ccf.data.2$lag)[Best.2],as.numeric(ccf.data.3$lag)[Best.3],
                 as.numeric(ccf.data.4$lag)[Best.4],as.numeric(ccf.data.5$lag)[Best.5]),
               c(as.numeric(ccf.data.1$acf)[Best.1],as.numeric(ccf.data.2$acf)[Best.2],as.numeric(ccf.data.3$acf)[Best.3],
                 as.numeric(ccf.data.4$acf)[Best.4],as.numeric(ccf.data.5$acf)[Best.5]))
  
  colnames(mat) <- c("Monatlich","2-W?chentlich","W?chentlich","T?glich","Durchschnitt")
  rownames(mat) <- c("Lag","Correlation")
  mat[1,] <- paste(mat[1,],"Tage") 
  mat[2,] <- round(as.numeric(mat[2,]),2)
  p <- plot_ly(x = as.numeric(ccf.data.1$lag), y = as.numeric(ccf.data.1$acf),
               type = "bar", name = "Monatlich") %>%
    layout(legend = list(),  xaxis = list(title = 'Lag'), yaxis = list(title = '')) %>%
    add_bars(y = as.numeric(ccf.data.2$acf),name = '2-W?chentlich') %>%
    add_bars(y = as.numeric(ccf.data.3$acf),name = 'W?chentlich') %>%
    add_bars(y = as.numeric(ccf.data.4$acf),name = "T?glich") %>%
    add_bars(y = as.numeric(ccf.data.5$acf),name = "Durchschnitt")
  
  list(p,mat)
}

### Performance Matrix
Performance.Table <- function(y,dates,date,OS){
  
  if(!OS){
    Daten <- as.Date(dates$Date)
  } else {
    dates <- subset(dates,as.Date(Date)>=as.Date(date))
    Daten <- as.Date(dates$Date)
  }
  
  y <- subset(y,as.Date(dates) %in% Daten)
  
  if(nrow(y)<=1){
    y[1:2,] <- NA 
  }
  
  Differences <- as.data.frame(diff(y[,1]))
  if(length(Daten) > 0){
    Differences$Trend <- dates$Trend[1:(length(Daten)-1)]
  } else {
    Differences$Trend <- 0
  }
  colnames(Differences) <- c("Difference","Trend")
  
  Hitratio.pos <- nrow(subset(Differences,Trend == "Positive" & as.numeric(Difference)>0))/nrow(subset(Differences,Trend == "Positive"))*100
  Hitratio.neg <- nrow(subset(Differences,Trend == "Negative" & as.numeric(Difference)>0))/nrow(subset(Differences,Trend == "Negative"))*100
  Hitratio.total <- nrow(subset(Differences,as.numeric(Difference)>0))/nrow(Differences)*100
  
  Hitratio.pos <- paste0(round(as.numeric(Hitratio.pos),2),"%")
  Hitratio.neg <- paste0(round(as.numeric(Hitratio.neg),2),"%")
  Hitratio.total <- paste0(round(as.numeric(Hitratio.total),2),"%")
  
  Mat <- as.data.frame(matrix(c(Hitratio.pos,Hitratio.neg,Hitratio.total),nrow = 1))
  colnames(Mat) <- c("Positive Trends","Negative Trends","Total")
  Mat
}

# Performance Matrix Box & Chart with Boxes
Performance.Table.Box <- function(x,y,dates,date,OS,Variance){
  
  if(!OS){
    Daten <- as.Date(dates$Date)
  } else {
    dates <- subset(dates,as.Date(Date)>=as.Date(date))
    Daten <- as.Date(dates$Date)
  }
  y <- as.data.frame(y)
  y$dates <- x
  y1 <- subset(y,as.Date(dates) %in% Daten)
  
  Start <- as.data.frame(y1[,1])
  Start$Trend <- dates$Trend[1:(length(Daten))]
  
  colnames(Start) <- c("Val","Trend")
  
  Thresholds <- as.data.frame(cbind(y1[,1]+Variance,y1[,1]-Variance))
  Thresholds$Trend <- dates$Trend[1:(length(Daten))]
  Thresholds$Date <- Daten
  colnames(Thresholds) <- c("Upper","Lower","Trend","Date")
  Thresholds$Val <- y1[,1]
  Thresholds$Value <- NA
  Thresholds$Compare <- NA
  Thresholds$Compare[which(Thresholds$Trend == "Positive")] <- 1
  Thresholds$Compare[which(Thresholds$Trend == "Negative")] <- -1
  
  for (i in 1:(nrow(Start)-1)) {
    subset.y <- subset(y,as.Date(dates) >= Thresholds$Date[i] & as.Date(dates) <= Thresholds$Date[i+1])
    k <- 1
    while (subset.y[k,1]<Thresholds$Upper[i] & subset.y[k,1]>Thresholds$Lower[i] & !is.na(subset.y[k,1])) {
        k <- k+1
    } 
    Thresholds$Value[i] <- ifelse(is.na(subset.y[k,1]),0,ifelse(subset.y[k,1]<=Thresholds$Lower[i],-1,1))
    
  }
  a <- y$dates[nrow(y)]
  subset.y <- subset(y,as.Date(dates) >= Thresholds$Date[i+1] & as.Date(dates) <= as.Date(a))
  k <- 1
  while (subset.y[k,1]<Thresholds$Upper[i+1] & subset.y[k,1]>Thresholds$Lower[i+1] & !is.na(subset.y[k,1])) {
    k <- k+1
  } 
  Thresholds$Value[i+1] <- ifelse(is.na(subset.y[k,1]),0,ifelse(subset.y[k,1]<=Thresholds$Lower[i+1],-1,1))
  
  Thresholds$Final <-  Thresholds$Compare + Thresholds$Value
  
  Hitratio.pos <- nrow(subset(Thresholds,Trend == "Positive" & as.numeric(Final)==2))/nrow(subset(Thresholds,Trend == "Positive"))*100
  Hitratio.neg <- nrow(subset(Thresholds,Trend == "Negative" & as.numeric(Final)==-2))/nrow(subset(Thresholds,Trend == "Negative"))*100
  Hitratio.pos.n <- nrow(subset(Thresholds,Trend == "Positive" & as.numeric(Final)==1))/nrow(subset(Thresholds,Trend == "Positive"))*100
  Hitratio.neg.n <- nrow(subset(Thresholds,Trend == "Negative" & as.numeric(Final)==-1))/nrow(subset(Thresholds,Trend == "Negative"))*100
  Hitratio.pos.neg <- nrow(subset(Thresholds,Trend == "Positive" & as.numeric(Final)==0))/nrow(subset(Thresholds,Trend == "Positive"))*100
  Hitratio.neg.neg <- nrow(subset(Thresholds,Trend == "Negative" & as.numeric(Final)==0))/nrow(subset(Thresholds,Trend == "Negative"))*100
  Hitratio.total <- nrow(subset(Thresholds,abs(as.numeric(Final))==2))/nrow(Thresholds)*100
  Hitratio.total.n <- nrow(subset(Thresholds,abs(as.numeric(Final))==1))/nrow(Thresholds)*100
  Hitratio.total.neg <- nrow(subset(Thresholds,abs(as.numeric(Final))==0))/nrow(Thresholds)*100
  
  Hitratio.pos <- paste0(round(as.numeric(Hitratio.pos),2),"%")
  Hitratio.neg <- paste0(round(as.numeric(Hitratio.neg),2),"%")
  Hitratio.total <- paste0(round(as.numeric(Hitratio.total),2),"%")
  Hitratio.pos.n <- paste0(round(as.numeric(Hitratio.pos.n),2),"%")
  Hitratio.neg.n <- paste0(round(as.numeric(Hitratio.neg.n),2),"%")
  Hitratio.total.n <- paste0(round(as.numeric(Hitratio.total.n),2),"%")
  Hitratio.pos.neg <- paste0(round(as.numeric(Hitratio.pos.neg),2),"%")
  Hitratio.neg.neg <- paste0(round(as.numeric(Hitratio.neg.neg),2),"%")
  Hitratio.total.neg <- paste0(round(as.numeric(Hitratio.total.neg),2),"%")
  
  Mat <- as.data.frame(matrix(c(Hitratio.pos,Hitratio.pos.n,Hitratio.pos.neg),ncol = 1,nrow = 3))
  Mat$`Negative Trends` <- c(Hitratio.neg,Hitratio.neg.n,Hitratio.neg.neg)
  Mat$Total <- c(Hitratio.total,Hitratio.total.n,Hitratio.total.neg)
  colnames(Mat) <- c("Positive Trends","Negative Trends","Total")
  rownames(Mat) <- c("Korrekter Trend","Kein Trend","Entgegengesetzter Trend")
  
  line1 <- list(
    type = "rect",
    fillcolor = "grey",
    line = list(color = "grey", width = 0.1),
    opacity = 0.3,
    xref = "x",
    yref = "y",
    name = "Performance Box"
  )
  
  lines1 <- list()
  for (i in 1:(nrow(Thresholds)-1)) {
    line1[["x0"]] <- Thresholds$Date[i]
    line1[["x1"]] <- Thresholds$Date[i+1]
    line1[c("y0", "y1")] <- c(Thresholds$Lower[i], Thresholds$Upper[i])
    lines1 <- c(lines1, list(line1))
  }
  line1[["x0"]] <- Thresholds$Date[i+1]
  line1[["x1"]] <- as.Date(a)
  line1[c("y0", "y1")] <- c(Thresholds$Lower[i+1], Thresholds$Upper[i+1])
  lines1 <- c(lines1, list(line1))
  
  list(Mat,Thresholds,lines1)
}

### Result Matrix
Result <- function(data,data.1,data.2,data.3,data.4,data.5,x){
  
  CCor <- Cross.Corr.akt.R(data,data.1,data.2,data.3,data.4,data.5)[[2]]
  row1 <- c(ifelse(data.1[length(data.1)]>0,"Steigend","Fallend"),
            ifelse(data.2[length(data.2)]>0,"Steigend","Fallend"),
            ifelse(data.3[length(data.3)]>0,"Steigend","Fallend"),
            ifelse(data.4[length(data.4)]>0,"Steigend","Fallend"),
            ifelse(data.5[length(data.5)]>0,"Steigend","Fallend"))
  
  res <- c(NA,data.1)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  indx <- which(res1!=0)
  indx.1 <- as.Date(x)[indx[length(indx)]]
  
  res <- c(NA,data.2)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  indx <- which(res1!=0)
  indx.2 <- as.Date(x)[indx[length(indx)]]
  
  res <- c(NA,data.3)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  indx <- which(res1!=0)
  indx.3 <- as.Date(x)[indx[length(indx)]]
  
  res <- c(NA,data.4)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  indx <- which(res1!=0)
  indx.4 <- as.Date(x)[indx[length(indx)]]
  
  res <- c(NA,data.5)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  indx <- which(res1!=0)
  indx.5 <- as.Date(x)[indx[length(indx)]]
  
  row2 <- paste0(c(indx.1,indx.2,indx.3,indx.4,indx.5)," (",as.Date(x)[length(as.Date(x))]-c(indx.1,indx.2,indx.3,indx.4,indx.5)," Tage)")
  mat <- rbind(row1,row2,CCor[1,])
  colnames(mat) <- c("Monatlich","2-W?chentlich","W?chentlich","T?glich","Durchschnitt")
  rownames(mat) <- c("Aktueller Trend","Trend l?uft seit","Durchschnittlicher Lag")
  
  mat
}

### Zinskurve
Zinskurve <- function(x.3,x.12,x.60,x.120,x.240,x.360,y.3,y.12,y.60,y.120,y.240,y.360,
                      double = FALSE,y2.3 = 0,y2.12 = 0,y2.60 = 0,y2.120 = 0,y2.240 = 0,y2.360 = 0,Variate = ""){
  x <- c(3,12,60,120,240,360)
  
  y.3 <- ifelse(y.3>=0,x.3+0.1,x.3-0.1)
  y.12 <- ifelse(y.12>=0,x.12+0.1,x.12-0.1)
  y.60 <- ifelse(y.60>=0,x.60+0.1,x.60-0.1)
  y.120 <- ifelse(y.120>=0,x.120+0.1,x.120-0.1)
  y.240 <- ifelse(y.240>=0,x.240+0.1,x.240-0.1)
  y.360 <- ifelse(y.360>=0,x.360+0.1,x.360-0.1)
  
  if(double){
    y2.3 <- ifelse(y2.3>=0,x.3+0.1,x.3-0.1)
    y2.12 <- ifelse(y2.12>=0,x.12+0.1,x.12-0.1)
    y2.60 <- ifelse(y2.60>=0,x.60+0.1,x.60-0.1)
    y2.120 <- ifelse(y2.120>=0,x.120+0.1,x.120-0.1)
    y2.240 <- ifelse(y2.240>=0,x.240+0.1,x.240-0.1)
    y2.360 <- ifelse(y2.360>=0,x.360+0.1,x.360-0.1)
  }
  
  p <- plot_ly(x = x, y = c(x.3,x.12,x.60,x.120,x.240,x.360),
          type = 'scatter', mode = 'lines', name = "Aktuelle Zinskurve", line = list(shape = "spline")) %>% 
    layout(xaxis = list(title = "Duration"),yaxis = list(title = "Yield")) %>% 
    add_trace(y = c(y.3,y.12,y.60,y.120,y.240,y.360), name = paste("Forecast",Variate),line = list(shape = "spline"),
              mode = 'lines',opacity = 0.8)
  
  if(double){
    p <- p %>%
      add_trace(y = c(y2.3,y2.12,y2.60,y2.120,y2.240,y2.360), name = "Forecast Multivariate",line = list(shape = "spline"),
                mode = 'lines',opacity = 0.8)
  }
  
  p
  
}


# Zinskurve & Resultat Table
Zinskurve.Resultate.Table.DFA <- function(data,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                      L2,k2,lambda2,eta2,cutoff2,Constraint2,Frequence,Table.Bool,Explanatory,Double,L.Vect,k.Vect){
  
  # Durations
  Durations <- c(3,12,60,120,240,360)
  L.Vect <- as.numeric(L.Vect)
  k.Vect <- as.numeric(k.Vect)
  # Data Preparation
  data.Dur <- vector(mode = "list", length = 6)
  data.DurD <- vector(mode = "list", length = 6)
  data.DurW <- vector(mode = "list", length = 6)
  data.Dur2W <- vector(mode = "list", length = 6)
  data.DurM <- vector(mode = "list", length = 6)
  
  CC.Mat <- vector(mode = "list", length = 6)
  CC.MatD <- vector(mode = "list", length = 6)
  CC.MatW <- vector(mode = "list", length = 6)
  CC.Mat2W <- vector(mode = "list", length = 6)
  CC.MatM <- vector(mode = "list", length = 6)
  
  Zinskurve.current <- rep(NA,6)
  
  if(Frequence != "Durchschnitt"){
    # DFA on each Duration
    for(i in 1:6){
      data.Interim <- Choice(data,Explanatory,Durations[i],"Multivariate",Frequence,"1900-01-01",Sys.Date())
      Zinskurve.current[i] <- data.Interim[[1]][nrow(data.Interim[[1]]),2]
      data.single <- Univariate.Filter.cust(data.Interim[[2]][,1],0,L1,k1,lambda1,eta1,cutoff1,Constraint1)[[1]]
      
      data.Dur[[i]] <- if(!Double){data.single} else {
        Univariate.Filter.cust.conv(data.Interim[[2]][,1],0,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                    0,L2,k2,lambda2,eta2,cutoff2,Constraint2)[[1]]}
      
      if(Table.Bool){ 
        CC.Mat[[i]] <- Result.2.TS(data.Interim[[2]][,1],data.Dur[[i]],data.Interim[[1]][,1],Frequence)
      }
      
    }
  } else {
    
    # DFA on each Duration
    for(i in 1:6){
      # Data Preparation
      data.Freq <- vector(mode = "list", length = 4)
      date.Freq <- vector(mode = "list", length = 4)
      jj <- 1
      # DFA on each Frequence
      for(j in c("T?glich","W?chentlich","2-W?chentlich","Monatlich")){
        data.Interim <- Choice(data,Explanatory,Durations[i],"Multivariate",j,"1900-01-01",Sys.Date())
        
        if(j == "T?glich"){
          data.Interim.T <- data.Interim
        }
        
        Zinskurve.current[i] <- data.Interim[[1]][nrow(data.Interim[[1]]),2]
        data.single <- Univariate.Filter.cust(data.Interim[[2]][,1],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1)[[1]]
        date.Freq[[jj]] <- as.character(data.Interim[[1]][-1,1])
        data.Freq[[jj]] <- if(!Double){data.single} else {
          Univariate.Filter.cust.conv(data.Interim[[2]][,1],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1,
                                      0,L2,1/k.Vect[jj],lambda2,eta2,pi/k.Vect[jj],Constraint2)[[1]]}
        jj <- jj+1
      }
      
      data.Mean <- Mean.of.Vectors(date.Freq[[1]],date.Freq[[2]],date.Freq[[3]],date.Freq[[4]],Skalieren(data.Freq[[1]]),Skalieren(data.Freq[[2]]),Skalieren(data.Freq[[3]]),Skalieren(data.Freq[[4]]))
      data.Dur[[i]] <- data.Mean[[1]]
      data.DurD[[i]] <- data.Mean[[5]]
      data.DurW[[i]] <- data.Mean[[4]]
      data.Dur2W[[i]] <- data.Mean[[3]]
      data.DurM[[i]] <- data.Mean[[2]]
      if(Table.Bool){ 
        CC.Mat[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Mean[[1]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatD[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurD[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatW[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurW[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2W[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur2W[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatM[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurM[[i]],data.Interim.T[[1]][,1],"T?glich")
      }
      
    }
  }
  
  # Zinskurve
  p <- Zinskurve(Zinskurve.current[1],Zinskurve.current[2],Zinskurve.current[3],
                 Zinskurve.current[4],Zinskurve.current[5],Zinskurve.current[6],
                 data.Dur[[1]][length(data.Dur[[1]])],data.Dur[[2]][length(data.Dur[[2]])],
                 data.Dur[[3]][length(data.Dur[[3]])],data.Dur[[4]][length(data.Dur[[4]])],
                 data.Dur[[5]][length(data.Dur[[5]])],data.Dur[[6]][length(data.Dur[[6]])])
  
  mat <- matrix(NA,nrow = 5,ncol = 6)
  
  if(Table.Bool){
    
    mat <- Table.End(CC.Mat[[1]],CC.Mat[[2]],CC.Mat[[3]],CC.Mat[[4]],CC.Mat[[5]],CC.Mat[[6]])
    
    if(Frequence == "Durchschnitt"){
      matD <- Table.End(CC.MatD[[1]],CC.MatD[[2]],CC.MatD[[3]],CC.MatD[[4]],CC.MatD[[5]],CC.MatD[[6]])
      
      matW <- Table.End(CC.MatW[[1]],CC.MatW[[2]],CC.MatW[[3]],CC.MatW[[4]],CC.MatW[[5]],CC.MatW[[6]])
      
      mat2W <- Table.End(CC.Mat2W[[1]],CC.Mat2W[[2]],CC.Mat2W[[3]],CC.Mat2W[[4]],CC.Mat2W[[5]],CC.Mat2W[[6]])
      
      matM <- Table.End(CC.MatM[[1]],CC.MatM[[2]],CC.MatM[[3]],CC.MatM[[4]],CC.MatM[[5]],CC.MatM[[6]])
      
      
    } else {
      
      matD <- data.frame(matrix(0))
      
      matW <- data.frame(matrix(0))
      
      mat2W <- data.frame(matrix(0))
      
      matM <- data.frame(matrix(0))
      
    }
    
  }
  
  if(Frequence != "Durchschnitt"){
    a.list <- list(mat)
  } else {
    a.list <- list(matD,matW,mat2W,matM,mat)
  }
  
  
  list(p,mat,mat,a.list)
  
}

# Zinskurve & Resultat Table
Zinskurve.Resultate.Table.Beide <- function(data,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                          L2,k2,lambda2,eta2,cutoff2,Constraint2,Frequence,Table.Bool,Explanatory,Double,L.Vect,k.Vect){
  
  # Durations
  Durations <- c(3,12,60,120,240,360)
  L.Vect <- as.numeric(L.Vect)
  k.Vect <- as.numeric(k.Vect)
  # Data Preparation
  data.Dur <- vector(mode = "list", length = 6)
  data.DurD <- vector(mode = "list", length = 6)
  data.DurW <- vector(mode = "list", length = 6)
  data.Dur2W <- vector(mode = "list", length = 6)
  data.DurM <- vector(mode = "list", length = 6)
  data.Dur2 <- vector(mode = "list", length = 6)
  data.Dur2D <- vector(mode = "list", length = 6)
  data.Dur2W2 <- vector(mode = "list", length = 6)
  data.Dur22W <- vector(mode = "list", length = 6)
  data.Dur2M <- vector(mode = "list", length = 6)
  CC.Mat <- vector(mode = "list", length = 6)
  CC.MatD <- vector(mode = "list", length = 6)
  CC.MatW <- vector(mode = "list", length = 6)
  CC.Mat2W <- vector(mode = "list", length = 6)
  CC.MatM <- vector(mode = "list", length = 6)
  CC.Mat2 <- vector(mode = "list", length = 6)
  CC.Mat2D <- vector(mode = "list", length = 6)
  CC.Mat2W2 <- vector(mode = "list", length = 6)
  CC.Mat22W <- vector(mode = "list", length = 6)
  CC.Mat2M <- vector(mode = "list", length = 6)
  Zinskurve.current <- rep(NA,6)
  
  if(Frequence != "Durchschnitt"){
    # DFA on each Duration
    for(i in 1:6){
      data.Interim <- Choice(data,Explanatory,Durations[i],"Multivariate",Frequence,"1900-01-01",Sys.Date())
      Zinskurve.current[i] <- data.Interim[[1]][nrow(data.Interim[[1]]),2]
      data.single <- Multivariate.Filter.cust(data.Interim[[2]],0,L1,k1,lambda1,eta1,cutoff1,Constraint1)[[1]]
      data.single1 <- Univariate.Filter.cust(data.Interim[[2]][,1],0,L1,k1,lambda1,eta1,cutoff1,Constraint1)[[1]]
      
      data.Dur[[i]] <- if(!Double){data.single} else {
        Multivariate.Filter.cust.conv(data.Interim[[2]],0,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                    0,L2,k2,lambda2,eta2,cutoff2,Constraint2)[[1]]}
      data.Dur2[[i]] <- if(!Double){data.single1} else {
        Univariate.Filter.cust.conv(data.Interim[[2]][,1],0,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                      0,L2,k2,lambda2,eta2,cutoff2,Constraint2)[[1]]}
      
      if(Table.Bool){ 
        CC.Mat[[i]] <- Result.2.TS(data.Interim[[2]][,1],data.Dur[[i]],data.Interim[[1]][,1],Frequence)
        CC.Mat2[[i]] <- Result.2.TS(data.Interim[[2]][,1],data.Dur2[[i]],data.Interim[[1]][,1],Frequence)
      }
      
    }
  } else {
    
    # DFA on each Duration
    for(i in 1:6){
      # Data Preparation
      data.Freq <- vector(mode = "list", length = 4)
      data.Freq2 <- vector(mode = "list", length = 4)
      date.Freq <- vector(mode = "list", length = 4)
      jj <- 1
      # DFA on each Frequence
      for(j in c("T?glich","W?chentlich","2-W?chentlich","Monatlich")){
        data.Interim <- Choice(data,Explanatory,Durations[i],"Multivariate",j,"1900-01-01",Sys.Date())
        
        if(j == "T?glich"){
          data.Interim.T <- data.Interim
        }
        
        Zinskurve.current[i] <- data.Interim[[1]][nrow(data.Interim[[1]]),2]
        data.single <- Multivariate.Filter.cust(data.Interim[[2]],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1)[[1]]
        data.single1 <- Univariate.Filter.cust(data.Interim[[2]][,1],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1)[[1]]
        date.Freq[[jj]] <- as.character(data.Interim[[1]][-1,1])
        data.Freq[[jj]] <- if(!Double){data.single} else {
          Multivariate.Filter.cust.conv(data.Interim[[2]],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1,
                                      0,L2,1/k.Vect[jj],lambda2,eta2,pi/k.Vect[jj],Constraint2)[[1]]}
        data.Freq2[[jj]] <- if(!Double){data.single1} else {
          Univariate.Filter.cust.conv(data.Interim[[2]][,1],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1,
                                        0,L2,1/k.Vect[jj],lambda2,eta2,pi/k.Vect[jj],Constraint2)[[1]]}
        jj <- jj+1
      }
      
      data.Mean <- Mean.of.Vectors(date.Freq[[1]],date.Freq[[2]],date.Freq[[3]],date.Freq[[4]],Skalieren(data.Freq[[1]]),Skalieren(data.Freq[[2]]),Skalieren(data.Freq[[3]]),Skalieren(data.Freq[[4]]))
      data.Mean2 <- Mean.of.Vectors(date.Freq[[1]],date.Freq[[2]],date.Freq[[3]],date.Freq[[4]],Skalieren(data.Freq2[[1]]),Skalieren(data.Freq2[[2]]),Skalieren(data.Freq2[[3]]),Skalieren(data.Freq2[[4]]))
      data.Dur[[i]] <- data.Mean[[1]]
      data.DurD[[i]] <- data.Mean[[5]]
      data.DurW[[i]] <- data.Mean[[4]]
      data.Dur2W[[i]] <- data.Mean[[3]]
      data.DurM[[i]] <- data.Mean[[2]]
      data.Dur2[[i]] <- data.Mean2[[1]]
      data.Dur2D[[i]] <- data.Mean2[[5]]
      data.Dur2W2[[i]] <- data.Mean2[[4]]
      data.Dur22W[[i]] <- data.Mean2[[3]]
      data.Dur2M[[i]] <- data.Mean2[[2]]
      
      if(Table.Bool){ 
        CC.Mat[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Mean[[1]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatD[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurD[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatW[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurW[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2W[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur2W[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatM[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurM[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Mean2[[1]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2D[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur2D[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2W2[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur2W2[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat22W[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur22W[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2M[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur2M[[i]],data.Interim.T[[1]][,1],"T?glich")
      }
      
    }
  }
  
  # Zinskurve
  p <- Zinskurve(Zinskurve.current[1],Zinskurve.current[2],Zinskurve.current[3],
                 Zinskurve.current[4],Zinskurve.current[5],Zinskurve.current[6],
                 data.Dur2[[1]][length(data.Dur2[[1]])],data.Dur2[[2]][length(data.Dur2[[2]])],
                 data.Dur2[[3]][length(data.Dur2[[3]])],data.Dur2[[4]][length(data.Dur2[[4]])],
                 data.Dur2[[5]][length(data.Dur2[[5]])],data.Dur2[[6]][length(data.Dur2[[6]])],double = TRUE,
                 data.Dur[[1]][length(data.Dur[[1]])],data.Dur[[2]][length(data.Dur[[2]])],
                 data.Dur[[3]][length(data.Dur[[3]])],data.Dur[[4]][length(data.Dur[[4]])],
                 data.Dur[[5]][length(data.Dur[[5]])],data.Dur[[6]][length(data.Dur[[6]])],Variate = "Univariate")
  
  mat <- matrix(NA,nrow = 5,ncol = 6)
  
  if(Table.Bool){
    
    mat <- Table.End(CC.Mat[[1]],CC.Mat[[2]],CC.Mat[[3]],CC.Mat[[4]],CC.Mat[[5]],CC.Mat[[6]])
    
    if(Frequence == "Durchschnitt"){
    matD <- Table.End(CC.MatD[[1]],CC.MatD[[2]],CC.MatD[[3]],CC.MatD[[4]],CC.MatD[[5]],CC.MatD[[6]])
    
    matW <- Table.End(CC.MatW[[1]],CC.MatW[[2]],CC.MatW[[3]],CC.MatW[[4]],CC.MatW[[5]],CC.MatW[[6]])
    
    mat2W <- Table.End(CC.Mat2W[[1]],CC.Mat2W[[2]],CC.Mat2W[[3]],CC.Mat2W[[4]],CC.Mat2W[[5]],CC.Mat2W[[6]])
    
    matM <- Table.End(CC.MatM[[1]],CC.MatM[[2]],CC.MatM[[3]],CC.MatM[[4]],CC.MatM[[5]],CC.MatM[[6]])
    
    mat2D <- Table.End(CC.Mat2D[[1]],CC.Mat2D[[2]],CC.Mat2D[[3]],CC.Mat2D[[4]],CC.Mat2D[[5]],CC.Mat2D[[6]])
    
    mat2W2 <- Table.End(CC.Mat2W2[[1]],CC.Mat2W2[[2]],CC.Mat2W2[[3]],CC.Mat2W2[[4]],CC.Mat2W2[[5]],CC.Mat2W2[[6]])
    
    mat22W <- Table.End(CC.Mat22W[[1]],CC.Mat22W[[2]],CC.Mat22W[[3]],CC.Mat22W[[4]],CC.Mat22W[[5]],CC.Mat22W[[6]])
    
    mat2M <- Table.End(CC.Mat2M[[1]],CC.Mat2M[[2]],CC.Mat2M[[3]],CC.Mat2M[[4]],CC.Mat2M[[5]],CC.Mat2M[[6]])
    
    } else {
      
      matD <- data.frame(matrix(0))
      
      matW <- data.frame(matrix(0))
      
      mat2W <- data.frame(matrix(0))
      
      matM <- data.frame(matrix(0))
      
      mat2D <- data.frame(matrix(0))
      
      mat2W2 <- data.frame(matrix(0))
      
      mat22W <- data.frame(matrix(0))
      
      mat2M <- data.frame(matrix(0)) 
      
    }
    mat2 <- Table.End(CC.Mat2[[1]],CC.Mat2[[2]],CC.Mat2[[3]],CC.Mat2[[4]],CC.Mat2[[5]],CC.Mat2[[6]])
    
    

  }
  
  list(p,mat2,mat,list(mat2D,mat2W2,mat22W,mat2M,mat2,matD,matW,mat2W,matM,mat))
  
}

# Build Table
Table.End <- function(a,b,c,d,e,f){
  mat <- cbind(a,b,c,d,e,f)
  mat <- as.data.frame(mat)
  colnames(mat) <- c("3 Monate","12 Monate","60 Monate","120 Monate","240 Monate","360 Monate")
  rownames(mat) <- c("Aktueller Trend","Der Trend l?uft seit","Durchschnittlicher Lag")
  mat
}


# Zinskurve & Resultat Table
Zinskurve.Resultate.Table.MDFA <- function(data,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                           L2,k2,lambda2,eta2,cutoff2,Constraint2,Frequence,Table.Bool,Explanatory,Double,L.Vect,k.Vect){
  
  # Durations
  Durations <- c(3,12,60,120,240,360)
  L.Vect <- as.numeric(L.Vect)
  k.Vect <- as.numeric(k.Vect)
  # Data Preparation
  data.Dur <- vector(mode = "list", length = 6)
  data.DurD <- vector(mode = "list", length = 6)
  data.DurW <- vector(mode = "list", length = 6)
  data.Dur2W <- vector(mode = "list", length = 6)
  data.DurM <- vector(mode = "list", length = 6)
  
  CC.Mat <- vector(mode = "list", length = 6)
  CC.MatD <- vector(mode = "list", length = 6)
  CC.MatW <- vector(mode = "list", length = 6)
  CC.Mat2W <- vector(mode = "list", length = 6)
  CC.MatM <- vector(mode = "list", length = 6)
  Zinskurve.current <- rep(NA,6)
  
  if(Frequence != "Durchschnitt"){
    # DFA on each Duration
    for(i in 1:6){
      data.Interim <- Choice(data,Explanatory,Durations[i],"Multivariate",Frequence,"1900-01-01",Sys.Date())
      Zinskurve.current[i] <- data.Interim[[1]][nrow(data.Interim[[1]]),2]
      data.single <- Multivariate.Filter.cust(data.Interim[[2]],0,L1,k1,lambda1,eta1,cutoff1,Constraint1)[[1]]
      
      data.Dur[[i]] <- if(!Double){data.single} else {
        Multivariate.Filter.cust.conv(data.Interim[[2]],0,L1,k1,lambda1,eta1,cutoff1,Constraint1,
                                      0,L2,k2,lambda2,eta2,cutoff2,Constraint2)[[1]]}
      
      if(Table.Bool){ 
        CC.Mat[[i]] <- Result.2.TS(data.Interim[[2]][,1],data.Dur[[i]],data.Interim[[1]][,1],Frequence)
      }
      
    }
  } else {
    
    # DFA on each Duration
    for(i in 1:6){
      # Data Preparation
      data.Freq <- vector(mode = "list", length = 4)
      date.Freq <- vector(mode = "list", length = 4)
      jj <- 1
      # DFA on each Frequence
      for(j in c("T?glich","W?chentlich","2-W?chentlich","Monatlich")){
        data.Interim <- Choice(data,Explanatory,Durations[i],"Multivariate",j,"1900-01-01",Sys.Date())
        
        if(j == "T?glich"){
          data.Interim.T <- data.Interim
        }
        
        Zinskurve.current[i] <- data.Interim[[1]][nrow(data.Interim[[1]]),2]
        data.single <- Multivariate.Filter.cust(data.Interim[[2]],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1)[[1]]
        date.Freq[[jj]] <- as.character(data.Interim[[1]][-1,1])
        data.Freq[[jj]] <- if(!Double){data.single} else {
          Multivariate.Filter.cust.conv(data.Interim[[2]],0,L.Vect[jj],1/k.Vect[jj],lambda1,eta1,pi/k.Vect[jj],Constraint1,
                                        0,L2,1/k.Vect[jj],lambda2,eta2,pi/k.Vect[jj],Constraint2)[[1]]}
        jj <- jj+1
      }
      
      data.Mean <- Mean.of.Vectors(date.Freq[[1]],date.Freq[[2]],date.Freq[[3]],date.Freq[[4]],Skalieren(data.Freq[[1]]),Skalieren(data.Freq[[2]]),Skalieren(data.Freq[[3]]),Skalieren(data.Freq[[4]]))
      data.Dur[[i]] <- data.Mean[[1]]
      data.Dur[[i]] <- data.Mean[[1]]
      data.DurD[[i]] <- data.Mean[[5]]
      data.DurW[[i]] <- data.Mean[[4]]
      data.Dur2W[[i]] <- data.Mean[[3]]
      data.DurM[[i]] <- data.Mean[[2]]
      if(Table.Bool){ 
        CC.Mat[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Mean[[1]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatD[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurD[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatW[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurW[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.Mat2W[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.Dur2W[[i]],data.Interim.T[[1]][,1],"T?glich")
        CC.MatM[[i]] <- Result.2.TS(data.Interim.T[[2]][,1],data.DurM[[i]],data.Interim.T[[1]][,1],"T?glich")
      }
      
    }
  }
  
  # Zinskurve
  p <- Zinskurve(Zinskurve.current[1],Zinskurve.current[2],Zinskurve.current[3],
                 Zinskurve.current[4],Zinskurve.current[5],Zinskurve.current[6],
                 data.Dur[[1]][length(data.Dur[[1]])],data.Dur[[2]][length(data.Dur[[2]])],
                 data.Dur[[3]][length(data.Dur[[3]])],data.Dur[[4]][length(data.Dur[[4]])],
                 data.Dur[[5]][length(data.Dur[[5]])],data.Dur[[6]][length(data.Dur[[6]])])
  
  mat <- matrix(NA,nrow = 5,ncol = 6)
  
  if(Table.Bool){
    
    mat <- Table.End(CC.Mat[[1]],CC.Mat[[2]],CC.Mat[[3]],CC.Mat[[4]],CC.Mat[[5]],CC.Mat[[6]])
    
    if(Frequence == "Durchschnitt"){
      matD <- Table.End(CC.MatD[[1]],CC.MatD[[2]],CC.MatD[[3]],CC.MatD[[4]],CC.MatD[[5]],CC.MatD[[6]])
      
      matW <- Table.End(CC.MatW[[1]],CC.MatW[[2]],CC.MatW[[3]],CC.MatW[[4]],CC.MatW[[5]],CC.MatW[[6]])
      
      mat2W <- Table.End(CC.Mat2W[[1]],CC.Mat2W[[2]],CC.Mat2W[[3]],CC.Mat2W[[4]],CC.Mat2W[[5]],CC.Mat2W[[6]])
      
      matM <- Table.End(CC.MatM[[1]],CC.MatM[[2]],CC.MatM[[3]],CC.MatM[[4]],CC.MatM[[5]],CC.MatM[[6]])
      
      
    } else {
      
      matD <- data.frame(matrix(0))
      
      matW <- data.frame(matrix(0))
      
      mat2W <- data.frame(matrix(0))
      
      matM <- data.frame(matrix(0))
      
    }
  }
  
  if(Frequence != "Durchschnitt"){
    a.list <- list(mat)
  } else {
    a.list <- list(matD,matW,mat2W,matM,mat)
  }
  
  
  list(p,mat,mat,a.list)
}

# Result for two Time-Series
Result.2.TS <- function(data,data.1,x,Frequence){
  
  Frequenz <- ifelse(Frequence == "T?glich","Tage",
                     ifelse(Frequence == "W?chentlich","Wochen",
                             ifelse(Frequence == "2-W?chentlich","Wochen","Monate")))
  
  non_na_vec <- which(!is.na(data.1))
  first_non_na.1 <- min(non_na_vec)-1
  ccf.data.1 <- ccf(data.1[-(1:first_non_na.1)],data[-(1:first_non_na.1)],plot = FALSE,na.action = na.pass,lag.max = 260)
  Best.1 <- which(as.numeric(ccf.data.1$acf) == max(as.numeric(ccf.data.1$acf),na.rm = T))
  CCor <- rbind(as.numeric(ccf.data.1$lag)[Best.1],
               as.numeric(ccf.data.1$acf)[Best.1])
  if(Frequence == "2-W?chentlich"){
    CCor[1,] <- 2*CCor[1,]
  }
  
  colnames(CCor) <- c("Frequence")
  rownames(CCor) <- c("Lag","Correlation")
  
  row1 <- ifelse(data.1[length(data.1)]>0,"Steigend","Fallend")
  
  res <- c(NA,data.1)
  res[res == 0] <- NA
  res_min <- c(NA,res)[1:(length(res))]
  res1 <- ifelse(res>res_min & res>0 & res_min < 0,1 ,
                 ifelse(res_min>res & res<0 & res_min > 0,-1,0))
  indx <- which(res1!=0)
  indx.1 <- as.Date(x)[indx[length(indx)]]
  
  row2 <- paste0(indx.1," (",as.Date(x)[length(as.Date(x))]-indx.1," Tage)")
  mat <- rbind(row1,row2,paste(CCor[1,],Frequenz))
  colnames(mat) <- c("Frequence")
  rownames(mat) <- c("Aktueller Trend","Trend l?uft seit","Durchschnittlicher Lag")
  
  mat
}

# Trendwenden
Trendwenden.aktuell <- function(mat.list,date,Variate){
  
  n <- length(mat.list)
  
  mat.final <- as.data.frame(matrix(NA,1,6))
  colnames(mat.final) <- c("Approach","Daten-Frequenz","Trend","Trend l?uft seit","Lag","Duration")
  Frequenz <- rep(c("T?glich","W?chentlich","2-W?chentlich","Monatlich","Durchschnitt"),2)
  for (i in 1:n) {
    
    Approach <- ifelse(i > 5 | Variate == "Multivariate","Multivariate","Univariate")
    
    columns.selected <- apply(mat.list[[i]][2,], 2, function(r) any(substr(r,1,10) >= date))
    c.names <- colnames(mat.list[[i]])[columns.selected]
    mat <- cbind(t(as.data.frame(mat.list[[i]][,columns.selected])),c.names)
    k <- nrow(mat)
    mat <- cbind(rep(Approach,k),rep(Frequenz[i],k),mat)
    colnames(mat) <- c("Approach","Daten-Frequenz","Trend","Trend l?uft seit","Lag","Duration")
    mat.final <- rbind(mat.final,mat)
  }
  mat1 <- mat.final[-1,]
  mat1 <- as.data.frame(mat1)
  colnames(mat1) <- c("Approach","Daten-Frequenz","Trend","Trend l?uft seit","Lag","Duration")
  mat1[order(as.numeric(gsub( " .*$", "",mat1$Duration))),]
}
