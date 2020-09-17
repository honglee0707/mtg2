subroutine trf_h2 (h2,se,k,p11,indc)

!**********************************************************************
! Expected prediction accuracy (also see ccgpa5.r in /R/Daetwyler/)
! S. Hong Lee  (2016)
!************************************************************************

implicit none

  real ( kind = 8 ) bound,mean,p,q,sd,x
  integer ( kind = 4 ) status,which

  integer::N,M
  real::h2,se,k,p11,p12,ckv,iv,iv2,thd,zv,cv,theta,lamda,vgcc,Q_r2,r2
  real::cv2,theta2,h2o,h2o2,v102,v103,qv,cthd,czv,civ,p_case_s_int
  real::var_l_s,p_case_s,p_case_s11,var_l_s_,p_case_s_int_,p_case_s11_
  real::p_case_s12,p_case_s12_,estim_OR2,estim_OR3
  character(len=64)::indc

  if (indc=="obs") then
   print*,"*************************************************************"
   print*,"h2     : estimated h2 on the observed scale using LMM"
   print*,"se     : its standard error"
   print*,"k      : population prevalence"
   print*,"p      : proportion of cases in training sample"
   print*,"method : liability transformation (AJHG 88: 294-305)"
   print*,"*************************************************************"
   print*,""
  elseif (indc=="liab") then
   print*,"*************************************************************"
   print*,"h2     : estimated h2 on the liability scale"
   print*,"se     : its standard error"
   print*,"k      : population prevalence"
   print*,"p      : proportion of cases in training sample"
   print*,"method : liability transformation (AJHG 88: 294-305)"
   print*,"*************************************************************"
   print*,""

  elseif (indc=="R2CS") then
   print*,"*************************************************************"
   print*,"R2     : Cox & Snell R2 (equivalent to R2 on the observed scale)"
   print*,"se     : its standard error"
   print*,"k      : population prevalence"
   print*,"p      : proportion of cases in training sample"
   print*,"method : liability transformation (Genet Epidemiol 36: 214-224)"
   print*,"*************************************************************"
   print*,""
  elseif (indc=="R2N") then
   print*,"*************************************************************"
   print*,"R2     : Nagelkerke R2"
   print*,"se     : its standard error"
   print*,"k      : population prevalence"
   print*,"p      : proportion of cases in training sample"
   print*,"method : liability transformation (Genet Epidemiol 36: 214-224)"
   print*,"*************************************************************"
   print*,""

  else
   print*,"please indicate what scale the heritablity is >>> 'liab' or 'obs'"
  end if

!CDF functions *****************************************************
      p = 1D+00 - k !0.011135D+00
      q = 1.0D+00 - p
      !p = huge ( p )
      !q = huge ( q )
      x = huge ( x )
      !x=-1.28155 !2 3.0D+00
      mean = 0.0D+00 ! 5.0D+00
      sd = 1.0D+00 !0.875D+00
!which=1  !Calculate P and Q from X, MEAN and SD;
which=2  !Calculate X from P, Q, MEAN and SD;
!which=3  !Calculate MEAN from P, Q, X and SD;
!which=4  !Calculate SD from P, Q, X and MEAN;

call cdfnor ( which, p, q, x, mean, sd, status, bound )
!print*,p,q,x,mean,sd,status,bound
thd=x
zv=1/(sd*sqrt(2*3.141593)) * exp(-0.5*((x-mean)/sd)**2)
!print*,thd,zv

   iv=zv/k          !#mean liability for cases
   iv2=-iv*k/(1-k)  !#mean liability for controls

   cv=(k*(1-k))**2/(zv**2*p11*(1-p11))    !#the spread sheet
   theta=iv*((p11-k)/(1-k))*(iv*((p11-k)/(1-k))-thd)

 if (indc=="obs") then   
  h2=h2*cv
  se=se*cv
  print*,'Transformed h2 on the liability scale    : ',h2
  print*,'Transformed SE on the liability scale    : ',se
  print*,""

 elseif (indc=="liab") then
  h2=h2/cv
  se=se/cv
  print*,'Transformed h2 on the observed scale    : ',h2
  print*,'Transformed SE on the observed scale    : ',se
  print*,""

 elseif (indc=="R2CS") then
  
  se=se*cv/(cv*theta*h2+1)**2  !Taylor expansion
  !se=se*(h2*cv/(1+h2*theta*cv))/h2
  h2=h2*cv/(1+h2*theta*cv)
  print*,'Transformed R2 on the liabiltiy scale    : ',h2
  print*,'Transformed SE on the liability scale    : ',se
  print*,""

 elseif (indc=="R2N") then
  se=se*(1-p11**(2*p11)*(1-p11)**(2*(1-p11)))*cv/(cv*theta*h2*(1-p11**(2*p11)*(1-p11)**(2*(1-p11)))+1)**2  !Taylor expansion
  !se=se*(h2*(1-p11**(2*p11)*(1-p11)**(2*(1-p11)))*cv/(1+h2*(1-p11**(2*p11)*(1-p11)**(2*(1-p11)))*theta*cv))/h2
  h2=h2*(1-p11**(2*p11)*(1-p11)**(2*(1-p11)))
  h2=h2*cv/(1+h2*theta*cv)
  print*,'Transformed R2 on the liabiltiy scale    : ',h2
  print*,'Transformed SE on the liability scale    : ',se
  print*,""

 end if


end subroutine


