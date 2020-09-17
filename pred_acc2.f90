subroutine pred_acc (h2,N,M,Me,k,p11,p12,ckv)

!**********************************************************************
! Expected prediction accuracy (also see ccgpa5.r in /R/Daetwyler/)
! S. Hong Lee  (2016)
!************************************************************************

implicit none

  real ( kind = 8 ) bound,mean,p,q,sd,x
  integer ( kind = 4 ) status,which

  integer::N,M
  real::h2,Me,k,p11,p12,ckv,iv,iv2,thd,zv,cv,theta,lamda,vgcc,Q_r2,r2
  real::cv2,theta2,h2o,h2o2,v102,v103,qv,cthd,czv,civ,p_case_s_int
  real::var_l_s,p_case_s,p_case_s11,var_l_s_,p_case_s_int_,p_case_s11_
  real::p_case_s12,p_case_s12_,estim_OR2,estim_OR3

   print*,"*************************************************************"
   print*,"h2  : Proportion of variance due to SNPs on the liability scale"
   print*,"N   : Sample size"
   print*,"M   : Number of SNPs"
   print*,"Me  : Effective number of chromosome segments"
   print*,"k   : Population prevalence"
   print*,"p   : Proportion of cases in training sample"
   print*,"p2  : Proportion of cases in validation sample"
   print*,"ckv : Top prop. of genetic profile scores in validation sample"
   print*,"*************************************************************"
   print*,""


   h2=h2*(real(M)/(real(Me)+real(M))) !adjusted heritability due to maker density
   !print*,h2

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
   lamda=real(N)/real(Me)
   !print*,lamda
   vgcc=h2*(1-h2*theta)     !#g variance on the liability in CC

   Q_r2=h2/(h2+1/lamda)     !R2 in quantitative traits

   r2=h2*zv**2/(h2*zv**2+(k*(1-k))**2/(lamda*p11*(1-p11))) !R2 in CC 

   cv2=(k*(1-k))**2/(zv**2*p12*(1-p12))    !#the spread sheet
   theta2=iv*((p12-k)/(1-k))*(iv*((p12-k)/(1-k))-thd)
   h2o=h2/(cv2-h2*theta2*cv2)
   h2o2=h2/cv2
   !print*,h2o2,h2,cv2

   !#auc from (r2*h2o2*cv2), i.e. r2 (0,1) on the liability scale
   v102=r2*h2o2*cv2
   qv=(iv-iv2)*v102/sqrt(v102*((1-v102*iv*(iv-thd))+(1-v102*iv2*(iv2-thd))))

p = huge ( p )
q = huge ( q )
x = qv
which=1  !Calculate P and Q from X, MEAN and SD;
call cdfnor ( which, p, q, x, mean, sd, status, bound )
   v103=p
   !print*,v102,qv,v103,r2,h2o2,cv2


p = 1D+00 - ckv 
q = 1.0D+00 - p
x = huge(x)
which=2  !Calculate X
call cdfnor ( which, p, q, x, mean, sd, status, bound )
  cthd=x
  czv=1/(sd*sqrt(2*3.141593)) * exp(-0.5*((x-mean)/sd)**2)
  civ=czv/ckv

  var_l_s=(1-civ*(civ-cthd))*h2*r2+(1-h2*r2)

  p_case_s_int=(thd-civ*(h2*r2)**.5)/var_l_s**.5

p = huge ( p )
q = huge ( q )
x = p_case_s_int
which=1  !Calculate P and Q from X, MEAN and SD;
call cdfnor ( which, p, q, x, mean, sd, status, bound )

  p_case_s11=(1-p)*ckv
  var_l_s_=(1+civ*(-civ+cthd))*(h2*r2)+(1-(h2*r2))
  !print*,var_l_s

  p_case_s_int_=(thd+civ*(h2*r2)**.5)/var_l_s_**.5

p = huge ( p )
q = huge ( q )
x = p_case_s_int_
which=1  !Calculate P and Q from X, MEAN and SD;
call cdfnor ( which, p, q, x, mean, sd, status, bound )

  p_case_s11_=(1-p)*ckv


p = 1D+00 - p_case_s11/ckv 
q = 1.0D+00 - p
x = huge(x)
which=2  !Calculate X
call cdfnor ( which, p, q, x, mean, sd, status, bound )

p = huge ( p )
q = huge ( q )
x = (thd-(x*var_l_s**.5-thd))/var_l_s**.5
which=1  !Calculate P and Q from X, MEAN and SD;
call cdfnor ( which, p, q, x, mean, sd, status, bound )

  p_case_s11_=(1-p)*ckv
  !print*,p_case_s11,p_case_s11_

  !#this is for top and bottom x % ************************************
  p_case_s12=ckv-p_case_s11
  p_case_s12_=ckv-p_case_s11_

  estim_OR2=(p_case_s11/p_case_s12)/(p_case_s11_/p_case_s12_)
  estim_OR3=(p_case_s11/p_case_s12)/(k/(1-k))

  print*,'Cor(g,g-hat) in quantitative trait       : ',Q_r2**.5
  print*,'Cor(y,g-hat) in quantitative trait       : ',(Q_r2*h2)**.5
  print*,'Cor(u,u-hat) in case-control data        : ',r2**.5
  print*,'Cor(y,u-hat) in case-control data        : ',(r2*h2o2)**.5
  print*,'Cor(y,u-hat) on the liability scale      : ',(r2*h2o2*cv2)**.5
  print*,'AUC in case-control data                 : ',v103
  print*,'OR contrasting top/bottom percentile     : ',estim_OR2
  print*,'   *** pA and pB for this OR             : ',p_case_s11/ckv,p_case_s11_/ckv
  print*,'OR contrasting top/general population    : ',estim_OR3
  print*,'   *** pA and pB for this OR             : ',p_case_s11/ckv,k

  print*,""

end subroutine


