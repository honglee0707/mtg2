subroutine pred_acc2 (fl_pred)
!subroutine pred_acc2 (h2,N,M,Me,k,p11,p12,ckv)

!**********************************************************************
! Expected prediction accuracy (also see ccgpa5.r in /R/Daetwyler/)
! S. Hong Lee  (2016)
!************************************************************************

implicit none

  real ( kind = 8 ) bound,mean,p,q,sd,x
  integer ( kind = 4 ) status,which

  integer::N,M,refn,zi,i,j
  real::h2,t_h2,t_k,Me,k,p11,p12,ckv,iv,iv2,thd,zv,cv,theta,lamda,vgcc,Q_r2,r2
  real::cv2,theta2,h2o,h2o2,v102,v103,qv,cthd,czv,civ,p_case_s_int
  real::var_l_s,p_case_s,p_case_s11,var_l_s_,p_case_s_int_,p_case_s11_
  real::p_case_s12,p_case_s12_,estim_OR2,estim_OR3
  character(len=64)::fl_pred

  integer,allocatable::r_N(:)
  real,allocatable::r_Me(:),r_h2(:),r_p11(:),r_Q_r2(:),r_r2(:),rG(:,:),r_k(:)

  real::v1,v2
  double precision,allocatable::pmat(:,:)
  double precision::dv

  open (unit=51,file=fl_pred,status='old')
  read(51,*)refn
  read(51,*)ckv
  read(51,*)M

  !print'(*(G0,2X))',' # ref sets  : ',refn
  print*,           'ckv         : ',ckv
  print'(*(G0,2X))',' M           : ',M
  print*,''

  read(51,*)t_h2,t_k,p12

  if (t_k==0 .and. p12==0) then
    print'(*(G0,2X))','Target set (QT)'
  elseif (t_k.ne.0 .and. p12.ne.0) then
    print'(*(G0,2X))','Target set (BT)'
  end if

  print*,           'h2          : ',t_h2
  print*,           'k           : ',t_k
  print*,           'P           : ',p12

  allocate(r_Me(refn),r_h2(refn),r_p11(refn),r_Q_r2(refn),r_r2(refn),r_N(refn))
  allocate(rG((1+refn),(1+refn)),pmat(refn,refn),r_k(refn))

  do i=1,refn
    read(51,*)r_h2(i),r_N(i),r_Me(i),r_k(i),r_p11(i)
    print*,''
    if (r_k(i)==0 .and. r_p11(i)==0) then
      print'(*(G0,2X))','Reference set',i,'(QT)'
    elseif (r_k(i).ne.0 .and. r_p11(i).ne.0) then
      print'(*(G0,2X))','Reference set',i,'(BT)'
    end if

    print*,           'h2          : ',r_h2(i)
    print'(*(G0,2X))',' N           : ',r_N(i)
    print'(*(G0,2X))',' Me          : ',r_Me(i)
    print*,           'k           : ',r_k(i)
    print*,           'P           : ',r_p11(i)
  end do

  print*,''
  print*,'lower triangle rG among target and reference sets'
  do i=1,refn
    read(51,*)(rG((1+i),j),j=1,i)
    print*,(rG((1+i),j),j=1,i)
  end do
 

  do zi=1,refn

    Me=r_Me(zi)
    N=r_N(zi)
    p11=r_p11(zi)
    h2=r_h2(zi)*(real(M)/(Me+real(M))) !adjusted h2 
    k=r_k(zi)
    lamda=real(N)/Me
    !print*,h2,p11

    if (k.ne.0 .and. p11.ne.0) then
      !CDF functions *****************************************************
      p = 1D+00 - k !0.011135D+00
      q = 1.0D+00 - p
      !p = huge ( p )
      !q = huge ( q )
      x = huge ( x )
      !x=-1.28155 !2 3.0D+00
      mean = 0.0D+00 ! 5.0D+00
      sd = 1.0D+00 !0.875D+00
      which=2  !Calculate X from P, Q, MEAN and SD;

      call cdfnor ( which, p, q, x, mean, sd, status, bound )
      thd=x
      zv=1/(sd*sqrt(2*3.141593)) * exp(-0.5*((x-mean)/sd)**2)
      !print*,thd,zv

      iv=zv/k          !#mean liability for cases
      iv2=-iv*k/(1-k)  !#mean liability for controls

      cv=(k*(1-k))**2/(zv**2*p11*(1-p11))    !#the spread sheet
      theta=iv*((p11-k)/(1-k))*(iv*((p11-k)/(1-k))-thd)
      vgcc=h2*(1-h2*theta)     !#g variance on the liability in CC

      r2=h2*zv**2/(h2*zv**2+(k*(1-k))**2/(lamda*p11*(1-p11))) !R2 in CC 
      r_r2(zi)=r2

    elseif (k==0 .and. p11==0) then
      Q_r2=h2/(h2+1/lamda)     !R2 in quantitative traits
      r_r2(zi)=Q_r2
    else
      print*,'for QT, k=P=0, and for BT, k>0 and P>0 >>> check'
    end if

    !print*,r_Q_r2(zi),r_r2(zi)
  end do !zi

  !P matrix for QT + BT 
  do i=1,refn
    pmat(i,i)=rG((1+i),1)*rG((1+i),1)*r_r2(i)
    do j=1,(i-1)
      pmat(i,j)=rG((1+i),(1+j))*rG((1+i),1)*rG((1+j),1)*r_r2(i)*r_r2(j)
      pmat(j,i)=pmat(i,j)
    end do
    !print*,pmat(i,:)
  end do
  call cholesky_inv (pmat,refn,dv) 

  v2=0
  do i=1,refn
    !print*,pmat(i,:)
    v1=0
    do j=1,refn
      v1=v1+rG((1+j),1)*rG((1+j),1)*r_r2(j)*pmat(j,i)
    end do
    v2=v2+v1*rG((1+i),1)*rG((1+i),1)*r_r2(i)
  end do
  !print*,v2
  r2=v2


  !In target ****************************************************
   h2=t_h2*(real(M)/(Me+real(M))) !adjusted h2 
   k=t_k

   if (k.ne.0 .and. p12.ne.0) then
      !CDF functions *****************************************************
      p = 1D+00 - k !0.011135D+00
      q = 1.0D+00 - p
      x = huge ( x )
      mean = 0.0D+00 ! 5.0D+00
      sd = 1.0D+00 !0.875D+00
      which=2  !Calculate X from P, Q, MEAN and SD;
      call cdfnor ( which, p, q, x, mean, sd, status, bound )
      !print*,p,q,x,mean,sd,status,bound
      thd=x
      zv=1/(sd*sqrt(2*3.141593)) * exp(-0.5*((x-mean)/sd)**2)

      iv=zv/k          !#mean liability for cases
      iv2=-iv*k/(1-k)  !#mean liability for controls

      !print*,thd
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

      !OR percentile
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


      print*,''
      print*,'Cor(u,u-hat) in the target data set      : ',r2**.5
      print*,'Cor(y,u-hat) in the target data set      : ',(r2*h2o2)**.5
      print*,'Cor(y,u-hat) on the liability scale      : ',(r2*h2o2*cv2)**.5
      print*,'AUC in case-control data                 : ',v103
      print*,'OR contrasting top/bottom percentile     : ',estim_OR2
      print*,'   *** pA and pB for this OR             : ',p_case_s11/ckv,p_case_s11_/ckv
      print*,'OR contrasting top/general population    : ',estim_OR3
      print*,'   *** pA and pB for this OR             : ',p_case_s11/ckv,k
      print*,""

    elseif (k==0 .and.p12==0) then

      print*,''
      print*,'Cor(g,g-hat) in the target data set      : ',r2**.5
      print*,'Cor(y,g-hat) in the target data set      : ',(r2*h2)**.5
      print*,""

    else
      print*,'for QT, k=P=0, and for BT, k>0 and P>0 >>> check'
    end if


end subroutine


