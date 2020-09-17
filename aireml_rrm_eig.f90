
subroutine aireml_rrm_eig (mbin,eval,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl12,xmm,nit,conv)        

!***********************************************************************
!aireml_rrm: random regression analysis  
!S. Hong Lee

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nb2,nb3,nbc,nbc2,xi,yi,yi2,ri,zi,vi,ui,vi2,ui2,wi,wj
!integer, parameter::tnk=4    !no. polynomial components
integer::tnk(mn)    !no. polynomial components
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,yidx1(rn),yidx2(rn),nit
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1),eval(pedn,mn)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),pig_vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
double precision::blup_ebv2(pedn,tn,mn) !,cor
integer::i,j,m,k,k2,l,io
double precision::sum_h2,sum_v(tn),flegendre
double precision::d_beta(tfn,1),d_res(trn,1),res(trn,1)

double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v3,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,1000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn),xvmy(trn,1)
double precision::py(trn,1),xb(trn,1)
!double precision::pig(trn,trn),pig_tmp(trn,trn)
double precision::pig(pedn,tn,tn),pig_tmp(pedn,tn,tn)

double precision::v_tmp(trn),v2_tmp(trn),tmp_pm(tn,tn)
double precision::dpy(trn,1),xvmdpy(tfn,1),pdpy(trn,1)

double precision::vmi(pedn,tn,tn),pm(pedn,tn,tn)

!for BLUP ebv
!double precision::vmat(mn,pedn,rn)

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

character(len=7)::cdum
character(len=64)::fl4,fl5,fl11,fl12
logical::file_exist


double precision::vcva2(tn,tn),am2(20,20),mm2(tn,20)
double precision,allocatable::am(:,:),mm(:,:),am_2(:,:),mm_2(:,:)
double precision::tkm(mn,20,20),tkm_rsv(mn,20,20)
double precision,allocatable::km(:,:),km2(:,:),phi(:,:),km2_2(:,:),phi_2(:,:)
double precision::wk(tn),xk(tn)
double precision,allocatable::km_phi(:,:)

double precision,allocatable::aim(:,:)
double precision,allocatable::sdmI(:,:)
double precision,allocatable::up(:,:),dldv(:,:)

integer::idum,idum2,ii,jj

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which

open (UNIT=46,FILE=fl12,STATUS='old')
read(46,*)tnk(:)
read(46,*)wk(:)
close(46)

do i=1,pedn
  yidx2(yidx(1,i))=i   !index for y to the fam order
end do

idum=tn
do ii=1,mn
  idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
end do

!print*,idum,tn+mn*(tnk(1)+(tnk(1)**2-tnk(1))/2)
allocate(aim(idum,idum),sdmI(idum,idum),up(idum,1),dldv(idum,1))


!xm=0;trn=0;tfn=0
!do xi=1,tn
!  tfn=tfn+fnm(xi)
!  do i=1,rn
!    if (yv(i,xi).ne.-99999) then
!      trn=trn+1
!      tyv(trn,1)=yv(i,xi)

!      do j=1,fnm(xi)
!        xm(trn,tfn-fnm(xi)+j)=xmm(xi,i,j)
!      end do
!    end if
!  end do
!end do

trn=0
do xi=1,tn
  do i=1,rnm(xi)
    trn=trn+1
    v1=0
    do j=1,rnm(xi)
      !v1=v1+mbin(1,yidx(xi,j),yidx(xi,i))*yv(yidx(xi,j),xi)
      v1=v1+mbin(1,yidx(xi,j),yidx(xi,i))*yv(j,xi)
    end do
    tyv(trn,1)=v1                !check
  end do
end do

xm=0;xi=1    !assuming all traits have the same n size, fn size
do j=1,fnm(xi)
  do i=1,rnm(xi)
    v1=0
    do k=1,rnm(xi)
      !v1=v1+mbin(1,yidx(xi,k),yidx(xi,i))*xmm(xi,yidx(xi,k),j)
      v1=v1+mbin(1,yidx(xi,k),yidx(xi,i))*xmm(xi,k,j)
    end do
    xm(i,j)=v1       !check
  end do
end do

tfn=fnm(1);trn=rnm(1)    !from second trait
do xi=2,tn
  tfn=tfn+fnm(xi)
  do i=1,rnm(xi)
    trn=trn+1
    xm(trn,(tfn-fnm(xi)+1):tfn)=xm(i,1:fnm(1))
    !print*,trn,tfn-fnm(xi)+1,tfn,fnm(1)
  end do
end do

print*,'*** number of records used ***'
do i=1,tn
  print*,'site',i,':',rnm(i)
end do
print*,''

  ! y = Xb + Zu + e ---- fitting effects without QTL

  !arbitary starting value (half of the total variance)
  vcve=0;vcva=0;tkm=0   !;km=0
  do xi=1,tn
    x1=0;x2=0
    do i=1,rn
      if (yv(i,xi).ne.-99999) then
        x1=x1+yv(i,xi)
        x2=x2+yv(i,xi)**2
      end if
    end do
    vcve(xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/2
    vcva(:,xi,xi)=((x2-(x1**2/rnm(xi)))/(rnm(xi)-1))/(2*mn)
  end do
  
  inquire(file=fl4,exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=fl4,STATUS='old')
    do i=1,tn
      read(44,*,iostat=io)cdum,vcve(i,i)
      if (io.ne.0) exit
    end do
    do wi=1,mn
      allocate(km(tnk(wi),tnk(wi)))

      do xi=1,tnk(wi)
        read(44,*,iostat=io)cdum,km(xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,tnk(wi)
        do yi=1,xi-1
          read(44,*,iostat=io)cdum,km(xi,yi)
          km(yi,xi)=km(xi,yi)
          if (io.ne.0) exit
        end do
      end do
      tkm(wi,1:tnk(wi),1:tnk(wi))=km

      deallocate(km)

    end do
    close(44)
  end if

  !Legendre polynomial function (k=4) ***********************


  am2=0

am2(1,1)=0.7071068

am2(2,2)=1.224745

am2(1,3)=-0.7905694
am2(3,3)=2.371708

am2(2,4)=-2.806243
am2(4,4)=4.677072

am2(1,5)=0.7954951
am2(3,5)=- 7.954951
am2(5,5)=9.280777

am2(2,6)=4.397265
am2(4,6)=- 20.52057
am2(6,6)=18.46851

am2(1,7)=-0.7967218
am2(3,7)=16.731162
am2(5,7)=- 50.19347
am2(7,7)= 36.80855

am2(2,8)=-5.990715
am2(4,8)= 53.91644
am2(6,8)=- 118.6162
am2(8,8)=+ 73.42906

am2(1,9)=0.7972005
am2(3,9)=- 28.69922
am2(5,9)=+ 157.8457
am2(7,9)=- 273.5992
am2(9,9)=+ 146.571

am2(2,10)=7.585119
am2(4,10)=- 111.2484
am2(6,10)=+ 433.8688
am2(8,10)=- 619.8126
am2(10,10)=+ 292.6893

am2(1,11)=-0.7974349
am2(3,11)=+ 43.85892
am2(5,11)=- 380.1106
am2(7,11)=+ 1140.332
am2(9,11)=- 1384.689
am2(11,11)=+584.6464


am2(2,12)=-9.17999
am2(4,12)=+ 198.8998
am2(6,12)=- 1193.399
am2(8,12)=+ 2898.254
am2(10,12)=- 3059.268
am2(12,12)=+ 1168.084

am2(1,13)=0.7975667
am2(3,13)=- 62.2102
am2(5,13)=+ 777.6276
am2(7,13)=- 3525.245
am2(9,13)=+ 7176.392
am2(11,13)=- 6697.965
am2(13,13)=+ 2334.139

am2(2,14)=10.77512
am2(4,14)=- 323.2537
am2(6,14)=+ 2747.657
am2(8,14)=- 9943.9
am2(10,14)=+ 17401.82
am2(12,14)=- 14554.25
am2(14,14)=+ 4664.825

am2(1,15)=-0.7976481
am2(3,15)=+ 83.75305
am2(5,15)=- 1423.802
am2(7,15)=+ 9017.412
am2(9,15)=- 27052.24
am2(11,15)=+ 41480.09
am2(13,15)=- 31424.31
am2(15,15)=+ 9323.698

am2(2,16)=-12.37042
am2(4,16)=+ 490.6933
am2(6,16)=- 5593.904
am2(8,16)=+ 27969.52
am2(10,16)=- 71477.66
am2(12,16)=+ 97469.54
am2(14,16)=- 67478.91
am2(16,16)=+ 18637.03

am2(1,17)=0.7977018
am2(3,17)=- 108.4874
am2(5,17)=+ 2404.805
am2(7,17)=- 20200.36
am2(9,17)=+ 82965.78
am2(11,17)=- 184368.4
am2(13,17)=+ 226270.3
am2(15,17)=- 144216.2
am2(17,17)=+ 37255.86


am2(2,18)=13.96582
am2(4,18)=- 707.6017
am2(6,18)=+ 10401.75
am2(8,18)=- 68354.33
am2(10,18)=+ 237341.4
am2(12,18)=- 466052.2
am2(14,18)=+ 519827.5
am2(16,18)=- 306945.8
am2(18,18)=+ 74479.49

am2(1,19)=-0.7977391
am2(3,19)=+ 136.4134
am2(5,19)=- 3819.575
am2(7,19)=+ 40996.77
am2(9,19)=- 219625.6
am2(11,19)=+ 658876.7
am2(13,19)=- 1158026
am2(15,19)=+ 1183477
am2(17,19)=- 650912.2
am2(19,19)=+ 148901.5

am2(2,20)=-15.5613
am2(4,20)=+ 980.362
am2(6,20)=- 18038.66
am2(8,20)=+ 150322.2
am2(10,20)=- 676449.8
am2(12,20)=+ 1783368
am2(14,20)=- 2835097
am2(16,20)=+ 2673092
am2(18,20)=- 1375856
am2(20,20)=+ 297699.8

  wk=wk-wk(1)
  xk=(wk-wk(tn)/2)/(wk(tn)/2)
 ! print*,xk
  !pause

  mm2=1
  do i=2,20
    mm2(:,i)=xk**(i-1)
  end do

  !polynomial components not same for all multiple random effects
  do wi=1,mn

    allocate(am(tnk(wi),tnk(wi)),mm(tn,tnk(wi)),phi(tn,tnk(wi)),km(tnk(wi),tnk(wi)))
    am=am2(1:tnk(wi),1:tnk(wi))
    mm=mm2(1:tn,1:tnk(wi))

    phi=matmul(mm,am)   !***************************************

    !do i=1,tn
    !  print*,mm(i,1:tnk(wi))
    !end do
    !print*,''
    !do i=1,tnk(wi)
    !  print*,am(i,1:tnk(wi))
    !end do
    !print*,''
    !do i=1,tn
    !  print*,phi(i,1:tnk(wi))
    !end do
    !pause


    km=tkm(wi,1:tnk(wi),1:tnk(wi))
    vcva2=matmul(matmul(phi,km),transpose(phi))
    vcva(wi,:,:)=vcva2
    deallocate(am,mm,phi,km)
  end do

  itit=1   !iteration in the iteration

  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000

  vmi=0;pm=0

  !print*,'iteration start'
  !do zi=1,200
  !do zi=1,nit
  do zi=1,10000000
    LKHP=LKH

  call cpu_time (t2_cpu)
  !polynomial components not same for all multiple random effects
  do wi=1,mn

    allocate(am(tnk(wi),tnk(wi)),mm(tn,tnk(wi)),phi(tn,tnk(wi)),km(tnk(wi),tnk(wi)))
    am=am2(1:tnk(wi),1:tnk(wi))
    mm=mm2(1:tn,1:tnk(wi))

    phi=matmul(mm,am)   !***************************************

    km=tkm(wi,1:tnk(wi),1:tnk(wi))
    vcva2=matmul(matmul(phi,km),transpose(phi))
    vcva(wi,:,:)=vcva2

    deallocate(am,mm,phi,km)
  end do

  LA=0
  do i=1,pedn     !assuming all ID have record
    do xi=1,tn
      do yi=1,xi
        v1=vcve(xi,yi)
        do wi=1,mn
          !v1=v1+vcva(wi,xi,yi)*eval(i,1)
          !v1=v1+vcva(wi,xi,yi)*eval(i,wi)
          v1=v1+vcva(wi,xi,yi)*eval(yidx(xi,i),wi)
        end do
        tmp_pm(xi,yi)=v1
        tmp_pm(yi,xi)=v1
      end do
    end do
    ! get inverse and determinant
    call cholesky_inv (tmp_pm,tn,x1)
    LA=LA+x1

    do xi=1,tn
      do yi=1,xi
        !vmi(pedn*xi-pedn+i,pedn*yi-pedn+i)=tmp_pm(xi,yi)
        !vmi(pedn*yi-pedn+i,pedn*xi-pedn+i)=tmp_pm(xi,yi)
        vmi(i,xi,yi)=tmp_pm(xi,yi)
        vmi(i,yi,xi)=tmp_pm(xi,yi)
      end do
    end do
  end do !i

  !call cpu_time (t2_cpu)
  !call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA

    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)
    do i=1,tfn
      do yi=1,tn
        do j=1,pedn
          v1=0
          do xi=1,tn
            !v1=v1+xm(pedn*xi-pedn+j,i)*vmi(pedn*xi-pedn+j,pedn*yi-pedn+j)
            v1=v1+xm(pedn*xi-pedn+j,i)*vmi(j,xi,yi)
          end do
          xvm(i,pedn*yi-pedn+j)=v1
        end do
      end do
    end do

    !xvmx=MATMUL(xvm,xm)
    do i=1,tfn
      do j=1,tfn
        v1=0
        do k=1,trn
          v1=v1+xvm(i,k)*xm(k,j)    !xm*pm
        end do
        xvmx(i,j)=v1
      end do
    end do

    !call cholesky_inv (xvmx,tfn,xvmxI,LG)
    call cholesky_inv (xvmx,tfn,LG)
    !BLUE = (X'VI X)I XVI y **********************
    !beta=matmul(matmul(matmul(xvmx,transpose(xm)),pm),tyv)
    do i=1,tfn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvmx(i,k)*xm(j,k)    !xvmx*xm
        end do
        xvm2(i,j)=v1
      end do
    end do

    xvmy=0
    do i=1,trn
      v1=0
      do j=1,tfn
        xvmy(j,1)=xvmy(j,1)+xvm(j,i)*tyv(i,1)
      end do
    end do
   beta=matmul(xvmx,xvmy)
   !print*,beta

    !residual
    res=tyv-matmul(xm,beta)

    !py
    do yi=1,tn
      do i=1,pedn
        v1=0
        do xi=1,tn
          !v1=v1+vmi(pedn*yi-pedn+i,pedn*xi-pedn+i)*res(pedn*xi-pedn+i,1)
          v1=v1+vmi(i,yi,xi)*res(pedn*xi-pedn+i,1)
        end do
        py(pedn*yi-pedn+i,1)=v1
      end do
    end do

    !(XV-1X)-1XV-1
    do i=1,tfn
      do j=1,pedn
        do xi=1,tn
          v1=0
          do k=1,tn
            !v1=v1+xvm2(i,k*pedn-pedn+j)*vmi2(k*pedn-pedn+j,xi*pedn-pedn+j)
            v1=v1+xvm2(i,k*pedn-pedn+j)*vmi(j,k,xi)
          end do
          xvm3(xi*pedn-pedn+j,i)=v1
        end do
      end do
    end do

    do yi=1,tn
      do i=1,pedn
        do xi=1,tn
          v1=0
          do k=1,tfn
            v1=v1+xvm3(yi*pedn-pedn+i,k)*xvm(k,xi*pedn-pedn+i)
          end do
          !pm2(yi*pedn-pedn+i,xi*pedn-pedn+i)=vmi2(yi*pedn-pedn+i,xi*pedn-pedn+i)-v1
          pm(i,yi,xi)=vmi(i,yi,xi)-v1
        end do
      end do
    end do

    !ypy estimation
    !ypyf=MATMUL(MATMUL(transpose(tyv),pm),tyv)
    v2=0
    do i=1,trn
      v2=v2+py(i,1)*tyv(i,1)
    end do

    LKH=-0.5*(LA+LG+v2)
    !print*,LKH,LA,LG,v2  !ypyf
    !pause
    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    if ((LKH>-99999999 .and. LKH<99999999) .and. (LKH>=LKHP)) then
      itit=1
      !print*,'ok',itit,LKHP
    else
      itit=itit+1
      LKH=LKHP
      !print*,'likelihood NaN >> reduce by 0.7^',itit-1
      print*,itit-1,'likelihood nan >> update reduced by the factor'
      !print*,'not ok',itit,LKHP
      goto 111
    end if

    PRINT '(a7,100f14.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f14.4)','Ve',vcve(xi,xi)
    end do
    print*,''
    do wi=1,mn
      do xi=1,tnk(wi)
        PRINT '(a7,100f14.4)','Vk',tkm(wi,xi,xi)
      end do
      do xi=2,tnk(wi)
        do yi=1,xi-1
          PRINT '(a7,100f14.4)','cov',tkm(wi,xi,yi)
        end do
      end do
      print*,''
    end do
    !pause

  !print*,'trasnformed ***'
  !do wi=1,mn
  !  do xi=1,tn
  !    print '(a7,100f14.4)','Va',vcva(wi,xi,xi)!,sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
  !    sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
  !  end do
  !  nbc=0
  !  do xi=2,tn
  !    do yi=1,xi-1
  !      nbc=nbc+1
  !      print '(a7,100f14.4)','cov',vcva(wi,xi,yi)!,sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
  !    end do
  !  end do
  !  write(41,*)''
  !end do



    ! AI matrix
    !py=MATMUL(pm,tyv)       !for py

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    nb=tnk(1)+(tnk(1)**2-tnk(1))/2  !no block for random effects (Vg, cov) check
    !Ve (diagonal)  **********************************
    !eig simplificaiton
    !making dPy (see thompson3.xlsx) 
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)
      dpy=0
      xvmdpy=0  !xvm %*% dpy
      do i=1,rnm(xi)
        dpy(trnx-rnm(xi)+i,1)=py(trnx-rnm(xi)+i,1)
        v1=0
        do j=1,tfn
          xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,trnx-rnm(xi)+i)*dpy(trnx-rnm(xi)+i,1)
        end do
      end do
      d_beta=matmul(xvmx,xvmdpy)   !beta with dpy as response variable
      d_res=dpy-matmul(xm,d_beta)  !res with dpy as repsonse variable

      !making PdPy
      v3=0
      do ui=1,tn
        do i=1,pedn
          v2=0
          do vi=1,tn
            !v2=v2+vmi(ui*pedn-pedn+i,vi*pedn-pedn+i)*d_res(vi*pedn-pedn+i,1)
            v2=v2+vmi(i,ui,vi)*d_res(vi*pedn-pedn+i,1)
          end do
          pdpy(ui*pedn-pedn+i,1)=v2
          v3=v3+dpy(ui*pedn-pedn+i,1)*v2
        end do
      end do
      aim(xi,xi)=v3
      !print*,v3

      !Ve (off diagonal) aim(2,1) *************************************
      trny=0
      do yi=1,(xi-1)
        trny=trny+rnm(yi)  !total rn upto yi
        v3=0
        do i=1,rnm(yi)
          !dpy(trnx-rnm(xi)+i,1)=py(trnx-rnm(xi)+i,1)
          v3=v3+py(trny-rnm(yi)+i,1)*pdpy(trny-rnm(yi)+i,1)
        end do
        aim(xi,yi)=v3
        !print*,v3
      end do !yi
    end do !xi

    !dldv Ve ***********************************************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      tr1=0
      do i=1,rnm(xi)
        !tr1=tr1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
        tr1=tr1+pm(i,xi,xi)
      end do
      tr1=-0.5*tr1
      v2=0
      do i=1,rnm(xi)
        v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
      end do
      dldv(xi,1)=v2*(0.5)+tr1
      !print*,v2,tr1
    end do

    !print*,'Vg'
    !for Vg
    nb2=tn
    do k=1,mn !no. random effects
      !for Vg (doagonal) ******************************************************* 
      allocate(am(tnk(k),tnk(k)),mm(tn,tnk(k)),km2(tnk(k),tnk(k)),phi(tn,tnk(k)))
      am=am2(1:tnk(k),1:tnk(k))
      mm=mm2(1:tn,1:tnk(k))

      phi=matmul(mm,am)  

      do xi=1,tnk(k)   !check
        nb2=nb2+1

        km2=0;km2(xi,xi)=1;pig_vcva=0
        vcva2=matmul(matmul(phi,km2),transpose(phi))
        pig_vcva(k,:,:)=vcva2

        !do i=1,tnk(k)
        !  print*,km2(i,1:tnk(k))
        !end do
        !print*,''
        !do i=1,tn
        !  print*,phi(i,1:tnk(k))
        !end do
        !print*,''
        !do i=1,tn
        !  print*,vcva2(i,:)
        !end do
        !pause

        !for d (i.e. dV/di) *py
        !pig=0
        !do i=1,pedn     !assuming all ID have record
        !  do ui=1,tn
        !    v2=0
        !    do yi=1,tn
        !      v1=pig_vcva(k,ui,yi)*eval(i,k)
        !      pig(pedn*ui-pedn+i,pedn*yi-pedn+i)=v1
        !      !pig(i,ui,yi)=v1
        !      v2=v2*py(yi*pedn-pedn+i,1)
        !    end do
        !    dpy(pedn*ui-pedn+i,1)=v1
        !  end do
        !end do

        !dpy=matmul(pig,py)
        !print*,dpy(1:10,1)

        pig=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do yi=1,tn
              !v2=pig_vcva(k,ui,yi)*eval(i,k)
              v2=pig_vcva(k,ui,yi)*eval(yidx(ui,i),k)
              !pig(pedn*ui-pedn+i,pedn*yi-pedn+i)=v2
              pig(i,ui,yi)=v2
              v1=v1+v2*py(yi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do
        !print*,dpy(1:10,1)
        !pause

        xvmdpy=0
        do i=1,trn 
          do j=1,tfn
            xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,i)*dpy(i,1)
          end do
        end do
        d_beta=matmul(xvmx,xvmdpy)   !beta with dpy as response variable
        d_res=dpy-matmul(xm,d_beta)  !res with dpy as repsonse variable

        !making PdPy
        v3=0
        do ui=1,tn
          do i=1,pedn
            v2=0
            do vi=1,tn
              v2=v2+vmi(i,ui,vi)*d_res(vi*pedn-pedn+i,1)
            end do
            pdpy(ui*pedn-pedn+i,1)=v2
            v3=v3+dpy(ui*pedn-pedn+i,1)*v2
          end do
        end do
        !print*,v3

        idum=tn
        do ii=1,k-1
          idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
        end do
        aim(idum+xi,idum+xi)=v3
        !print*,tmpm(1,1) !idum+xi,tn+nb*k-nb+xi  !'tmpm:',tmpm(1,1)/2
        !pause

        !for Vg x Ve (off diagonal) !******************************************
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)     !*dpy
          end do
          aim(idum+xi,vi)=v2
        end do

        !for Vg(1~[i-1]) x Vg(i)
        do vi=1,xi-1
          km2=0;km2(vi,vi)=1;pig_vcva=0
          vcva2=matmul(matmul(phi,km2),transpose(phi))
          pig_vcva(k,:,:)=vcva2

          !for d (i.e. dV/di) for Vg(1~[i-1])
          !pig_tmp=0
          !do i=1,pedn     !assuming all ID have record
          !  do ui=1,tn
          !    do yi=1,ui
          !      v1=0 !vcve(xi,yi)
          !      do k2=1,mn   !check
          !        v1=v1+pig_vcva(k2,ui,yi)*eval(i,1)
          !      end do
          !      pig_tmp(pedn*ui-pedn+i,pedn*yi-pedn+i)=v1
          !      pig_tmp(pedn*yi-pedn+i,pedn*ui-pedn+i)=v1
          !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
          !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
          !    end do
          !  end do
          !end do
          !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do yi=1,tn
              !v2=pig_vcva(k,ui,yi)*eval(i,k)
              v2=pig_vcva(k,ui,yi)*eval(yidx(ui,i),k)
              !pig_tmp(pedn*ui-pedn+i,pedn*yi-pedn+i)=v2
              pig_tmp(i,ui,yi)=v2
              v1=v1+v2*py(yi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do


          v2=0
          do i=1,trn
            v2=v2+pdpy(i,1)*dpy(i,1)
          end do
          aim(idum+xi,idum+vi)=v2
        end do


        !off diagonal *********************************************************
        nb3=tn
        do wi=1,k-1  !check in case multiple random effects
          allocate(am_2(tnk(wi),tnk(wi)),mm_2(tn,tnk(wi)))
          allocate(km2_2(tnk(wi),tnk(wi)),phi_2(tn,tnk(wi)))

          am_2=am2(1:tnk(wi),1:tnk(wi))
          mm_2=mm2(1:tn,1:tnk(wi))

          phi_2=matmul(mm_2,am_2)   
          do vi=1,tnk(wi)
            nb3=nb3+1

            km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
            vcva2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
            pig_vcva(wi,:,:)=vcva2

            !for d (i.e. dV/di)
            !pig_tmp=0
            !do i=1,pedn     !assuming all ID have record
            !  do ui=1,tn
            !    do yi=1,ui
            !      v1=0 !vcve(xi,yi)
            !      do k2=1,mn   !check
            !        v1=v1+pig_vcva(k2,ui,yi)*eval(i,1)
            !      end do
            !      pig_tmp(pedn*ui-pedn+i,pedn*yi-pedn+i)=v1
            !      pig_tmp(pedn*yi-pedn+i,pedn*ui-pedn+i)=v1
            !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
            !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
            !    end do
            !  end do
            !end do

            !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do yi=1,tn
              !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
              v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
              !pig_tmp(pedn*ui-pedn+i,pedn*yi-pedn+i)=v2
              pig_tmp(i,ui,yi)=v2
              v1=v1+v2*py(yi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do





            v2=0
            do i=1,trn
              v2=v2+pdpy(i,1)*dpy(i,1)
            end do

            idum2=tn
            do ii=1,wi-1
              idum2=idum2+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
            end do
            aim(idum+xi,idum2+vi)=v2   !with the xi th Vg for first trait
          end do

          nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
          do vi2=2,tnk(wi)
            do ui2=1,(vi2-1)
              nbc=nbc+1

              km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
              vcva2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
              pig_vcva(wi,:,:)=vcva2

              !for d (i.e. dV/di)
              !pig_tmp=0
              !do i=1,pedn     !assuming all ID have record
              !  do ui=1,tn
              !    do yi=1,ui
              !      v1=0 !vcve(xi,yi)
              !      do k2=1,mn   !check
              !        v1=v1+pig_vcva(k2,ui,yi)*eval(i,1)
              !      end do
              !      pig_tmp(pedn*ui-pedn+i,pedn*yi-pedn+i)=v1
              !      pig_tmp(pedn*yi-pedn+i,pedn*ui-pedn+i)=v1
              !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
              !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
              !    end do
              !  end do
              !end do

              !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do yi=1,tn
              !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
              v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
              !pig_tmp(pedn*ui-pedn+i,pedn*yi-pedn+i)=v2
              pig_tmp(i,ui,yi)=v2
              v1=v1+v2*py(yi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do




              v2=0
              do i=1,trn
                v2=v2+pdpy(i,1)*dpy(i,1)
              end do

              aim(idum+xi,idum2+tnk(wi)+nbc)=v2  !check
            end do
          end do
          deallocate(km2_2,phi_2,am_2,mm_2)
        end do !wi

        !dldv Vg *********************************************************
        tr1=0
        do i=1,pedn
          do ui=1,tn
            do yi=1,tn
              !tr1=tr1+pm(pedn*ui-pedn+i,pedn*yi-pedn+i)*pig(pedn*ui-pedn+i,pedn*yi-pedn+i)
              !tr1=tr1+pm(i,ui,yi)*pig(pedn*ui-pedn+i,pedn*yi-pedn+i)
              !tr1=tr1+pm(i,ui,yi)*pig(i,ui,yi)
              tr1=tr1+pm(i,ui,yi)*pig(i,yi,ui)
            end do
          end do
        end do
        tr1=-0.5*tr1

        v2=0
        do i=1,trn
          v2=v2+pdpy(i,1)*tyv(i,1)
        end do

        dldv(idum+xi,1)=v2*0.5+tr1 !************************************
      end do !xi


      !for cov (diagonal) *************************************************
      nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
      do xi=2,tnk(k)
        do yi=1,xi-1
          nbc=nbc+1

          km2=0;km2(xi,yi)=1;km2(yi,xi)=1;pig_vcva=0
          vcva2=matmul(matmul(phi,km2),transpose(phi))
          pig_vcva(k,:,:)=vcva2

          !for d (i.e. dV/di)
          !pig=0
          !do i=1,pedn     !assuming all ID have record
          !  do ui=1,tn
          !    do vi=1,ui
          !      v1=0 !vcve(xi,yi)
          !      do k2=1,mn   !check
          !        v1=v1+pig_vcva(k2,ui,vi)*eval(i,1)
          !      end do
          !      pig(pedn*ui-pedn+i,pedn*vi-pedn+i)=v1
          !      pig(pedn*vi-pedn+i,pedn*ui-pedn+i)=v1
          !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
          !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
          !    end do
          !  end do
          !end do

          !dpy=matmul(pig,py)

        pig=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do vi=1,tn
              !v2=pig_vcva(k,ui,vi)*eval(i,k)
              v2=pig_vcva(k,ui,vi)*eval(yidx(ui,i),k)
              !pig(pedn*ui-pedn+i,pedn*vi-pedn+i)=v2
              pig(i,ui,vi)=v2
              v1=v1+v2*py(vi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do





          xvmdpy=0
          do i=1,trn
            do j=1,tfn
              xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,i)*dpy(i,1)
            end do
          end do
          d_beta=matmul(xvmx,xvmdpy)   !beta with dpy as response variable
          d_res=dpy-matmul(xm,d_beta)  !res with dpy as repsonse variable

          !making PdPy
          v3=0
          do ui=1,tn
            do i=1,pedn
              v2=0
              do vi=1,tn
                v2=v2+vmi(i,ui,vi)*d_res(vi*pedn-pedn+i,1)
              end do
              pdpy(ui*pedn-pedn+i,1)=v2
              v3=v3+dpy(ui*pedn-pedn+i,1)*v2
            end do
          end do
          !print*,v3

          idum=tn
          do ii=1,k-1
            idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
          end do
          aim(idum+tnk(k)+nbc,idum+tnk(k)+nbc)=v3     !check

          !for cov x Ve (off diagonal) ************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)     !*dpy
            end do
            aim(idum+tnk(k)+nbc,vi)=v2          !check with Ve1
          end do

          !for cov x Vg **************************************************
          do vi=1,tnk(k)

            km2=0;km2(vi,vi)=1;pig_vcva=0
            vcva2=matmul(matmul(phi,km2),transpose(phi))
            pig_vcva(k,:,:)=vcva2

            !for d (i.e. dV/di)
            !pig_tmp=0
            !do i=1,pedn     !assuming all ID have record
            !  do ui=1,tn
            !    do wi=1,ui
            !      v1=0 !vcve(xi,yi)
            !      do k2=1,mn   !check
            !        v1=v1+pig_vcva(k2,ui,wi)*eval(i,1)
            !      end do
            !      pig_tmp(pedn*ui-pedn+i,pedn*wi-pedn+i)=v1
            !      pig_tmp(pedn*wi-pedn+i,pedn*ui-pedn+i)=v1
            !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
            !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
            !    end do
            !  end do
            !end do

            !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do wi=1,tn
              !v2=pig_vcva(k,ui,wi)*eval(i,k)
              v2=pig_vcva(k,ui,wi)*eval(yidx(ui,i),k)
              !pig_tmp(pedn*ui-pedn+i,pedn*wi-pedn+i)=v2
              pig_tmp(i,ui,wi)=v2
              v1=v1+v2*py(wi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do




            v2=0
            do i=1,trn
              v2=v2+pdpy(i,1)*dpy(i,1)
            end do

            idum=tn
            do ii=1,k-1
              idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
            end do              
            aim(idum+tnk(k)+nbc,idum+vi)=v2  !with the same Vg1 for first trait
          end do

          !for cov x cov *************************************************
          trnv=rnm(1)
          nbc2=0    !no block for cov (vi,ui)
          do vi=2,tnk(k)
            trnv=trnv+rnm(vi)
            trnu=0
            do ui=1,(vi-1)
              nbc2=nbc2+1

              if (nbc2<nbc) then  !off diagonal within cova ************************
                trnu=trnu+rnm(ui)

                km2=0;km2(vi,ui)=1;km2(ui,vi)=1;pig_vcva=0
                vcva2=matmul(matmul(phi,km2),transpose(phi))
                pig_vcva(k,:,:)=vcva2

                !for d (i.e. dV/di)
                !pig_tmp=0
                !do i=1,pedn     !assuming all ID have record
                !  do wj=1,tn
                !    do wi=1,wj
                !      v1=0 !vcve(xi,yi)
                !      do k2=1,mn   !check
                !        v1=v1+pig_vcva(k2,wj,wi)*eval(i,1)
                !      end do
                !      pig_tmp(pedn*wj-pedn+i,pedn*wi-pedn+i)=v1
                !      pig_tmp(pedn*wi-pedn+i,pedn*wj-pedn+i)=v1
                !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
                !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
                !    end do
                !  end do
                !end do

                !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do wj=1,tn
          do i=1,pedn
            v1=0
            do wi=1,tn
              !v2=pig_vcva(k,wj,wi)*eval(i,k)
              v2=pig_vcva(k,wj,wi)*eval(yidx(wj,i),k)
              !pig_tmp(pedn*wj-pedn+i,pedn*wi-pedn+i)=v2
              pig_tmp(i,wj,wi)=v2
              v1=v1+v2*py(wi*pedn-pedn+i,1)
            end do
            dpy(pedn*wj-pedn+i,1)=v1
          end do
        end do



                v2=0
                do i=1,trn
                  v2=v2+pdpy(i,1)*dpy(i,1)
                end do

                idum=tn
                do ii=1,k-1
                  idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
                end do
                aim(idum+tnk(k)+nbc,idum+tnk(k)+nbc2)=v2  !with the wi th cov for first trait
              end if ! *************************************************************
            end do !ui
          end do !vi


          !off diagonal *********************************************************
          do wi=1,k-1   !check in case of multiple random effects
            allocate(am_2(tnk(wi),tnk(wi)),mm_2(tn,tnk(wi)),km2_2(tnk(wi),tnk(wi)),phi_2(tn,tnk(wi)))

            am_2=am2(1:tnk(wi),1:tnk(wi))
            mm_2=mm2(1:tn,1:tnk(wi))

            phi_2=matmul(mm_2,am_2)   !***************************************
            trnv=0
            do vi=1,tnk(wi)
              trnv=trnv+rnm(vi)

              km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
              vcva2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
              pig_vcva(wi,:,:)=vcva2

              !for d (i.e. dV/di)
              !pig_tmp=0
              !do i=1,pedn     !assuming all ID have record
              !  do ui=1,tn
              !    do yi2=1,ui
              !      v1=0 !vcve(xi,yi)
              !      do k2=1,mn   !check
              !        v1=v1+pig_vcva(k2,ui,yi2)*eval(i,1)
              !      end do
              !      pig_tmp(pedn*ui-pedn+i,pedn*yi2-pedn+i)=v1
              !      pig_tmp(pedn*yi2-pedn+i,pedn*ui-pedn+i)=v1
              !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
              !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
              !    end do
              !  end do
              !end do

              !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do yi2=1,tn
              !v2=pig_vcva(wi,ui,yi2)*eval(i,wi)
              v2=pig_vcva(wi,ui,yi2)*eval(yidx(ui,i),wi)
              !pig_tmp(pedn*ui-pedn+i,pedn*yi2-pedn+i)=v2
              pig_tmp(i,ui,yi2)=v2
              v1=v1+v2*py(yi2*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do





              v2=0
              do i=1,trn
                v2=v2+pdpy(i,1)*dpy(i,1)
              end do

              idum=tn
              do ii=1,k-1
                idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
              end do
              idum2=tn
              do ii=1,wi-1
                idum2=idum2+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
              end do
              aim(idum+tnk(k)+nbc,idum2+vi)=v2  !with the wi th Vg for first trait
            end do

            nbc2=0    !no block for cov (vi,ui)
            do vi2=2,tnk(wi)
              do ui2=1,(vi2-1)
                nbc2=nbc2+1

                km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                vcva2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                pig_vcva(wi,:,:)=vcva2

                !for d (i.e. dV/di)
                !pig_tmp=0
                !do i=1,pedn     !assuming all ID have record
                !  do ui=1,tn
                !    do vi=1,ui
                !      v1=0 !vcve(xi,yi)
                !      do k2=1,mn   !check
                !        v1=v1+pig_vcva(k2,ui,vi)*eval(i,1)
                !      end do
                !      pig_tmp(pedn*ui-pedn+i,pedn*vi-pedn+i)=v1
                !      pig_tmp(pedn*vi-pedn+i,pedn*ui-pedn+i)=v1
                !      !vmi(i,xi,yi)=tmp_pm(xi,yi)
                !      !vmi(i,yi,xi)=tmp_pm(xi,yi)
                !    end do
                !  end do
                !end do

                !dpy=matmul(pig_tmp,py)

        pig_tmp=0
        do ui=1,tn
          do i=1,pedn
            v1=0
            do vi=1,tn
              !v2=pig_vcva(wi,ui,vi)*eval(i,wi)
              v2=pig_vcva(wi,ui,vi)*eval(yidx(ui,i),wi)
              !pig_tmp(pedn*ui-pedn+i,pedn*vi-pedn+i)=v2
              pig_tmp(i,ui,vi)=v2
              v1=v1+v2*py(vi*pedn-pedn+i,1)
            end do
            dpy(pedn*ui-pedn+i,1)=v1
          end do
        end do





                v2=0
                do i=1,trn
                  v2=v2+pdpy(i,1)*dpy(i,1)
                end do

                idum=tn
                do ii=1,k-1
                  idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
                end do
                idum2=tn
                do ii=1,wi-1
                  idum2=idum2+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
                end do

                aim(idum+tnk(k)+nbc,idum2+tnk(wi)+nbc2)=v2  !with the wi th cov for first trait
              end do !ui
            end do !vi

            deallocate(am_2,mm_2,phi_2,km2_2)

          end do !wi

          !dldv for cov *****************************
          tr1=0
          do i=1,pedn
            do ui=1,tn
              do vi=1,tn
                !tr1=tr1+pm(pedn*ui-pedn+i,pedn*yi-pedn+i)*pig(pedn*ui-pedn+i,pedn*yi-pedn+i)
                !tr1=tr1+pm(i,ui,vi)*pig(pedn*ui-pedn+i,pedn*vi-pedn+i)
                !tr1=tr1+pm(i,ui,vi)*pig(i,ui,vi)
                tr1=tr1+pm(i,ui,vi)*pig(i,vi,ui)
              end do
            end do
          end do
          tr1=-0.5*tr1

          v2=0
          do i=1,trn
            v2=v2+pdpy(i,1)*tyv(i,1)
          end do

          idum=tn
          do ii=1,k-1
            idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
          end do
          dldv(idum+tnk(k)+nbc,1)=v2*0.5+tr1
        end do !yi
      end do !xi

      deallocate(am,mm,phi,km2)
    end do !k

    call cpu_time (t2_cpu)
    print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    print*,''
 
    idum=tn
    do ii=1,mn
      idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
    end do

    do i=1,idum
      aim(i,i)=aim(i,i)/2
      do j=1,i-1
        aim(i,j)=aim(i,j)/2
        aim(j,i)=aim(i,j)
      end do
    end do

    !do i=1,idum
    !  print*,i,real(aim(i,1:i))
    !end do
    !print*,''

    !call cholesky_inv (aim,(2+mn*3),x1)
    !call cholesky_inv (aim,(tn+nb*mn),x1)
    call cholesky_inv (aim,idum,x1)
    print *,''

    if (LKH-MLKH.gt.0) then
      MLKH=LKH
      mva=va
      mva2=va2
      mve=ve
      mve2=ve2
      mcov=cov
      sdmI=aim
    end if

    !if (ABS(LKH-LKHP).lt.0.001) goto 1000
    if (ABS(LKH-LKHP).lt.conv) goto 1000


    up=MATMUL(aim,dldv)
    !do i=1,idum
    !  print*,real(up(i,1)),real(dldv(i,1))
    !end do
    !pause

111 continue
    if (itit==1) then
      vcve_rsv=vcve
      tkm_rsv=tkm
      if (zi.ge.nit) then
        print*,'Likelihood not converged >>> may need a longer iterations'
        print*,''
        goto 1000
      end if

    else
      !up=up*0.1**(itit-1)
      up=up*0.7**(itit-1)
      vcve=vcve_rsv
      tkm=tkm_rsv
    end if

    do xi=1,tn
      vcve(xi,xi)=vcve(xi,xi)+up(xi,1)
      !if (vcve(xi,xi)<0) vcve(xi,xi)=0.000001
    end do
    do wi=1,mn

      idum=tn
      do ii=1,wi-1
        idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
      end do
      do xi=1,tnk(wi)
        !tkm(wi,xi,xi)=tkm(wi,xi,xi)+up(tn+nb*wi-nb+xi,1)
        tkm(wi,xi,xi)=tkm(wi,xi,xi)+up(idum+xi,1)
        !if (tkm(wi,xi,xi)<0) tkm(wi,xi,xi)=0.000001
      end do
      nbc=0
      do xi=2,tnk(wi)
        do yi=1,xi-1
          nbc=nbc+1
          !tkm(wi,xi,yi)=tkm(wi,xi,yi)+up(tn+nb*wi-nb+tnk(wi)+nbc,1)  !check
          tkm(wi,xi,yi)=tkm(wi,xi,yi)+up(idum+tnk(wi)+nbc,1)  !check
          tkm(wi,yi,xi)=tkm(wi,xi,yi)

        end do
      end do
    end do

    !print*,ve,ve2,va,va2,cov
    !print*,vcve,vcva

    !do i=1,tnk
    !  print*,km(i,:)
    !end do
    !pause

  END do !zi



  !PRINT*,'LKH not converged na na na'
  !WRITE(41,*)'LKH not converged na na na'
  va=mva;va2=mva2;ve=mve;ve2=mve2;LKH=MLKH;cov=mcov


1000 continue 


  open (UNIT=41,FILE=fl5,STATUS='unknown')

  do xi=1,tn
    write(41,'(a7,100f14.4)')'Ve',vcve(xi,xi),sqrt(sdmI(xi,xi))
  end do
  write(41,*)''

  do xi=1,tn
    sum_v(xi)=vcve(xi,xi)
  end do
  do wi=1,mn
    idum=tn
    do ii=1,wi-1
      idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
    end do

    do xi=1,tnk(wi)
      !write(41,'(a7,100f14.4)')'Va',vcva(wi,xi,xi),sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      !write(41,'(a7,100f14.4)')'Vk',km(xi,xi),sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      write(41,'(a7,100f14.4)')'Vk',tkm(wi,xi,xi),sqrt(sdmI(idum+xi,idum+xi))
      !sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tnk(wi)
      do yi=1,xi-1
        nbc=nbc+1
        !write(41,'(a7,100f14.4)')'cov',vcva(wi,xi,yi),sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
        !write(41,'(a7,100f14.4)')'cov',km(xi,yi),sqrt(sdmI(tn+wi*nb-nb+tnk(wi)+nbc,tn+wi*nb-nb+tnk(wi)+nbc))
        write(41,'(a7,100f14.4)')'cov',tkm(wi,xi,yi),sqrt(sdmI(idum+tnk(wi)+nbc,idum+tnk(wi)+nbc))
      end do
    end do
    write(41,*)''
  end do

  !write(41,'(a7,100f14.4)')'LKH',LKH
  write(41,'(a7,100f14.4)')'LKH',LKH !,real(zi) !indicate convergence
  write(41,*)''

  write(41,*)'transformed ***'
  do wi=1,mn
    do xi=1,tn
      write(41,'(a7,100f14.4)')'Va',vcva(wi,xi,xi)!,sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f14.4)')'cov',vcva(wi,xi,yi)!,sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
      end do
    end do
    write(41,*)''
  end do



  !!total Var(V) *********************
  !do xi=1,tn
  !  vvc(xi)=sdmI(xi,xi)
  !  do wi=1,mn
  !    vvc(xi)=vvc(xi)+sdmI(tn+nb*wi-nb+xi,xi)*2
  !    do wj=1,mn
  !      vvc(xi)=vvc(xi)+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
  !    end do
  !  end do
  !end do
  !!print*,vvc

  do wi=1,mn
    do xi=1,tn
      x1=vcva(wi,xi,xi)/sum_v(xi)
   !   x11=sdmI(tn+nb*wi-nb+xi,xi) !for the row ************

  !    do wj=1,mn
  !      x11=x11+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
  !    end do
      !print*,x11

  !    !SE **********************************
  !    x2=x1**2*(sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      write(41,'(a7,100f14.4)')'h2',x1!,sqrt(x2)

    end do
    write(41,*)''

    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
  !      !var(cor) *************************************
  !      z2=sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/(4*vcva(wi,xi,xi)**2)
  !      z2=z2+sdmI(tn+nb*wi-nb+yi,tn+nb*wi-nb+yi)/(4*vcva(wi,yi,yi)**2)
  !      z2=z2+sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+tn+nbc)/(vcva(wi,xi,yi)**2)
  !      z2=z2+2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+yi)/(4*vcva(wi,xi,xi)*vcva(wi,yi,yi))
  !      z2=z2-2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+tn+nbc)/(2*vcva(wi,xi,xi)*vcva(wi,xi,yi))
  !      z2=z2-2*sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+yi)/(2*vcva(wi,xi,yi)*vcva(wi,yi,yi))
  !      z2=z2*(z1**2)
        write(41,'(a7,100f14.4)')'cor',z1!,sqrt(z2)
      end do
    end do
    write(41,*)''
  end do


  idum=tn
  do ii=1,mn
    idum=idum+tnk(ii)+(tnk(ii)**2-tnk(ii))/2
  end do
  write(41,*)''
  do i=1,idum
    write(41,'(100g18.8)')(sdmI(i,1:i))
  end do

  close(41)

if (fl11.ne.'null') then

  print*,'for BLUP solution after convergence *******************'

  !make it triangle -> full ****************
  do xi=1,tn
    do yi=1,xi
      do k=1,mn
        vcva(k,yi,xi)=vcva(k,xi,yi)
        !print*,vcva
      end do
    end do
  end do !***********************************


  do k=1,mn
    allocate(am(tnk(k),tnk(k)),mm(tn,tnk(k)),phi(tn,tnk(k)),km(tnk(k),tnk(k)))
    allocate(km_phi(tnk(k),tn))
    am=am2(1:tnk(k),1:tnk(k))
    mm=mm2(1:tn,1:tnk(k))
    phi=matmul(mm,am)   !***************************************
    km=tkm(k,1:tnk(k),1:tnk(k))
    km_phi=matmul(km,transpose(phi))

    do xi=1,tnk(k)    !check
      do zi=1,pedn
        v1=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          !do j=1,rnm(yi)
            !v1=v1+mbin(k,yidx(yi,j),zi)*km_phi(xi,yi)*py(trny-rnm(yi)+j,1)
          !end do
          !v1=v1+eval(zi,k)*km_phi(xi,yi)*py(trny-rnm(yi)+zi,1)
          v1=v1+eval(zi,k)*km_phi(xi,yi)*py(trny-rnm(yi)+yidx2(zi),1)
        end do
        blup_ebv(zi,xi,k)=v1
      end do !zi
    end do !xi
    deallocate(am,mm,phi,km,km_phi)
  end do !k

  !transformed back
  trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    do zi=1,pedn
      do k=1,mn
        v1=0
          do j=1,pedn
            !v1=v1+mbin(k,zi,j)*blup_ebv(j,xi,k)
            v1=v1+mbin(1,zi,j)*blup_ebv(j,xi,k)    !rrm e use the same evec
          end do
        blup_ebv2(zi,xi,k)=v1
      end do !k
    end do !zi
  end do !xi

  open (unit=47,file=trim(fl11)//".fsl",status='unknown')
  write(47,'(a3,2a14,a12,a14)')'#','BETA','SE','CHI','P'

  do xi=1,tfn
    !CDF functions *****************************************************
    p = huge ( p )
    q = huge ( q )
    x = beta(xi,1)**2/xvmx(xi,xi)
    sd = 1.0D+00 !degree of greedom for cdfchi
    which=1  !Calculate P and Q from X, MEAN and SD;
    call cdfchi ( which, p, q, x, sd, status, bound )
    write(47,'(i3,2f14.5,f12.3,e14.5)')xi,beta(xi,1),sqrt(xvmx(xi,xi)),x,q
  end do
  close(47)

  open (unit=43,file=fl11,status='unknown')
  write(43,'(a6,a14,a9,a14,a24)')'trait','random eff.#','RR order','ordered ind#','EBVs'

  !do ii=1,1
    do k=1,mn
      do xi=1,tnk(k)    !check
        do i=1,pedn
          write(43,'(i6,i14,i9,i14,f24.16)')ii,k,xi,i,blup_ebv2(i,xi,k)
        end do
      end do
    end do
  !end do

  close(43)


end if



end subroutine

!end module
