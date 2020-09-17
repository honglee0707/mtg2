
subroutine aireml_m_eig (mbin,eval,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,xmm,nit,conv)        

!***********************************************************************
!aireml_m: multivariate analysis  
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

INTEGER::n,nped,ix,iz,tn,nb,nbc,nbc2,xi,yi,ri,zi,vi,ui,wi,wj,nit
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,yidx1(rn),yidx2(rn)
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1),eval(pedn,1)
!double precision::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1),eval(pedn,1)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::vcva_rsv(mn,tn,tn),vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
double precision::blup_ebv2(pedn,tn,mn),blup_py2(pedn,tn,mn)
double precision::d_beta(tfn,1),d_res(trn,1)
integer::i,j,m,k,l,io
double precision::sum_h2,sum_v(tn),res(trn,1)


double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v3,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,1000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn),xvmy(trn,1) 
double precision::py(trn,1),xb(trn,1),tmp_pm(tn,tn)
double precision::dpy(trn,1),xvmdpy(tfn,1),pdpy(trn,1)

double precision::vmi(pedn,tn,tn),pm(pedn,tn,tn)

double precision::v_tmp(trn),v2_tmp(trn)

!double precision::aim(2+mn*3,2+mn*3),aiminv(2+mn*3,2+mn*3),up(2+mn*3,1),dldv(2+mn*3,1)
!double precision::sdm(2+mn*3,2+mn*3),sdmI(2+mn*3,2+mn*3)
double precision::aim(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
double precision::sdmI(tn+mn*(tn+(tn**2-tn)/2),tn+mn*(tn+(tn**2-tn)/2))
double precision::up(tn+mn*(tn+(tn**2-tn)/2),1),dldv(tn+mn*(tn+(tn**2-tn)/2),1)

!for BLUP ebv
!double precision::vmat(mn,pedn,rn)

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs,t0_cpu,t3_cpu,t4_cpu
!***********************************************************************

character(len=7)::cdum
character(len=64)::fl4,fl5,fl11,fl11r
logical::file_exist

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which

do i=1,pedn
  yidx2(yidx(1,i))=i   !index for y to the fam order
end do

!xm=0;trn=0;tfn=0
!do xi=1,tn
!  tfn=tfn+fnm(xi)
!  do i=1,rnm(xi)
!    trn=trn+1
!    v1=0
!    do j=1,rnm(xi)
!      !v1=v1+mbin(1,yidx(xi,i),yidx(xi,j))*yv(yidx(xi,j),xi)
!      v1=v1+mbin(1,yidx(xi,j),yidx(xi,i))*yv(yidx(xi,j),xi)
!    end do
!    !tyv(trn,1)=yv(yidx(xi,i),xi)
!    tyv(trn,1)=v1                !check

!    do j=1,fnm(xi)
!      v1=0
!      do k=1,rnm(xi)
!        !v1=v1+mbin(1,yidx(xi,i),yidx(xi,k))*xmm(xi,yidx(xi,k),j)
!        v1=v1+mbin(1,yidx(xi,k),yidx(xi,i))*xmm(xi,yidx(xi,k),j)
!      end do
!      !xm(trn,tfn-fnm(xi)+j)=xmm(xi,yidx(xi,i),j)
!      xm(trn,tfn-fnm(xi)+j)=v1       !check
!    end do
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
      !print*,yidx(xi,j),j,yv(yidx(xi,j),xi),yv(j,xi)
      !pause
    end do
    tyv(trn,1)=v1                !check
    !print*,tyv(trn,1)
    !pause
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
    !print*,v1
    !pause
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
  print*,'trait',i,':',rnm(i)
end do
print*,''



  ! y = Xb + Zu + e ---- fitting effects without QTL

  !arbitary starting value (half of the total variance)
  vcve=0;vcva=0
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
      do xi=1,tn
        read(44,*,iostat=io)cdum,vcva(wi,xi,xi)
        if (io.ne.0) exit
      end do
      do xi=1,tn
        do yi=1,xi-1
          read(44,*,iostat=io)cdum,vcva(wi,xi,yi)
          if (io.ne.0) exit
        end do
      end do
    end do
    close(44)
  end if

  itit=1
  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000

  vmi=0;pm=0 !just for allocation

  ! V matrix    V= ZAZ'Va+Z2GZ2'Vpe+IVe*************************
  !do zi=1,200
  !do zi=1,nit
  do zi=1,10000000
    call cpu_time (t3_cpu)
    LKHP=LKH
  !LKH and AIREML (without sparce technique)
  ! V matrix    V = ZAZ'Va + IVe *************************

  call cpu_time (t2_cpu)

  LA=0
  do i=1,pedn     !assuming all ID have record
    do xi=1,tn
      do yi=1,xi
        v1=vcve(xi,yi)
        do wi=1,mn
          !v1=v1+vcva(wi,xi,yi)*eval(i,1)
          v1=v1+vcva(wi,xi,yi)*eval(yidx(xi,i),1)
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
            !print*,xm(pedn*xi-pedn+j,i),vmi(j,xi,yi)
          end do
          !pause
          xvm(i,pedn*yi-pedn+j)=v1
          !print*,v1
          !pause
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

    !do i=1,tfn
    !  v2=0
    !  do j=1,trn
    !    v1=0
    !    do k=1,trn
    !      v1=v1+xvm2(i,k)*vmi(k,j)    !*pm
    !    end do
    !    v2=v2+v1*tyv(j,1)           !*tyv
    !    xvm3(j,i)=v1
    !  end do
    !  beta(i,1)=v2
    !end do
    xvmy=0
    do i=1,trn
      v1=0
      do j=1,tfn
        xvmy(j,1)=xvmy(j,1)+xvm(j,i)*tyv(i,1)
      end do
    end do
   beta=matmul(xvmx,xvmy)  
   !*********************************************
   !print*,'beta:',beta
    call cpu_time (t4_cpu)
    !print*,'test  - time:',real(t4_cpu-t3_cpu) !,LA

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

    call cpu_time (t4_cpu)
    !print*,'test  - time:',real(t4_cpu-t3_cpu) !,LA

    v2=0
    do i=1,trn
      v2=v2+py(i,1)*tyv(i,1)
    end do

    LKH=-0.5*(LA+LG+v2)
    !print*,LA,LG,v2  !ypyf
    call cpu_time (t4_cpu)
    !print*,'test  - time:',real(t4_cpu-t3_cpu) !,LA
    !pause
    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

    if ((LKH>-99999999 .and. LKH<99999999) .and. (LKH>=LKHP)) then
      itit=1
      !print*,'ok',itit,LKHP
    else
      itit=itit+1
      !LKH=-1000000000
      LKH=LKHP
      !print*,'likelihood NaN >> reduce by 0.7^',itit-1
      print*,itit-1,'likelihood nan >> update reduced by the factor'
      !print*,'not ok',itit,LKHP,LKH
      !pause
      goto 111
    end if


    PRINT '(a7,100f12.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f12.4)','Ve',vcve(xi,xi)
    end do
    print*,''
    do wi=1,mn
      do xi=1,tn
        PRINT '(a7,100f12.4)','Va',vcva(wi,xi,xi)
      end do
      do xi=2,tn
        do yi=1,xi-1
          PRINT '(a7,100f12.4)','cov',vcva(wi,xi,yi)
        end do
      end do
      print*,''
    end do
    !pause


    ! AI matrix
    !py=MATMUL(pm,tyv)       !for py

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    nb=tn+(tn**2-tn)/2  !no block for random effects (Vg, cov)
    !Ve (diagonal)  **********************************
    !trnx=0
    !do xi=1,tn
    !  trnx=trnx+rnm(xi)   !total rn upto xi
    !  v2=0
    !  do i=1,rnm(xi)
    !    v1=0
    !    do j=1,rnm(xi)
    !      v1=v1+py(trnx-rnm(xi)+j,1)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
    !    end do
    !    v2=v2+v1*py(trnx-rnm(xi)+i,1)
    !  end do 
    !  aim(xi,xi)=v2
    !  print*,v2
    !end do

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
    end do

    !print*,'Vg'
    !for Vg
    do k=1,mn !no. random effects
      !for Vg (doagonal) ******************************************************* 
      trnx=0
      do xi=1,tn
        trnx=trnx+rnm(xi)
        dpy=0
        !v2=0
        xvmdpy=0    !xvm %*% dpy
        do i=1,rnm(xi)
          !eig simplification
          v10=py(trnx-rnm(xi)+i,1)*eval(yidx(xi,i),1)     !py*pig
          !v_tmp(trnx-rnm(xi)+i)=v10
          dpy(trnx-rnm(xi)+i,1)=v10
          !xvm %*% dpy
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
        !print*,v3
        aim(tn+nb*k-nb+xi,tn+nb*k-nb+xi)=v3

        !do i=1,trn
        !  v10=0
        !  do j=1,rnm(xi)
        !    v10=v10+v_tmp(trnx-rnm(xi)+j)*pm(trnx-rnm(xi)+j,i)                      !*pm
        !  end do
        !  v2_tmp(i)=v10
        !end do

        !do i=1,rnm(xi)
        !  !eig simplification
        !  v1=v2_tmp(trnx-rnm(xi)+i)*eval(yidx(xi,i),1)     !*pig
        !  v2=v2+v1*py(trnx-rnm(xi)+i,1)                                !*py
        !end do
        !print*,"v2",v2
        !aim(tn+nb*k-nb+xi,tn+nb*k-nb+xi)=v2
        !pause

        !for Vg x Ve (off diagonal) !******************************************
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)                                !*py
            v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)                                !*py
          end do
          aim(tn+nb*k-nb+xi,vi)=v2
        end do

        !for Vg(1~[i-1]) x Vg(i)
        trnv=0
        do vi=1,xi-1
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            !eig simplification  >>> check when using multiple GRMs
            !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)   !*pig_tmp
            v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)   !*pig_tmp
            v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py

          end do
          aim(tn+nb*k-nb+xi,tn+nb*k-nb+vi)=v2
        end do


        !off diagonal *********************************************************
        do wi=1,k-1
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              !eig simplification
              !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)     !*pig_tmp
              v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)     !*pig_tmp
              v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
            end do
            aim(tn+nb*k-nb+xi,tn+nb*wi-nb+vi)=v2   !with the xi th Vg for first trait
          end do

          nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
          trnv=rnm(1)
          do vi=2,tn
            trnv=trnv+rnm(vi)   !total rn upto vi
            trnu=0
            do ui=1,(vi-1)
              nbc=nbc+1
              trnu=trnu+rnm(ui)
              v2=0
              do i=1,rnm(vi)
                !eig simplification  >>> check when using diff N for diff trait
                !v1=v2_tmp(trnu-rnm(ui)+i)*eval(yidx(ui,i),1)     !*pig_tmp
                v1=pdpy(trnu-rnm(ui)+i,1)*eval(yidx(ui,i),1)     !*pig_tmp
                v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
              end do
              do i=1,rnm(ui)
                !eig simplification >>> check when using diff M for diff trait
                !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)     !*pig_tmp
                v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)     !*pig_tmp
                v2=v2+v1*py(trnu-rnm(ui)+i,1)                                !*py
              end do
              aim(tn+nb*k-nb+xi,tn+nb*wi-nb+tn+nbc)=v2  !check
            end do

          end do

        end do !wi

      end do  !xi
      !pause


      !dldv Vg *********************************************************
      trnx=0
      do xi=1,tn
        trnx=trnx+rnm(xi)   !total rn upto xi
        tr1=0
        do i=1,rnm(xi)
          !eig simplification
          !v1=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)*eval(yidx(xi,i),1)
          v1=pm(i,xi,xi)*eval(yidx(xi,i),1)
          tr1=tr1+v1
        end do
        tr1=-0.5*tr1

        v2=0
        do i=1,rnm(xi)
          v1=0
          v1=py(trnx-rnm(xi)+i,1)*eval(yidx(xi,i),1)
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 
        dldv(tn+nb*k-nb+xi,1)=v2*0.5+tr1 !************************************
      end do


      !for cov (diagonal) *************************************************
      nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
      trnx=rnm(1)
      do xi=2,tn
        trnx=trnx+rnm(xi)
        trny=0
        do yi=1,xi-1
          nbc=nbc+1
          trny=trny+rnm(yi)
          dpy=0
          v2=0
          xvmdpy=0    !xvm %*% dpy
          do i=1,rnm(xi)
            !eig simplification
            v10=py(trny-rnm(yi)+i,1)*eval(yidx(yi,i),1)     !py*pig
            !v_tmp(trnx-rnm(xi)+i)=v10
            dpy(trnx-rnm(xi)+i,1)=v10
            !xvm %*% dpy
            do j=1,tfn
              xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,trnx-rnm(xi)+i)*dpy(trnx-rnm(xi)+i,1)
            end do
          end do

          do i=1,rnm(yi)
            !eig simplification
            v10=py(trnx-rnm(xi)+i,1)*eval(yidx(xi,i),1)     !py*pig
            !v_tmp(trny-rnm(yi)+i)=v10
            dpy(trny-rnm(yi)+i,1)=v10
            !xvm %*% dpy
            do j=1,tfn
              xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,trny-rnm(xi)+i)*dpy(trny-rnm(yi)+i,1)
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
          !print*,v3
          aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+tn+nbc)=v3     !check


          !do i=1,trn     !check (may not trn)
          !  v10=0
          !  do j=1,rnm(xi)
          !    v10=v10+v_tmp(trnx-rnm(xi)+j)*pm(trnx-rnm(xi)+j,i)              !*pm
          !  end do
          !  do j=1,rnm(yi)
          !    v10=v10+v_tmp(trny-rnm(yi)+j)*pm(trny-rnm(yi)+j,i)              !*pm
          !  end do
          !  v2_tmp(i)=v10
          !end do 


          !do i=1,rnm(xi)
          !  !eig simplification (assum full record without missing)
          !  v1=v2_tmp(trny-rnm(yi)+i)*eval(yidx(yi,i),1)     !*pig
          !  v2=v2+v1*py(trnx-rnm(xi)+i,1)                                !*py
          !end do
          !do i=1,rnm(yi)
          !  !eig simplification
          !  v1=v2_tmp(trnx-rnm(xi)+i)*eval(yidx(xi,i),1)     !*pig
          !  v2=v2+v1*py(trny-rnm(yi)+i,1)                                !*py
          !end do
          !print*,v2
          !aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+tn+nbc)=v2     !check


          !for cov x Ve (off diagonal) ************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)                                !*py
              v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)                                !*py
            end do 
            aim(tn+nb*k-nb+tn+nbc,vi)=v2          !check with Ve1
          end do

          !for cov x Vg **************************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)     !*pig_tmp
              v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)     !*pig_tmp
              v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
            end do 
            aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+vi)=v2  !with the same Vg1 for first trait
          end do

          !for cov x cov *************************************************
          trnv=rnm(1)
          nbc2=0    !no block for cov (vi,ui)
          do vi=2,tn
            trnv=trnv+rnm(vi)
            trnu=0
            do ui=1,(vi-1)
              nbc2=nbc2+1

              if (nbc2<nbc) then  !off diagonal within cova ************************
              trnu=trnu+rnm(ui)
              v2=0
              do i=1,rnm(vi)
                !eig simplification
                !v1=v2_tmp(trnu-rnm(ui)+i)*eval(yidx(ui,i),1)     !*pig_tmp
                v1=pdpy(trnu-rnm(ui)+i,1)*eval(yidx(ui,i),1)     !*pig_tmp
                v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
              end do 
              do i=1,rnm(ui)
                !eig simplification
                !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)     !*pig_tmp
                v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)     !*pig_tmp
                v2=v2+v1*py(trnu-rnm(ui)+i,1)                                !*py
              end do 
              aim(tn+nb*k-nb+tn+nbc,tn+nb*k-nb+tn+nbc2)=v2  !with the wi th cov for first trait

              end if ! *************************************************************
            end do !ui
          end do !vi


          !off diagonal *********************************************************
          do wi=1,k-1
            trnv=0
            do vi=1,tn
              trnv=trnv+rnm(vi)
              v2=0
              do i=1,rnm(vi)
                !eig simplificaiton
                !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)     !*pig_tmp
                v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)     !*pig_tmp
                v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
              end do 
              aim(tn+nb*k-nb+tn+nbc,tn+nb*wi-nb+vi)=v2  !with the wi th Vg for first trait
            end do

            trnv=rnm(1)
            nbc2=0    !no block for cov (vi,ui)
            do vi=2,tn
              trnv=trnv+rnm(vi)
              trnu=0
              do ui=1,(vi-1)
                nbc2=nbc2+1
                trnu=trnu+rnm(ui)
                v2=0
                do i=1,rnm(vi)
                  !eig simlification
                  !v1=v2_tmp(trnu-rnm(ui)+i)*eval(yidx(ui,i),1)     !*pig_tmp
                  v1=pdpy(trnu-rnm(ui)+i,1)*eval(yidx(ui,i),1)     !*pig_tmp
                  v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
                end do 
                do i=1,rnm(ui)
                  !eig simplification
                  !v1=v2_tmp(trnv-rnm(vi)+i)*eval(yidx(vi,i),1)     !*pig_tmp
                  v1=pdpy(trnv-rnm(vi)+i,1)*eval(yidx(vi,i),1)     !*pig_tmp
                  v2=v2+v1*py(trnu-rnm(ui)+i,1)                                !*py
                end do 
                aim(tn+nb*k-nb+tn+nbc,tn+nb*wi-nb+tn+nbc2)=v2  !with the wi th cov for first trait
              end do !ui
            end do !vi
          end do !wi
        end do !yi

      end do !xi **********************************************

      !dldv ***************************************************************
      nbc=0     !no block for cov (xi, yi)
      trnx=rnm(1)
      do xi=2,tn
        trnx=trnx+rnm(xi)
        trny=0
        do yi=1,xi-1
          nbc=nbc+1
          trny=trny+rnm(yi)
          tr1=0
          do i=1,rnm(xi)
            !v10=pm(trnx-rnm(xi)+i,trny-rnm(xi)+i)*eval(yidx(xi,i),1)     !pm*pig
            v10=pm(i,xi,yi)*eval(yidx(xi,i),1)     !pm*pig
            tr1=tr1+v10
          end do
          do i=1,rnm(yi)
            !eig simplification
            !v10=pm(trny-rnm(yi)+i,trnx-rnm(yi)+i)*eval(yidx(yi,i),1)     !pm*pig
            v10=pm(i,yi,xi)*eval(yidx(yi,i),1)     !pm*pig
            tr1=tr1+v10
          end do
          tr1=-0.5*tr1

          v2=0
          do i=1,rnm(xi)
            !eig simplification
            v10=py(trny-rnm(yi)+i,1)*eval(yidx(yi,i),1)     !py*pig
            v2=v2+v10*py(trnx-rnm(xi)+i,1)
          end do
          do i=1,rnm(yi)
            v10=py(trnx-rnm(xi)+i,1)*eval(yidx(xi,i),1)     !py*pig
            v2=v2+v10*py(trny-rnm(yi)+i,1)
          end do
          dldv(tn+nb*k-nb+tn+nbc,1)=v2*0.5+tr1
        end do !yi
      end do !xi


    end do !k

    call cpu_time (t2_cpu)
    print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    print*,''
 
    !do i=1,(2+mn*3)
    do i=1,tn+nb*mn
      aim(i,i)=aim(i,i)/2
      do j=1,i-1
        aim(i,j)=aim(i,j)/2
        aim(j,i)=aim(i,j)
      end do
    end do

    !do i=1,tn+nb*mn
    !  print*,real(aim(i,1:i))*2
    !end do
    !pause

    !call cholesky_inv (aim,(2+mn*3),x1)
    call cholesky_inv (aim,(tn+nb*mn),x1)
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
    !do i=1,tn+nb*mn
    !  print*,real(up(i,1)),real(dldv(i,1))
    !end do
    !print*,''
    !pause

111 continue
    if (itit==1) then
      vcve_rsv=vcve
      vcva_rsv=vcva
      if (zi.ge.nit) then
        print*,'Likelihood not converged >>> may need a longer iterations'
        print*,''
        goto 1000
      end if

    else
      up=up*0.7**(itit-1)
      vcve=vcve_rsv
      vcva=vcva_rsv
    end if

    do xi=1,tn
      vcve(xi,xi)=vcve(xi,xi)+up(xi,1)
    end do
    do wi=1,mn
      do xi=1,tn
        vcva(wi,xi,xi)=vcva(wi,xi,xi)+up(tn+nb*wi-nb+xi,1)
          !constrain out of parameter space
          !if (vcva(wi,xi,xi).lt.0) vcva(wi,xi,xi)=0.001
      end do
      nbc=0
      do xi=2,tn
        do yi=1,xi-1
          nbc=nbc+1
          vcva(wi,xi,yi)=vcva(wi,xi,yi)+up(tn+nb*wi-nb+tn+nbc,1)

          !constrain out of parameter space
          !if (vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))>=1) then
          !  vcva(wi,xi,yi)=sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))-0.001
          !end if

        end do
      end do
    end do

    !print*,ve,ve2,va,va2,cov
    !print*,vcve,vcva

    
  call cpu_time (t0_cpu)
  !print*,'iteration done - time:',real(t0_cpu-t3_cpu)

  END do !zi



  !PRINT*,'LKH not converged na na na'
  !WRITE(41,*)'LKH not converged na na na'
  va=mva;va2=mva2;ve=mve;ve2=mve2;LKH=MLKH;cov=mcov


1000 continue 


  open (UNIT=41,FILE=fl5,STATUS='unknown')

  do xi=1,tn
    write(41,'(a7,100f12.4)')'Ve',vcve(xi,xi),sqrt(sdmI(xi,xi))
  end do
  write(41,*)''

  do xi=1,tn
    sum_v(xi)=vcve(xi,xi)
  end do
  do wi=1,mn
    do xi=1,tn
      write(41,'(a7,100f12.4)')'Va',vcva(wi,xi,xi),sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f12.4)')'cov',vcva(wi,xi,yi),sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
      end do
    end do
    write(41,*)''
  end do



  !total Var(V) *********************
  do xi=1,tn
    vvc(xi)=sdmI(xi,xi)
    do wi=1,mn
      vvc(xi)=vvc(xi)+sdmI(tn+nb*wi-nb+xi,xi)*2
      do wj=1,mn
        vvc(xi)=vvc(xi)+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
      end do
    end do
  end do
  !print*,vvc

  do wi=1,mn
    do xi=1,tn
      !x1=va(i)/sum_v
      x1=vcva(wi,xi,xi)/sum_v(xi)
      !x11=sdmI(2+i*3-2,1) !for the row ************
      x11=sdmI(tn+nb*wi-nb+xi,xi) !for the row ************

      do wj=1,mn
        !x11=x11+sdmI(2+i*3-2,2+j*3-2)
        x11=x11+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
      end do
      !print*,x11

      !SE **********************************
      !x2=x1**2*(sdmI(2+i*3-2,2+i*3-2)/va(i)**2+x10/sum_v**2-2*x11/(va(i)*sum_v))
      x2=x1**2*(sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      write(41,'(a7,100f12.4)')'h2',x1,sqrt(x2)

    end do
    write(41,*)''

    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        !z1=cov(i)/sqrt(va(i)*va2(i))
        z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
        !var(cor) *************************************
        !z2=sdmI(2+i*3-2,2+i*3-2)/(4*va(i)**2)
        z2=sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/(4*vcva(wi,xi,xi)**2)
        !z2=z2+sdmI(2+i*3-1,2+i*3-1)/(4*va2(i)**2)
        z2=z2+sdmI(tn+nb*wi-nb+yi,tn+nb*wi-nb+yi)/(4*vcva(wi,yi,yi)**2)
        !z2=z2+sdmI(2+i*3,2+i*3)/(cov(i)**2)
        z2=z2+sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+tn+nbc)/(vcva(wi,xi,yi)**2)
        !z2=z2+2*sdmI(2+i*3-2,2+i*3-1)/(4*va(i)*va2(i))
        z2=z2+2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+yi)/(4*vcva(wi,xi,xi)*vcva(wi,yi,yi))
        !z2=z2-2*sdmI(2+i*3-2,2+i*3)/(2*va(i)*cov(i))
        z2=z2-2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+tn+nbc)/(2*vcva(wi,xi,xi)*vcva(wi,xi,yi))
        !z2=z2-2*sdmI(2+i*3,2+i*3-1)/(2*cov(i)*va2(i))
        z2=z2-2*sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+yi)/(2*vcva(wi,xi,yi)*vcva(wi,yi,yi))
        z2=z2*(z1**2)
        write(41,'(a7,100f12.4)')'cor',z1,sqrt(z2)
      end do
    end do
    write(41,*)''
  end do

  write(41,*)''
  write(41,'(a7,100f12.4)')'LKH',LKH

  write(41,*)''
  do i=1,(tn+nb*mn)
    !write(41,*)real(sdmI(i,1:i))
    !write(41,'(100f12.5)')real(sdmI(i,1:i))
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

  trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    do zi=1,pedn
      do k=1,mn
        v1=0;v2=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          !do j=1,rnm(yi)
          !  v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
          !  if (yidx(yi,j)==zi) then
          !    v2=v2+vcva(k,xi,yi)*py(trny-rnm(yi)+j,1)
          !  end if
          !end do
          !v1=v1+eval(zi,1)*vcva(k,xi,yi)*py(yi*pedn-pedn+zi,1)
          !v2=v2+vcva(k,xi,yi)*py(yi*pedn-pedn+zi,1)
          v1=v1+eval(zi,1)*vcva(k,xi,yi)*py(yi*pedn-pedn+yidx2(zi),1)
          v2=v2+vcva(k,xi,yi)*py(yi*pedn-pedn+yidx2(zi),1)
        end do
        blup_ebv(zi,xi,k)=v1
        blup_py(zi,xi,k)=v2
      end do !k
    end do !zi
  end do !xi

  !transformed back
  trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    do zi=1,pedn
      do k=1,mn
        v1=0;v2=0
          do j=1,pedn
            v1=v1+mbin(k,zi,j)*blup_ebv(j,xi,k)
            v2=v2+mbin(k,zi,j)*blup_py(j,xi,k)
          end do
        blup_ebv2(zi,xi,k)=v1
        blup_py2(zi,xi,k)=v2
      end do !k
    end do !zi
  end do !xi


  open (unit=45,file=trim(fl11)//".py",status='unknown')
  do xi=1,tn
    do i=1,pedn
      write(45,'(i3,100f24.16)')xi,blup_py2(i,xi,:)
    end do
  end do
  close(45)

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
  !write(43,'(100f24.16)')beta
  write(43,'(a5,a26)')'trait','EBVs ...'
  do xi=1,tn
    do i=1,pedn
      write(43,'(i3,100f24.16)')xi,blup_ebv2(i,xi,:)
    end do
  end do
  close(43)

end if

if (fl11r.ne.'null') then

  print*,'With -eig, -bvr (reliability) is not supported'
  print*,'Use -g instead of -eig '
  pause

end if


end subroutine

!end module
