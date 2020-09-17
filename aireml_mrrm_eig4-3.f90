
subroutine aireml_rrm_eig_m2 (mbin,eval,yidx,yv,rn,rnm,trn,pedn,ttn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl12,xmm,nit,conv)        

!***********************************************************************
!aireml_rrm_eig_m: multivariate random regression analysis  
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

INTEGER::n,nped,ix,iz,nb,nb3,nbc,nbc2,xi,xi2,yi,yi2,ri,zi,wj
INTEGER::ttn,tn,sub_tn(ttn),wi,wi2,vi,vi2,ui,ui2
integer::m_tnk(ttn,mn)  !,ctnk(ttn,ttn,mn)    !no. polynomial components order
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,yidx1(rn),yidx2(rn),nit
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1),eval(pedn,mn)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),pig_vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
double precision::blup_ebv2(pedn,tn,mn) !,cor
integer::i,i2,j,j2,m,k,k2,l,io
double precision::sum_h2,sum_v(tn),flegendre
double precision::d_beta(tfn,1),d_res(trn,1),res(trn,1)

double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v3,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,1000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn),xvmy(trn,1)
double precision::py(trn,1),xb(trn,1)
double precision::pig(pedn,tn,tn),pig_tmp(pedn,tn,tn)

double precision::v_tmp(trn),v2_tmp(trn),tmp_pm(tn,tn)
double precision::dpy(trn,1),xvmdpy(tfn,1),pdpy(trn,1)

double precision::vmi(pedn,tn,tn),pm(pedn,tn,tn)

!for BLUP ebv

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

double precision::m_tkm(mn,ttn,20,20),m_tkm_rsv(mn,ttn,20,20)
double precision::ctkm(mn,ttn,ttn,20,20),ctkm_rsv(mn,ttn,ttn,20,20)
double precision,allocatable::km(:,:),km2(:,:),phi(:,:),km2_2(:,:),phi_2(:,:)
double precision,allocatable::phi2(:,:),tmp_vcva(:,:),tmp_vcva_2(:,:),phi_1(:,:)
double precision::wk(tn),m_wk(ttn,tn),xk(tn)
double precision,allocatable::tmp_wk(:),tmp_wk2(:),tmp_wk_2(:),tmp_wk_1(:)
double precision,allocatable::km_phi(:,:),km_phi2(:,:),gr_km_phi(:,:)

double precision,allocatable::aim(:,:)
double precision,allocatable::sdmI(:,:)
double precision,allocatable::up(:,:),dldv(:,:)

integer::idum,idum2,ii,ii2,jj,jj2,kk,kk2,ll,mm,mm2,oo,ii3,ii12,jj12,ii22,jj22

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which


open (UNIT=46,FILE=fl12,STATUS='old')
read(46,*)(sub_tn(i),i=1,ttn)
do i=1,ttn                      !# traits for random regression coeff. (rrc)
  read(46,*)m_tnk(i,:)
  read(46,*)m_wk(i,1:sub_tn(i))
end do
close(46)

do i=1,pedn
  yidx2(yidx(1,i))=i   !index for y to the fam order
end do

!do i=2,ttn
!  if (sub_tn(i).ne.sub_tn(i-1)) then
!    print*,'Not the same # sites across traits, check'
!    pause
!  end if
!end do

ii3=0
do ii=2,ttn
  do jj=1,(ii-1)
    do ii2=1,sub_tn(ii)
      do jj2=1,sub_tn(jj)
        if (m_wk(jj,jj2)==m_wk(ii,ii2)) then
          ii3=ii3+1
          !print*,m_wk(ii,ii2),m_wk(jj,jj2),ii3,ii,ii2,jj,jj2
        end if
      end do
    end do
  end do
end do


idum=tn+ii3
!idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
do ll=1,mn
  do kk=1,ttn
    idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
  end do
  do kk=2,ttn
    do mm=1,kk-1
      idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
    end do
  end do
end do

!print*,idum,tn+mn*(tnk(1)+(tnk(1)**2-tnk(1))/2)
allocate(aim(idum,idum),sdmI(idum,idum),up(idum,1),dldv(idum,1))


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

  ! y = Xb + Zu + e ---- Linear mixed model
  !arbitary starting value (half of the total variance)
  vcve=0;vcva=0;m_tkm=0;ctkm=0   !;km=0
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

  nbc2=0
  !do ii2=1,sub_tn(1)            !because the same # across traits
  do ii2=2,ttn
    do jj2=1,(ii2-1)
      do ii22=1,sub_tn(ii2)
        do jj22=1,sub_tn(jj2)
          if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
            nbc2=nbc2+1
            ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
            vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

          !read(44,*,iostat=io)cdum,vcve(jj*sub_tn(1)-sub_tn(1)+xi,kk*sub_tn(1)-sub_tn(1)+xi)
          read(44,*,iostat=io)cdum,vcve(ui,vi)
          !print*,vcve(ui,vi)
        end if

        end do
      end do
    end do
  end do 
         
    do wi=1,mn
      do ii=1,ttn
        allocate(km(m_tnk(ii,wi),m_tnk(ii,wi)))

        do xi=1,m_tnk(ii,wi)
          read(44,*,iostat=io)cdum,km(xi,xi)
          if (io.ne.0) exit
        end do
        do xi=1,m_tnk(ii,wi)
          do yi=1,xi-1
            read(44,*,iostat=io)cdum,km(xi,yi)
            km(yi,xi)=km(xi,yi)
            if (io.ne.0) exit
          end do
        end do
        m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))=km
        deallocate(km)
      end do
      do ii=2,ttn
        do jj=1,ii-1
          allocate(km(m_tnk(ii,wi),m_tnk(jj,wi)))
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              read(44,*,iostat=io)cdum,km(xi,yi)
              !print*,km(xi,yi)
              if (io.ne.0) exit
            end do
          end do
          ctkm(wi,ii,jj,1:m_tnk(ii,wi),1:m_tnk(jj,wi))=km
          deallocate(km)
        end do
      end do

    end do !wi
    close(44)
  end if

  !Legendre polynomial function ****************************

  !polynomial components not same for all multiple random effects
  do wi=1,mn
    do ii=1,ttn
      allocate(phi(sub_tn(ii),m_tnk(ii,wi)),km(m_tnk(ii,wi),m_tnk(ii,wi)))
      allocate(tmp_wk(sub_tn(ii)),tmp_vcva(sub_tn(ii),sub_tn(ii)))

      tmp_wk=m_wk(ii,1:sub_tn(ii))
      call legendre (m_tnk(ii,wi),tmp_wk,sub_tn(ii),phi)

      km=m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))

      tmp_vcva=matmul(matmul(phi,km),transpose(phi))
      mm=sum(sub_tn(1:ii))-sub_tn(ii)
      do jj=1,sub_tn(ii)
        do kk=1,sub_tn(ii)
          vcva(wi,mm+jj,mm+kk)=tmp_vcva(jj,kk)
        end do
      end do
      deallocate(phi,km,tmp_wk,tmp_vcva)
    end do

    do ii=2,ttn
      do jj=1,ii-1
        allocate(km(m_tnk(ii,wi),m_tnk(jj,wi)),tmp_vcva(sub_tn(ii),sub_tn(jj)))
        allocate(phi(sub_tn(ii),m_tnk(ii,wi)),tmp_wk(sub_tn(ii)))
        tmp_wk=m_wk(ii,1:sub_tn(ii))
        call legendre (m_tnk(ii,wi),tmp_wk,sub_tn(ii),phi)

        allocate(phi2(sub_tn(jj),m_tnk(jj,wi)),tmp_wk2(sub_tn(jj)))
        tmp_wk2=m_wk(jj,1:sub_tn(jj))
        call legendre (m_tnk(jj,wi),tmp_wk2,sub_tn(jj),phi2)

        km=ctkm(wi,ii,jj,1:m_tnk(ii,wi),1:m_tnk(jj,wi))

        tmp_vcva=matmul(matmul(phi,km),transpose(phi2))
        mm=sum(sub_tn(1:ii))-sub_tn(ii)
        mm2=sum(sub_tn(1:jj))-sub_tn(jj)
        do ll=1,sub_tn(ii)
          do kk=1,sub_tn(jj)
            vcva(wi,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
            vcva(wi,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
          end do
        end do
        deallocate(phi,phi2,km,tmp_wk,tmp_wk2,tmp_vcva)
      end do
    end do
  end do  !wi

  itit=1   !iteration in the iteration

  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000

  vmi=0;pm=0
  !print*,'iteration start'
  !do zi=1,nit
  do zi=1,10000000
    LKHP=LKH

  call cpu_time (t2_cpu)
  !polynomial components not same for all multiple random effects
  do wi=1,mn
    do ii=1,ttn
      allocate(phi(sub_tn(ii),m_tnk(ii,wi)),km(m_tnk(ii,wi),m_tnk(ii,wi)))
      allocate(tmp_wk(sub_tn(ii)),tmp_vcva(sub_tn(ii),sub_tn(ii)))

      tmp_wk=m_wk(ii,1:sub_tn(ii))
      call legendre (m_tnk(ii,wi),tmp_wk,sub_tn(ii),phi)

      km=m_tkm(wi,ii,1:m_tnk(ii,wi),1:m_tnk(ii,wi))

      tmp_vcva=matmul(matmul(phi,km),transpose(phi))
      mm=sum(sub_tn(1:ii))-sub_tn(ii)
      do jj=1,sub_tn(ii)
        do kk=1,sub_tn(ii)
          vcva(wi,mm+jj,mm+kk)=tmp_vcva(jj,kk)
          !print*,mm,jj,kk,mm+jj,mm+kk
        end do
      end do
      !pause
      deallocate(phi,km,tmp_wk,tmp_vcva)
    end do

    do ii=2,ttn
      do jj=1,ii-1
        allocate(km(m_tnk(ii,wi),m_tnk(jj,wi)),tmp_vcva(sub_tn(ii),sub_tn(jj)))
        allocate(phi(sub_tn(ii),m_tnk(ii,wi)),tmp_wk(sub_tn(ii)))
        tmp_wk=m_wk(ii,1:sub_tn(ii))
        call legendre (m_tnk(ii,wi),tmp_wk,sub_tn(ii),phi)

        allocate(phi2(sub_tn(jj),m_tnk(jj,wi)),tmp_wk2(sub_tn(jj)))
        tmp_wk2=m_wk(jj,1:sub_tn(jj))
        call legendre (m_tnk(jj,wi),tmp_wk2,sub_tn(jj),phi2)

        km=ctkm(wi,ii,jj,1:m_tnk(ii,wi),1:m_tnk(jj,wi))

        tmp_vcva=matmul(matmul(phi,km),transpose(phi2))
        mm=sum(sub_tn(1:ii))-sub_tn(ii)
        mm2=sum(sub_tn(1:jj))-sub_tn(jj)
        do ll=1,sub_tn(ii)
          do kk=1,sub_tn(jj)
            vcva(wi,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
            vcva(wi,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
            !print*,mm,mm2,ll,kk,mm+ll,mm2+kk
          end do
        end do
        !pause
        deallocate(phi,phi2,km,tmp_wk,tmp_wk2,tmp_vcva)
      end do
    end do

    !do i=1,tn
    !  print*,vcva(wi,i,:)
    !end do
    !pause

  end do  !wi


  LA=0
  do i=1,pedn     !assuming all ID have record
    do xi=1,tn
      do yi=1,xi
        v1=vcve(xi,yi)
        do wi=1,mn
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
        vmi(i,xi,yi)=tmp_pm(xi,yi)
        vmi(i,yi,xi)=tmp_pm(xi,yi)
      end do
    end do
  end do !i

  !call cpu_time (t2_cpu)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA

    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)
    do i=1,tfn
      do yi=1,tn
        do j=1,pedn
          v1=0
          do xi=1,tn
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
          pm(i,yi,xi)=vmi(i,yi,xi)-v1
        end do
      end do
    end do

    !ypy estimation
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
      print*,itit-1,'likelihood nan >> update reduced by the factor'
      goto 111
    end if

    PRINT '(a7,100f12.4)','LKH',LKH
    do xi=1,tn
      PRINT '(a7,100f12.4)','Ve',vcve(xi,xi)
    end do

  nbc2=0
  !do ii2=1,sub_tn(1)            !because the same # across traits
  do ii2=2,ttn
    do jj2=1,(ii2-1)
      do ii22=1,sub_tn(ii2)
        do jj22=1,sub_tn(jj2)
          if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
            nbc2=nbc2+1
            ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
            vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

          !PRINT '(a7,100f12.4)','cove',vcve(jj*sub_tn(1)-sub_tn(1)+xi,kk*sub_tn(1)-sub_tn(1)+xi)
          PRINT '(a7,100f12.4)','cove',vcve(ui,vi)
          end if
        end do
      end do
    end do
  end do

    print*,''
    do wi=1,mn
      do ii=1,ttn
        do xi=1,m_tnk(ii,wi)
          PRINT '(a7,100f12.4)','Vk',m_tkm(wi,ii,xi,xi)
        end do
        do xi=2,m_tnk(ii,wi)
          do yi=1,xi-1
            PRINT '(a7,100f12.4)','cov',m_tkm(wi,ii,xi,yi)
          end do
        end do
        print*,''
      end do
      do ii=2,ttn
        do jj=1,ii-1
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              PRINT '(a7,100f12.4)','CVk',ctkm(wi,ii,jj,xi,yi)
            end do
          end do
          print*,''
        end do
      end do
      print*,''
    end do
    !pause


    ! AI matrix
    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

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
        tr1=tr1+pm(i,xi,xi)
      end do
      tr1=-0.5*tr1
      v2=0
      do i=1,rnm(xi)
        v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
      end do
      dldv(xi,1)=v2*(0.5)+tr1
      !print*,xi,dldv(xi,1)
    end do

    !cove (diagonal) for residual covariance
    nbc=0
    !do ii=1,sub_tn(1)            !because the same # across traits
    do ii=2,ttn
      do jj=1,(ii-1)
        do ii12=1,sub_tn(ii)
          do jj12=1,sub_tn(jj)
            if (m_wk(ii,ii12)==m_wk(jj,jj12)) then
              nbc=nbc+1
              xi=sum(sub_tn(1:ii))-sub_tn(ii)+ii12
              yi=sum(sub_tn(1:jj))-sub_tn(jj)+jj12

          trnx=0
          do i=1,xi
            trnx=trnx+rnm(i)
          end do
          trny=0
          do i=1,yi
            trny=trny+rnm(i)
          end do

          dpy=0
          xvmdpy=0   !xvm %*% dpy
          v2=0
          do i=1,rnm(xi)
            !eig simplification
            v10=py(trny-rnm(yi)+i,1)     !py*pig
            dpy(trnx-rnm(xi)+i,1)=v10
            !xvm %*% dpy
            do j=1,tfn
              xvmdpy(j,1)=xvmdpy(j,1)+xvm(j,trnx-rnm(xi)+i)*dpy(trnx-rnm(xi)+i,1)
            end do
          end do
          do i=1,rnm(yi)
            !eig simplification
            v10=py(trnx-rnm(xi)+i,1)     !py*pig
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
                v2=v2+vmi(i,ui,vi)*d_res(vi*pedn-pedn+i,1)
              end do
              pdpy(ui*pedn-pedn+i,1)=v2
              v3=v3+dpy(ui*pedn-pedn+i,1)*v2
            end do
          end do
          aim(tn+nbc,tn+nbc)=v3     !check
          !print*,tn+nbc,v3

          ! for cove x ve (off diagonal) ********************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              !v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)   !*py
              v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)   !*py
            end do
            aim(tn+nbc,vi)=v2          !check with Ve1
            !print*,tn+nbc,vi,v2
          end do

          !for cove x cove (off diagonal) **********************************

          nbc2=0
              !do ii2=1,sub_tn(1)            !because the same # across traits
              do ii2=2,ttn
                do jj2=1,(ii2-1)
                  do ii22=1,sub_tn(ii2)
                    do jj22=1,sub_tn(jj2)
                      if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
                        nbc2=nbc2+1
                        ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
                        vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

                if (nbc2<nbc) then  !off diagonal within cove ***************
                  trnu=0
                  do i=1,ui
                    trnu=trnu+rnm(i)
                  end do
                  trnv=0
                  do i=1,vi
                    trnv=trnv+rnm(i)
                  end do

                  v2=0
                  do i=1,rnm(vi)
                    v1=pdpy(trnu-rnm(ui)+i,1)     !*pdpy * d 
                    v2=v2+v1*py(trnv-rnm(vi)+i,1) !*py
                  end do
                  do i=1,rnm(ui)
                    v1=pdpy(trnv-rnm(vi)+i,1)     !*pdpy * d
                    v2=v2+v1*py(trnu-rnm(ui)+i,1)
                  end do
                  aim(tn+nbc,tn+nbc2)=v2

                end if

                      end if
                    end do !jj22
                  end do !ii22
                end do !jj2
              end do !ii2

          !dldv for cove !***********************************************
          tr1=0
          do i=1,rnm(xi)
            v10=pm(i,xi,yi)     !pm*pig
            tr1=tr1+v10
          end do
          do i=1,rnm(yi)
            v10=pm(i,yi,xi)     !pm*pig
            tr1=tr1+v10
          end do
          tr1=-0.5*tr1

          v2=0
          do i=1,rnm(xi)
            v2=v2+py(trny-rnm(yi)+i,1)*py(trnx-rnm(xi)+i,1)
          end do
          do i=1,rnm(yi)
            v10=py(trnx-rnm(xi)+i,1) !py*pig
            v2=v2+py(trnx-rnm(xi)+i,1)*py(trny-rnm(yi)+i,1)
          end do
          dldv(tn+nbc,1)=v2*0.5+tr1
          !print*,v2,tr1

            end if
          end do !jj12
        end do !ii12
      end do !jj
    end do !ii


    !for Vk
    do k=1,mn !no. random effects
      !for Vk (doagonal) ******************************************************* 
      do ii=1,ttn
        !print*,ii
        allocate(phi(sub_tn(ii),m_tnk(ii,k)),km2(m_tnk(ii,k),m_tnk(ii,k)))
        allocate(tmp_wk(sub_tn(ii)),tmp_vcva(sub_tn(ii),sub_tn(ii)))

        tmp_wk=m_wk(ii,1:sub_tn(ii))
        call legendre (m_tnk(ii,k),tmp_wk,sub_tn(ii),phi)

        do xi=1,m_tnk(ii,k)   !check

          km2=0;km2(xi,xi)=1;pig_vcva=0
          tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
          mm=sum(sub_tn(1:ii))-sub_tn(ii)
          do jj=1,sub_tn(ii)
            do kk=1,sub_tn(ii)
              pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
            end do
          end do

          pig=0
          do ui=1,tn
            do i=1,pedn
              v1=0
              do yi=1,tn
                !v2=pig_vcva(k,ui,yi)*eval(i,k)
                v2=pig_vcva(k,ui,yi)*eval(yidx(ui,i),k)
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

          idum=tn+ii3
          !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
          do ll=1,mn
            do kk=1,ttn
              if (ll==k .and. kk==ii) goto 50
              idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
            end do
            do kk=2,ttn
              do mm=1,kk-1
                !idum=idum+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
              end do
            end do
          end do
50        continue
          aim(idum+xi,idum+xi)=v3
          !print*,'Vk',idum+xi,v3 !'tmpm:',tmpm(1,1)/2
          !pause

          !for Vk x Ve (off diagonal) !******************************************
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)     !*dpy
            end do
            aim(idum+xi,vi)=v2
          end do

          !for Vk x cove ********************************************
          nbc2=0
          !do ii2=1,sub_tn(1)            !because the same # across traits
          do ii2=2,ttn
            do jj2=1,(ii2-1)
              do ii22=1,sub_tn(ii2)
                do jj22=1,sub_tn(jj2)
                  if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
                    nbc2=nbc2+1
                    ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
                    vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

                  trnu=0
                  do i=1,ui
                    trnu=trnu+rnm(i)
                  end do
                  trnv=0
                  do i=1,vi
                    trnv=trnv+rnm(i)
                  end do
  
                  v2=0
                  do i=1,rnm(vi)
                    v2=v2+pdpy(trnu-rnm(ui)+i,1)*py(trnv-rnm(vi)+i,1) !*pdpy * pdy (same n for each trait)
                  end do
                  do i=1,rnm(ui)
                    v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnu-rnm(ui)+i,1) !*pdpy * pdy (same n for each trait)
                  end do
                  aim(idum+xi,tn+nbc2)=v2

                  end if
                end do !jj22 
              end do !ii22
            end do !jj2
          end do !ii2

          !for Vg(1~[i-1]) x Vg(i)
          do vi=1,xi-1
            km2=0;km2(vi,vi)=1;pig_vcva=0
            tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
            mm=sum(sub_tn(1:ii))-sub_tn(ii)
            do jj=1,sub_tn(ii)
              do kk=1,sub_tn(ii)
                pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
              end do
            end do

            pig_tmp=0
            do ui=1,tn
              do i=1,pedn
                v1=0
                do yi=1,tn
                  !v2=pig_vcva(k,ui,yi)*eval(i,k)
                  v2=pig_vcva(k,ui,yi)*eval(yidx(ui,i),k)
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
            !print*,'Vk x Vk',idum+xi,idum+vi
          end do

          ! off-diagonal between other and same t (and in other random eff.) for Vk
          do wi=1,k
            oo=ttn                   !for other random effect
            if (wi==k) oo=ii-1       !for within random effect
            do ii2=1,oo

              allocate(phi_2(sub_tn(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(ii2,wi)))
              allocate(tmp_wk_2(sub_tn(ii2)),tmp_vcva_2(sub_tn(ii2),sub_tn(ii2)))

              tmp_wk_2=m_wk(ii2,1:sub_tn(ii2))
              call legendre (m_tnk(ii2,wi),tmp_wk_2,sub_tn(ii2),phi_2)

              do vi=1,m_tnk(ii2,wi)   !check

                km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
                tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                do jj=1,sub_tn(ii2)
                  do kk=1,sub_tn(ii2)
                    pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                  end do
                end do

                pig_tmp=0
                do ui=1,tn
                  do i=1,pedn
                    v1=0
                    do yi=1,tn
                      !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
                      v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
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

                idum2=tn+ii3
                !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
                do ll=1,mn
                  do kk=1,ttn
                    if (ll==wi .and. kk==ii2) goto 52
                    idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                  end do
                  do kk=2,ttn
                    do mm=1,kk-1
                      !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                      idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                    end do
                  end do
                end do
52              continue

                aim(idum+xi,idum2+vi)=v2  !check
                !print*,'Vk x other t',idum+xi,idum2+vi
              end do ! vi

              nbc=0  !for cov for other trait, i.e. 1: (2,1),2: (3,1), 3:(3,2),... 
              do vi2=2,m_tnk(ii2,wi)
                do ui2=1,(vi2-1)
                  nbc=nbc+1

                  km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                  tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                  mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                  do jj=1,sub_tn(ii2)
                    do kk=1,sub_tn(ii2)
                      pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                    end do
                  end do

                  pig_tmp=0
                  do ui=1,tn
                    do i=1,pedn
                      v1=0
                      do yi=1,tn
                        !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
                        v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
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
                  aim(idum+xi,idum2+m_tnk(ii2,wi)+nbc)=v2  !check
                  !print*,'Vk x other t cov',idum+xi,idum2+m_tnk(ii2,wi)+nbc,v2
                end do !ui2
              end do !vi2

              deallocate(phi_2,km2_2,tmp_wk_2,tmp_vcva_2)
            end do ! ii2

            !for Vk x CVk   *************************************
            if (wi<k) then

              do ii2=2,ttn
                do jj2=1,ii2-1

                  allocate(phi_1(sub_tn(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(jj2,wi)))
                  allocate(tmp_wk_1(sub_tn(ii2)),tmp_vcva_2(sub_tn(ii2),sub_tn(jj2)))
 
                  tmp_wk_1=m_wk(ii2,1:sub_tn(ii2))
                  call legendre (m_tnk(ii2,wi),tmp_wk_1,sub_tn(ii2),phi_1)
 
                  allocate(phi_2(sub_tn(jj2),m_tnk(jj2,wi)),tmp_wk_2(sub_tn(jj2)))
                  tmp_wk_2=m_wk(jj2,1:sub_tn(jj2))
                  call legendre (m_tnk(jj2,wi),tmp_wk_2,sub_tn(jj2),phi_2)

                  nbc2=0
                  do vi=1,m_tnk(ii2,wi)
                    do vi2=1,m_tnk(jj2,wi)
                      nbc2=nbc2+1
                      km2_2=0;km2_2(vi,vi2)=1;pig_vcva=0
                      tmp_vcva_2=matmul(matmul(phi_1,km2_2),transpose(phi_2))
                      mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                      mm2=sum(sub_tn(1:jj2))-sub_tn(jj2)
                      do ll=1,sub_tn(ii2)
                        do kk=1,sub_tn(jj2)
                          pig_vcva(wi,mm+ll,mm2+kk)=tmp_vcva_2(ll,kk)
                          pig_vcva(wi,mm2+kk,mm+ll)=tmp_vcva_2(ll,kk)
                        end do
                      end do

                      pig_tmp=0
                      do ui=1,tn
                        do i=1,pedn
                          v1=0
                          do yi=1,tn
                            !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
                            v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
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

                      idum2=tn+ii3
                      !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
                      do ll=1,mn
                        do kk=1,ttn
                          idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                        end do
                        do kk=2,ttn
                          do mm=1,kk-1
                            if (ll==wi .and. kk==ii2 .and. mm==jj2) goto 60
                            !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                            idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                          end do
                        end do
                      end do
60                    continue

                      aim(idum+xi,idum2+nbc2)=v2  !check
                      !print*,'Vk x CVk in prv',idum+xi,idum2+nbc2
                    end do ! vi2
                  end do ! vi

                  deallocate(phi_1,phi_2,km2_2)
                  deallocate(tmp_wk_1,tmp_wk_2,tmp_vcva_2)

                end do ! jj2
              end do ! ii2

            end if ! wk<k

          end do !wi

          !dldv Vk *********************************************************
          tr1=0
          do i=1,pedn
            do ui=1,tn
              do yi=1,tn
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
          !print*,idum+xi,dldv(idum+xi,1)
        end do !xi


        !for cov (diagonal) *************************************************
        nbc=0  !for no. block for covariance, i.e. 1: (2,1), 2: (3,1), 3:(3,2), ... 
        do xi=2,m_tnk(ii,k)
          do yi=1,xi-1
            nbc=nbc+1
            km2=0;km2(xi,yi)=1;km2(yi,xi)=1;pig_vcva=0
            tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
            mm=sum(sub_tn(1:ii))-sub_tn(ii)
            do jj=1,sub_tn(ii)
              do kk=1,sub_tn(ii)
                pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
              end do
            end do

            pig=0
            do ui=1,tn
              do i=1,pedn
                v1=0
                do vi=1,tn
                  !v2=pig_vcva(k,ui,vi)*eval(i,k)
                  v2=pig_vcva(k,ui,vi)*eval(yidx(ui,i),k)
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

            idum=tn+ii3
            !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
            do ll=1,mn
              do kk=1,ttn
                if (ll==k .and. kk==ii) goto 53
                idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
              end do
              do kk=2,ttn
                do mm=1,kk-1
                  idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
                end do
              end do
            end do
53          continue
            aim(idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc)=v3     !check
            !print*,'cov',idum+m_tnk(ii,k)+nbc,v3

            !for cov x Ve (off diagonal) ************************************
            trnv=0
            do vi=1,tn
              trnv=trnv+rnm(vi)
              v2=0
              do i=1,rnm(vi)
                v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)     !*dpy
              end do
              aim(idum+m_tnk(ii,k)+nbc,vi)=v2          !check with Ve1
            end do

            !for cov x cove ********************************************
            nbc2=0
            !do ii2=1,sub_tn(1)            !because the same # across traits
            do ii2=2,ttn
              do jj2=1,(ii2-1)
                do ii22=1,sub_tn(ii2)
                  do jj22=1,sub_tn(jj2)
                    if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
                      nbc2=nbc2+1
                      ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
                      vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

                  trnu=0
                  do i=1,ui
                    trnu=trnu+rnm(i)
                  end do
                  trnv=0
                  do i=1,vi
                    trnv=trnv+rnm(i)
                  end do

                  v2=0
                  do i=1,rnm(vi)
                    v2=v2+pdpy(trnu-rnm(ui)+i,1)*py(trnv-rnm(vi)+i,1) !*pdpy * pdy (same n for each trait)
                  end do
                  do i=1,rnm(ui)
                    v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnu-rnm(ui)+i,1) !*pdpy * pdy (same n for each trait)
                  end do
                  aim(idum+m_tnk(ii,k)+nbc,tn+nbc2)=v2

                    end if
                  end do !jj22
                end do !ii22
              end do !jj2
            end do !ii2

            !for cov x Vg **************************************************
            do vi=1,m_tnk(ii,k)

              km2=0;km2(vi,vi)=1;pig_vcva=0
              tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
              mm=sum(sub_tn(1:ii))-sub_tn(ii)
              do jj=1,sub_tn(ii)
                do kk=1,sub_tn(ii)
                  pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
                end do
              end do

              pig_tmp=0
              do ui=1,tn
                do i=1,pedn
                  v1=0
                  do wi=1,tn
                    !v2=pig_vcva(k,ui,wi)*eval(i,k)
                    v2=pig_vcva(k,ui,wi)*eval(yidx(ui,i),k)
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

              aim(idum+m_tnk(ii,k)+nbc,idum+vi)=v2 !with the same Vg1 for first trait
              !print*,'cov x Vk',idum+m_tnk(ii,k)+nbc,idum+vi
            end do

            !for cov x cov *************************************************
            trnv=rnm(1)
            nbc2=0    !no block for cov (vi,ui)
            do vi=2,m_tnk(ii,k)
              trnv=trnv+rnm(vi)
              trnu=0
              do ui=1,(vi-1)
                nbc2=nbc2+1

                if (nbc2<nbc) then  !off diagonal within cova **********************
                  trnu=trnu+rnm(ui)

                  km2=0;km2(vi,ui)=1;km2(ui,vi)=1;pig_vcva=0
                  tmp_vcva=matmul(matmul(phi,km2),transpose(phi))
                  mm=sum(sub_tn(1:ii))-sub_tn(ii)
                  do jj=1,sub_tn(ii)
                    do kk=1,sub_tn(ii)
                      pig_vcva(k,mm+jj,mm+kk)=tmp_vcva(jj,kk)
                    end do
                  end do

                  pig_tmp=0
                  do wj=1,tn
                    do i=1,pedn
                      v1=0
                      do wi=1,tn
                        !v2=pig_vcva(k,wj,wi)*eval(i,k)
                        v2=pig_vcva(k,wj,wi)*eval(yidx(wj,i),k)
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

                  aim(idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc2)=v2  !with the wi th cov for first trait
                  !print*,'cov x cov',idum+m_tnk(ii,k)+nbc,idum+m_tnk(ii,k)+nbc2
                end if ! *************************************************************
              end do !ui
            end do !vi


            ! off-diagonal between other traits for cov
            do wi=1,k
              oo=ttn                   !for other random effect
              if (wi==k) oo=ii-1       !for within random effect

              do ii2=1,oo
                allocate(phi_2(sub_tn(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(ii2,wi)))
                allocate(tmp_wk_2(sub_tn(ii2)),tmp_vcva_2(sub_tn(ii2),sub_tn(ii2)))

                tmp_wk_2=m_wk(ii2,1:sub_tn(ii2))
                call legendre (m_tnk(ii2,wi),tmp_wk_2,sub_tn(ii2),phi_2)

                do vi=1,m_tnk(ii2,wi)   !check

                  km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
                  tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                  mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                  do jj=1,sub_tn(ii2)
                    do kk=1,sub_tn(ii2)
                      pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                    end do
                  end do

                  pig_tmp=0
                  do ui=1,tn
                    do i=1,pedn
                      v1=0
                      do vi2=1,tn
                        !v2=pig_vcva(wi,ui,vi2)*eval(i,wi)
                        v2=pig_vcva(wi,ui,vi2)*eval(yidx(ui,i),wi)
                        pig_tmp(i,ui,vi2)=v2
                        v1=v1+v2*py(vi2*pedn-pedn+i,1)
                      end do
                      dpy(pedn*ui-pedn+i,1)=v1
                    end do
                  end do

                  v2=0
                  do i=1,trn
                    v2=v2+pdpy(i,1)*dpy(i,1)
                  end do

                idum2=tn+ii3
                !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
                do ll=1,mn
                  do kk=1,ttn
                    if (ll==wi .and. kk==ii2) goto 54
                    idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                  end do
                  do kk=2,ttn
                    do mm=1,kk-1
                      !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                      idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                    end do
                  end do
                end do
54              continue

                  aim(idum+m_tnk(ii,k)+nbc,idum2+vi)=v2
                  !print*,'cov x other trait',idum+m_tnk(ii,k)+nbc,idum2+vi,tmpm(1,1)
                end do ! vi

                nbc2=0  !for cov for other trait, i.e. 1: (2,1),2: (3,1),3:(3,2),... 
                do vi2=2,m_tnk(ii2,wi)
                  do ui2=1,(vi2-1)
                    nbc2=nbc2+1

                    km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                    tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                    mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                    do jj=1,sub_tn(ii2)
                      do kk=1,sub_tn(ii2)
                        pig_vcva(wi,mm+jj,mm+kk)=tmp_vcva_2(jj,kk)
                      end do
                    end do

                    pig_tmp=0
                    do ui=1,tn
                      do i=1,pedn
                        v1=0
                        do vi=1,tn
                          !v2=pig_vcva(wi,ui,vi)*eval(i,wi)
                          v2=pig_vcva(wi,ui,vi)*eval(yidx(ui,i),wi)
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

                    aim(idum+m_tnk(ii,k)+nbc,idum2+m_tnk(ii2,wi)+nbc2)=v2 !check
                    !print*,'cov x other t cov',idum+m_tnk(ii,k)+nbc,idum2+m_tnk(ii2,wi)+nbc2
                  end do !ui2
                end do !vi2

                deallocate(phi_2,km2_2,tmp_wk_2,tmp_vcva_2)
              end do ! ii2

            !for cov x CVk   *************************************
            if (wi<k) then

              do ii2=2,ttn
                do jj2=1,ii2-1

                  allocate(phi_1(sub_tn(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(jj2,wi)))
                  allocate(tmp_wk_1(sub_tn(ii2)),tmp_vcva_2(sub_tn(ii2),sub_tn(jj2)))

                  tmp_wk_1=m_wk(ii2,1:sub_tn(ii2))
                  call legendre (m_tnk(ii2,wi),tmp_wk_1,sub_tn(ii2),phi_1)

                  allocate(phi_2(sub_tn(jj2),m_tnk(jj2,wi)),tmp_wk_2(sub_tn(jj2)))
                  tmp_wk_2=m_wk(jj2,1:sub_tn(jj2))
                  call legendre (m_tnk(jj2,wi),tmp_wk_2,sub_tn(jj2),phi_2)

                  nbc2=0
                  do vi=1,m_tnk(ii2,wi)
                    do vi2=1,m_tnk(jj2,wi)
                      nbc2=nbc2+1

                    km2_2=0;km2_2(vi,vi2)=1;pig_vcva=0
                    tmp_vcva_2=matmul(matmul(phi_1,km2_2),transpose(phi_2))
                    mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                    mm2=sum(sub_tn(1:jj2))-sub_tn(jj2)
                    do ll=1,sub_tn(ii2)
                      do kk=1,sub_tn(jj2)
                        pig_vcva(wi,mm+ll,mm2+kk)=tmp_vcva_2(ll,kk)
                        pig_vcva(wi,mm2+kk,mm+ll)=tmp_vcva_2(ll,kk)
                      end do
                    end do

                    pig_tmp=0
                    do ui=1,tn
                      do i=1,pedn
                        v1=0
                        do yi2=1,tn
                          !v2=pig_vcva(wi,ui,yi2)*eval(i,wi)
                          v2=pig_vcva(wi,ui,yi2)*eval(yidx(ui,i),wi)
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

                    idum2=tn+ii3
                    !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
                    do ll=1,mn
                      do kk=1,ttn
                        idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                      end do
                      do kk=2,ttn
                        do mm=1,kk-1
                          if (ll==wi .and. kk==ii2 .and. mm==jj2) goto 61
                          idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                        end do
                      end do
                    end do
61                  continue

                    aim(idum+m_tnk(ii,k)+nbc,idum2+nbc2)=v2 !check
                    !print*,'cov x CVk in prv',idum+m_tnk(ii,k)+nbc,idum2+nbc2

                    end do !vi2
                  end do ! vi

                  deallocate(phi_1,phi_2,km2_2)
                  deallocate(tmp_wk_1,tmp_wk_2,tmp_vcva_2)
                  !print*,'??'

                end do ! jj2
              end do ! ii2

            end if ! wk<k


            end do ! wi

            !dldv for cov *****************************
            tr1=0
            do i=1,pedn
              do ui=1,tn
                do vi=1,tn
                  tr1=tr1+pm(i,ui,vi)*pig(i,vi,ui)
                end do
              end do
            end do
            tr1=-0.5*tr1

            v2=0
            do i=1,trn
              v2=v2+pdpy(i,1)*tyv(i,1)
            end do

            dldv(idum+m_tnk(ii,k)+nbc,1)=v2*0.5+tr1
            !print*,idum+m_tnk(ii,k)+nbc,dldv(idum+m_tnk(ii,k)+nbc,1)
          end do !yi
        end do !xi

        !print*,'???'
        !deallocate(am,mm,phi,km2)
        deallocate(phi,km2,tmp_wk,tmp_vcva)
        !print*,'????'

      end do ! ii

      !CVk (diagonal) i.e. covariance between mrrm traits  **************************
      do ii=2,ttn
        do jj=1,ii-1
          allocate(phi(sub_tn(ii),m_tnk(ii,k)),km2(m_tnk(ii,k),m_tnk(jj,k)))
          allocate(tmp_wk(sub_tn(ii)),tmp_vcva(sub_tn(ii),sub_tn(jj)))

          tmp_wk=m_wk(ii,1:sub_tn(ii))
          call legendre (m_tnk(ii,k),tmp_wk,sub_tn(ii),phi)

          allocate(phi2(sub_tn(jj),m_tnk(jj,k)),tmp_wk2(sub_tn(jj)))
          tmp_wk2=m_wk(jj,1:sub_tn(jj))
          call legendre (m_tnk(jj,k),tmp_wk2,sub_tn(jj),phi2)

          nbc=0
          do xi=1,m_tnk(ii,k)
          do xi2=1,m_tnk(jj,k)
            nbc=nbc+1

            km2=0;km2(xi,xi2)=1;pig_vcva=0
            tmp_vcva=matmul(matmul(phi,km2),transpose(phi2))
            mm=sum(sub_tn(1:ii))-sub_tn(ii)
            mm2=sum(sub_tn(1:jj))-sub_tn(jj)
            do ll=1,sub_tn(ii)
              do kk=1,sub_tn(jj)
                pig_vcva(k,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
                pig_vcva(k,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
              end do
            end do

            pig=0
            do ui=1,tn
              do i=1,pedn
                v1=0
                do vi=1,tn
                  !v2=pig_vcva(k,ui,vi)*eval(i,k)
                  v2=pig_vcva(k,ui,vi)*eval(yidx(ui,i),k)
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

            idum=tn+ii3
            !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
            do ll=1,mn
              do kk=1,ttn
                idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
              end do
              do kk=2,ttn
                  do mm=1,kk-1
                  if (ll==k .and. kk==ii .and. mm==jj) goto 65
                    !idum=idum+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                    idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
                  end do
                end do
              end do
65            continue

              !aim(idum+xi,idum+xi)=v3
              aim(idum+nbc,idum+nbc)=v3
              !print*,'CVk',idum+nbc

              !for CVg x Ve (off diagonal)
              trnv=0
              do vi=1,tn
                trnv=trnv+rnm(vi)
                v2=0
                do i=1,rnm(vi)
                  v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnv-rnm(vi)+i,1)     !*dpy
                end do
                !aim(idum+xi,vi)=v2
                aim(idum+nbc,vi)=v2
                !print*,'Cvg x ve:',idum+xi,vi,v2
              end do

              !for CVk x cove ****************************************
              nbc2=0
            !do ii2=1,sub_tn(1)            !because the same # across traits
            do ii2=2,ttn
              do jj2=1,(ii2-1)
                do ii22=1,sub_tn(ii2)
                  do jj22=1,sub_tn(jj2)
                    if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
                      nbc2=nbc2+1
                      ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
                      vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

                    trnu=0
                    do i=1,ui
                      trnu=trnu+rnm(i)
                    end do
                    trnv=0
                    do i=1,vi
                      trnv=trnv+rnm(i)
                    end do

                    v2=0
                    do i=1,rnm(vi)
                      v2=v2+pdpy(trnu-rnm(ui)+i,1)*py(trnv-rnm(vi)+i,1) !*pdpy * pdy (same n for each trait)
                    end do
                    do i=1,rnm(ui)
                      v2=v2+pdpy(trnv-rnm(vi)+i,1)*py(trnu-rnm(ui)+i,1) !*pdpy * pdy (same n for each trait)
                    end do
                    aim(idum+nbc,tn+nbc2)=v2

                    end if

                  end do !jj22
                end do !ii22
              end do !ii2
            end do !jj2


              !for CVk(1~[i-1]) x CVk(i) *********************************
              nbc2=0
              do vi=1,m_tnk(ii,k)
              do vi2=1,m_tnk(jj,k)
                nbc2=nbc2+1
                if (nbc2<nbc) then

                km2=0;km2(vi,vi2)=1;pig_vcva=0
                tmp_vcva=matmul(matmul(phi,km2),transpose(phi2))
                mm=sum(sub_tn(1:ii))-sub_tn(ii)
                mm2=sum(sub_tn(1:jj))-sub_tn(jj)
                do ll=1,sub_tn(ii)
                  do kk=1,sub_tn(jj)
                    pig_vcva(k,mm+ll,mm2+kk)=tmp_vcva(ll,kk)
                    pig_vcva(k,mm2+kk,mm+ll)=tmp_vcva(ll,kk)
                  end do
                end do

                pig_tmp=0
                do ui=1,tn
                  do i=1,pedn
                    v1=0
                    do yi=1,tn
                      !v2=pig_vcva(k,ui,yi)*eval(i,k)
                      v2=pig_vcva(k,ui,yi)*eval(yidx(ui,i),k)
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
                !aim(idum+xi,idum+vi)=v2
                aim(idum+nbc,idum+nbc2)=v2
                !print*,'CVk x CVk:',idum+nbc,idum+nbc2,v2
                end if
              end do !vi2
              end do !vi

              !off diagonal between diff random eff. for CVk
              !******************************
              !do wi=1,k-1  !check in case multiple random effects
              do wi=1,k     !check in case multiple random effects
                do ii2=2,ttn
                  do jj2=1,ii2-1
                    if (wi==k .and. ii2>ii) then
                    elseif (wi==k .and. ii2==ii .and. jj2>=jj) then
                    else


                      allocate(phi_1(sub_tn(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(jj2,wi)))
                      allocate(tmp_wk_1(sub_tn(ii2)),tmp_vcva_2(sub_tn(ii2),sub_tn(jj2)))
  
                      tmp_wk_1=m_wk(ii2,1:sub_tn(ii2))
                      call legendre (m_tnk(ii2,wi),tmp_wk_1,sub_tn(ii2),phi_1)

                      allocate(phi_2(sub_tn(jj2),m_tnk(jj2,wi)),tmp_wk_2(sub_tn(jj2)))
                      tmp_wk_2=m_wk(jj2,1:sub_tn(jj2))
                      call legendre (m_tnk(jj2,wi),tmp_wk_2,sub_tn(jj2),phi_2)

                      nbc2=0
                      do vi=1,m_tnk(ii2,wi)
                      do vi2=1,m_tnk(jj2,wi)
                        nbc2=nbc2+1

                        km2_2=0;km2_2(vi,vi2)=1;pig_vcva=0
                        tmp_vcva_2=matmul(matmul(phi_1,km2_2),transpose(phi_2))
                        mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                        mm2=sum(sub_tn(1:jj2))-sub_tn(jj2)
                        do ll=1,sub_tn(ii2)
                          do kk=1,sub_tn(jj2)
                            pig_vcva(wi,mm+ll,mm2+kk)=tmp_vcva_2(ll,kk)
                            pig_vcva(wi,mm2+kk,mm+ll)=tmp_vcva_2(ll,kk)
                          end do
                        end do

                        pig_tmp=0
                        do ui=1,tn
                          do i=1,pedn
                            v1=0
                            do yi=1,tn
                              !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
                              v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
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

                        idum2=tn+ii3
                        !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
                        do ll=1,mn
                          do kk=1,ttn
                            idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                          end do
                          do kk=2,ttn
                            do mm=1,kk-1
                              if (ll==wi .and. kk==ii2 .and. mm==jj2) goto 70
                              !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                              idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                            end do
                          end do
                        end do
70                      continue

                        !aim(idum+xi,idum2+vi)=v2
                        aim(idum+nbc,idum2+nbc2)=v2
                        !print*,'CVk x CVk in prv',idum+nbc,idum2+nbc2,idum2,nbc2,v2
                      end do  ! vi2
                      end do ! vi

                      deallocate (phi_1,phi_2,km2_2,tmp_wk_1,tmp_vcva_2,tmp_wk_2)

                    end if
                  end do ! jj2
                end do ! ii2
              end do  ! wi


              ! CVk x Vk, off-diagonal between other traits 
              do wi=1,k
                do ii2=1,ttn
                  allocate(phi_2(sub_tn(ii2),m_tnk(ii2,wi)),km2_2(m_tnk(ii2,wi),m_tnk(ii2,wi)))
                  allocate(tmp_wk_2(sub_tn(ii2)),tmp_vcva_2(sub_tn(ii2),sub_tn(ii2)))

                  tmp_wk_2=m_wk(ii2,1:sub_tn(ii2))
                  call legendre (m_tnk(ii2,wi),tmp_wk_2,sub_tn(ii2),phi_2)

                  do vi=1,m_tnk(ii2,wi)   !check
 
                    km2_2=0;km2_2(vi,vi)=1;pig_vcva=0
                    tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                    mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                    do jj2=1,sub_tn(ii2)
                      do kk=1,sub_tn(ii2)
                        pig_vcva(wi,mm+jj2,mm+kk)=tmp_vcva_2(jj2,kk)
                      end do
                    end do

                    pig_tmp=0
                    do ui=1,tn
                      do i=1,pedn
                        v1=0
                        do yi=1,tn
                          !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
                          v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
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

                    idum2=tn+ii3
                    !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
                    do ll=1,mn
                      do kk=1,ttn
                        if (ll==wi .and. kk==ii2) goto 58
                        idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                      end do
                      do kk=2,ttn
                        do mm=1,kk-1
                          !idum2=idum2+ctnk(kk,mm,ll)+(ctnk(kk,mm,ll)**2-ctnk(kk,mm,ll))/2
                          idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                        end do
                      end do
                    end do
58                  continue

                    !aim(idum+xi,idum2+vi)=v2
                    aim(idum+nbc,idum2+vi)=v2
                    !print*,'CVk x Vk',idum+nbc,idum2+vi,idum2,vi,v2
                  end do ! vi

                  nbc2=0  ! CVk x cov, i.e. 1: (2,1),2:(3,1),3:(3,2),... 
                  do vi2=2,m_tnk(ii2,wi)
                    do ui2=1,(vi2-1)
                      nbc2=nbc2+1

                      km2_2=0;km2_2(ui2,vi2)=1;km2_2(vi2,ui2)=1;pig_vcva=0
                      tmp_vcva_2=matmul(matmul(phi_2,km2_2),transpose(phi_2))
                      mm=sum(sub_tn(1:ii2))-sub_tn(ii2)
                      do jj2=1,sub_tn(ii2)
                        do kk=1,sub_tn(ii2)
                          pig_vcva(wi,mm+jj2,mm+kk)=tmp_vcva_2(jj2,kk)
                        end do
                      end do

                    pig_tmp=0
                    do ui=1,tn
                      do i=1,pedn
                        v1=0
                        do yi=1,tn
                          !v2=pig_vcva(wi,ui,yi)*eval(i,wi)
                          v2=pig_vcva(wi,ui,yi)*eval(yidx(ui,i),wi)
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

                    !aim(idum+xi,idum2+m_tnk(ii2,wi)+nbc2)=v2 !check
                    aim(idum+nbc,idum2+m_tnk(ii2,wi)+nbc2)=v2 !check
                    !print*,'CVk x cov',idum+nbc,idum2+m_tnk(ii2,wi)+nbc2
                  end do !ui2
                end do !vi2

                deallocate(phi_2,km2_2,tmp_wk_2,tmp_vcva_2)
              end do !wi
            end do ! ii2

            !dldv CVk *********************************************************
            tr1=0
            do i=1,pedn
              do ui=1,tn
                do yi=1,tn
                  tr1=tr1+pm(i,ui,yi)*pig(i,yi,ui)
                end do
              end do
            end do
            tr1=-0.5*tr1

            v2=0
            do i=1,trn
              v2=v2+pdpy(i,1)*tyv(i,1)
            end do

            !dldv(idum+xi,1)=v2*0.5+tr1 !************************************
            dldv(idum+nbc,1)=v2*0.5+tr1 !************************************
            !print*,idum+nbc,dldv(idum+nbc,1)
            !print*,'dldv:',tn
          end do !xi2
          end do !xi

          deallocate(phi,km2,tmp_wk,tmp_vcva,phi2,tmp_wk2)
        end do  ! jj
      end do  ! ii

    end do !k

    call cpu_time (t2_cpu)
    print*,'derivatives done - time:',real(t2_cpu-t1_cpu)
    print*,''
 
    idum=tn+ii3
    !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
    do ll=1,mn
      do kk=1,ttn
        idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
      end do
      do kk=2,ttn
        do mm=1,kk-1
          idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
        end do
      end do
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

    call cholesky_inv (aim,idum,x1)
    print *,''

    if (LKH-MLKH.gt.0) then
      MLKH=LKH
      sdmI=aim
    end if
    print*,''

    if (ABS(LKH-LKHP).lt.conv) goto 1000

    up=MATMUL(aim,dldv)
    !do i=1,idum
    !  print*,real(up(i,1)),real(dldv(i,1))
    !end do
    !pause

111 continue
    if (itit==1) then
      vcve_rsv=vcve
      m_tkm_rsv=m_tkm
      ctkm_rsv=ctkm
      if (zi.ge.nit) then
        print*,'Likelihood not converged >>> may need a longer iterations'
        print*,''
        goto 1000
      end if

    else
      up=up*0.7**(itit-1)
      vcve=vcve_rsv
      m_tkm=m_tkm_rsv
      ctkm=ctkm_rsv
    end if

    do xi=1,tn
      vcve(xi,xi)=vcve(xi,xi)+up(xi,1)
      !if (vcve(xi,xi)<0) vcve(xi,xi)=0.000001
    end do
    nbc2=0
  !do ii2=1,sub_tn(1)            !because the same # across traits
  do ii2=2,ttn
    do jj2=1,(ii2-1)
      do ii22=1,sub_tn(ii2)
        do jj22=1,sub_tn(jj2)
          if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
            nbc2=nbc2+1
            ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
            vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

          vcve(ui,vi)=vcve(ui,vi)+up(tn+nbc2,1)

          end if
        end do
      end do
    end do
  end do

    do wi=1,mn
      do ii=1,ttn

        idum=tn+ii3
        !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
        do ll=1,mn
          do kk=1,ttn
            if (ll==wi .and. kk==ii) goto 80
            idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
          end do
          do kk=2,ttn
            do mm=1,kk-1
              idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
            end do
          end do
        end do
80      continue

        do xi=1,m_tnk(ii,wi)
          m_tkm(wi,ii,xi,xi)=m_tkm(wi,ii,xi,xi)+up(idum+xi,1)
          !print*,idum+xi
        end do
        nbc=0
        do xi=2,m_tnk(ii,wi)
          do yi=1,xi-1
            nbc=nbc+1
            m_tkm(wi,ii,xi,yi)=m_tkm(wi,ii,xi,yi)+up(idum+m_tnk(ii,wi)+nbc,1)
            m_tkm(wi,ii,yi,xi)=m_tkm(wi,ii,xi,yi)
            !print*,idum+m_tnk(ii,wi)+nbc
          end do
        end do
      end do ! ii

      do ii=1,ttn
        do jj=1,ii-1
          nbc2=0
          do xi=1,m_tnk(ii,wi)
            do yi=1,m_tnk(jj,wi)
              nbc2=nbc2+1

              idum2=tn+ii3
              !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
              do ll=1,mn
                do kk=1,ttn
                   idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
                end do
                do kk=2,ttn
                  do mm=1,kk-1
                    if (ll==wi .and. kk==ii .and. mm==jj) goto 81
                    idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                  end do
                end do
              end do
81            continue

              ctkm(wi,ii,jj,xi,yi)=ctkm(wi,ii,jj,xi,yi)+up(idum2+nbc2,1)
              !print*,'*',idum+m_tnk(ii,wi)+nbc+nbc2,idum2+nbc2
            end do
          end do

        end do ! jj
      end do ! ii
    end do !wi

    !print*,vcve,vcva

    !do i=1,tnk
    !  print*,km(i,:)
    !end do
    !pause

  END do !zi


  !PRINT*,'LKH not converged na na na'
  !WRITE(41,*)'LKH not converged na na na'
  LKH=MLKH


1000 continue 

  open (UNIT=41,FILE=fl5,STATUS='unknown')

  do xi=1,tn
    write(41,'(a7,100f12.4)')'Ve',vcve(xi,xi),sqrt(sdmI(xi,xi))
  end do
  nbc2=0
  !do ii2=1,sub_tn(1)            !because the same # across traits
  do ii2=2,ttn
    do jj2=1,(ii2-1)
      do ii22=1,sub_tn(ii2)
        do jj22=1,sub_tn(jj2)
          if (m_wk(ii2,ii22)==m_wk(jj2,jj22)) then
            nbc2=nbc2+1
            ui=sum(sub_tn(1:ii2))-sub_tn(ii2)+ii22
            vi=sum(sub_tn(1:jj2))-sub_tn(jj2)+jj22

        write(41,'(a7,100f12.4)')'cove',vcve(ui,vi),sqrt(sdmI(tn+nbc2,tn+nbc2))

          end if
        end do
      end do
    end do
  end do
  write(41,*)''

  do xi=1,tn
    sum_v(xi)=vcve(xi,xi)
  end do
  do wi=1,mn
    do ii=1,ttn

          idum=tn+ii3
          !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
          do ll=1,mn
            do kk=1,ttn
              if (ll==wi .and. kk==ii) goto 95
              idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
            end do
            do kk=2,ttn
              do mm=1,kk-1
                idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
              end do
            end do
          end do
95        continue

      do xi=1,m_tnk(ii,wi)
        write(41,'(a7,100f12.4)')'Vk',m_tkm(wi,ii,xi,xi),sqrt(sdmI(idum+xi,idum+xi))
      end do

      nbc=0
      do xi=2,m_tnk(ii,wi)
        do yi=1,xi-1
          nbc=nbc+1
          write(41,'(a7,100f12.4)')'cov',m_tkm(wi,ii,xi,yi),sqrt(sdmI(idum+m_tnk(ii,wi)+nbc,idum+m_tnk(ii,wi)+nbc))
        end do
      end do

      write(41,*)''
    end do ! ii

    do ii=2,ttn
      do jj=1,ii-1
        nbc2=0
        do xi=1,m_tnk(ii,wi)
          do yi=1,m_tnk(jj,wi)

            idum2=tn+ii3
            !idum2=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
            do ll=1,mn
              do kk=1,ttn
                 idum2=idum2+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
              end do
              do kk=2,ttn
                do mm=1,kk-1
                  if (ll==wi .and. kk==ii .and. mm==jj) goto 96
                  idum2=idum2+m_tnk(kk,ll)*m_tnk(mm,ll)
                end do
              end do
            end do
96          continue

            nbc2=nbc2+1
            write(41,'(a7,100f12.4)')'CVk',ctkm(wi,ii,jj,xi,yi),sqrt(sdmI(idum2+nbc2,idum2+nbc2))
          end do
        end do
      end do !jj
      write(41,*)''
    end do !ii

  end do ! wi

  write(41,'(a7,100f12.4)')'LKH',LKH  !,real(zi) !indicate convergence
  write(41,*)''

  write(41,*)'transformed ***'
  do wi=1,mn
    do xi=1,tn
      write(41,'(a7,100f12.4)')'Va',vcva(wi,xi,xi)!,sqrt(sdmI(tn+wi*nb-nb+xi,tn+wi*nb-nb+xi))
      sum_v(xi)=sum_v(xi)+vcva(wi,xi,xi)
    end do
    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        write(41,'(a7,100f12.4)')'cov',vcva(wi,xi,yi)!,sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
      end do
    end do
    write(41,*)''
  end do

  do wi=1,mn
    do xi=1,tn
      x1=vcva(wi,xi,xi)/sum_v(xi)
      write(41,'(a7,100f12.4)')'h2',x1!,sqrt(x2)
    end do
    write(41,*)''

    nbc=0
    do xi=2,tn
      do yi=1,xi-1
        nbc=nbc+1
        z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
        write(41,'(a7,100f12.4)')'cor',z1!,sqrt(z2)
      end do
    end do
    write(41,*)''
  end do

  idum=tn+ii3
  !idum=tn+sub_tn(1)*(ttn**2/2-ttn/2)   !for residual correlation
  do ll=1,mn
    do kk=1,ttn
      idum=idum+m_tnk(kk,ll)+(m_tnk(kk,ll)**2-m_tnk(kk,ll))/2
    end do
    do kk=2,ttn
      do mm=1,kk-1
        idum=idum+m_tnk(kk,ll)*m_tnk(mm,ll)
      end do
    end do
  end do

  write(41,*)''
  do i=1,idum
    write(41,'(100g18.8)')(sdmI(i,1:i))
  end do

  close(41)

if (fl11.ne.'null') then

  print*,'for BLUP solution after convergence, check again *******************'

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
    allocate(gr_km_phi(sum(m_tnk(1:ttn,k)),tn))

    do ii=1,ttn
      allocate(phi(sub_tn(ii),m_tnk(ii,k)),km(m_tnk(ii,k),m_tnk(ii,k)))
      allocate(tmp_wk(sub_tn(ii)),km_phi(m_tnk(ii,k),sub_tn(ii)))

      tmp_wk=m_wk(ii,1:sub_tn(ii))
      call legendre (m_tnk(ii,k),tmp_wk,sub_tn(ii),phi)
      km=m_tkm(k,ii,1:m_tnk(ii,k),1:m_tnk(ii,k))

      km_phi=matmul(km,transpose(phi))

      do i=1,m_tnk(ii,k)
        do j=1,sub_tn(ii)
          i2=sum(m_tnk(1:ii,k))-m_tnk(ii,k)+i
          j2=sum(sub_tn(1:ii))-sub_tn(ii)+j
          gr_km_phi(i2,j2)=km_phi(i,j)
        end do
      end do
      deallocate(phi,km,tmp_wk,km_phi)

    end do ! ii

    do ii=2,ttn   !two traits only, check when > 2 trait 
      do jj=1,ii-1
        !allocate(phi(sub_tn(ii),m_tnk(ii,k)),km(m_tnk(ii,k),m_tnk(jj,k)))
        !allocate(tmp_wk(sub_tn(ii)),km_phi(m_tnk(ii,k),sub_tn(jj)))
        allocate(phi(sub_tn(jj),m_tnk(jj,k)),km(m_tnk(ii,k),m_tnk(jj,k)))
        allocate(tmp_wk(sub_tn(jj)),km_phi(m_tnk(ii,k),sub_tn(jj)))

        !tmp_wk=m_wk(ii,1:sub_tn(ii))
        tmp_wk=m_wk(jj,1:sub_tn(jj))
        !call legendre (m_tnk(ii,k),tmp_wk,sub_tn(ii),phi)
        call legendre (m_tnk(jj,k),tmp_wk,sub_tn(jj),phi)

        !allocate(phi2(sub_tn(jj),m_tnk(jj,k)),tmp_wk2(sub_tn(jj)))
        !allocate(km_phi2(m_tnk(jj,k),sub_tn(ii)))
        allocate(phi2(sub_tn(ii),m_tnk(ii,k)))
        allocate(tmp_wk2(sub_tn(ii)),km_phi2(m_tnk(jj,k),sub_tn(ii)))
        !tmp_wk2=m_wk(jj,1:sub_tn(jj))
        tmp_wk2=m_wk(ii,1:sub_tn(ii))
        !call legendre (m_tnk(jj,k),tmp_wk2,sub_tn(jj),phi2)
        call legendre (m_tnk(ii,k),tmp_wk2,sub_tn(ii),phi2)

        km=ctkm(k,ii,jj,1:m_tnk(ii,k),1:m_tnk(jj,k))

        km_phi=matmul(km,transpose(phi))

        do i=1,m_tnk(ii,k)
          do j=1,sub_tn(jj)
            i2=sum(m_tnk(1:ii,k))-m_tnk(ii,k)+i
            j2=sum(sub_tn(1:jj))-sub_tn(jj)+j
            gr_km_phi(i2,j2)=km_phi(i,j)
          end do
        end do

        km_phi2=matmul(transpose(km),transpose(phi2))

        do i=1,m_tnk(jj,k)
          do j=1,sub_tn(ii)
            i2=sum(m_tnk(1:jj,k))-m_tnk(jj,k)+i
            j2=sum(sub_tn(1:ii))-sub_tn(ii)+j
            gr_km_phi(i2,j2)=km_phi2(i,j)
          end do
        end do

        deallocate(phi,km,tmp_wk,km_phi,phi2,tmp_wk2,km_phi2)

      end do ! ii
    end do ! jj

      do xi=1,sum(m_tnk(1:ttn,k))    !check
        do zi=1,pedn
          v1=0;v2=0
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            !trny=trny+pedn     !becuase equal and complete design
            !v1=v1+eval(zi,k)*gr_km_phi(xi,yi)*py(trny-pedn+zi,1)
            !v1=v1+eval(zi,k)*gr_km_phi(xi,yi)*py(trny-rnm(yi)+zi,1)
            v1=v1+eval(zi,k)*gr_km_phi(xi,yi)*py(trny-rnm(yi)+yidx2(zi),1)
          end do
          blup_ebv(zi,xi,k)=v1
          !print*,zi,xi,k,blup_ebv(zi,xi,k)
        end do !zi
      end do !xi
      deallocate(gr_km_phi)
  end do !k

  !transformed back
    do xi=1,tn
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

    !do k=1,mn
    !  do xi=1,sum(m_tnk(1:ttn,k))    !check
    !    do zi=1,pedn
    !      v1=0
    !      do j=1,pedn
    !        v1=v1+mbin(1,zi,j)*blup_ebv(j,xi,k)    !rrm e use the same evec
    !      end do
    !      blup_ebv2(zi,xi,k)=v1
    !    end do !k
    !  end do !zi
    !end do !xi

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
  !write(43,'(a6,a14,a28)')'site','ordered ind#','EBVs ...'
  write(43,'(a6,a14,a9,a14,a24)')'trait','random eff.#','RR order','ordered ind#','EBVs'

  do ii=1,ttn
  do k=1,mn
  !do xi=1,tn    !check
  do xi=1,m_tnk(ii,k)    !check
    !print*,xi,m_tnk(ii,k),sum(m_tnk(1:ii,k))
    do i=1,pedn
      !write(43,'(i6,i14,100f24.16)')xi,i,blup_ebv(i,xi,:)
      write(43,'(i6,i14,i9,i14,f24.16)')ii,k,xi,i,blup_ebv2(i,sum(m_tnk(1:ii,k))-m_tnk(ii,k)+xi,k)
    end do
  end do
  end do
  end do

  close(43)

end if

end subroutine

