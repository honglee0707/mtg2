
subroutine aireml_b_h_rg (mbin,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,rg,nit,conv)        

!***********************************************************************
!aireml_b: bivariate analysis especially for case control 
!S. Hong Lee

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects should be = 1 (mean)
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,nb,nbc,nbc2,xi,yi,ri,zi,vi,ui,wi,wj,nit
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu
integer::yidx(tn,rn),rnm(tn),fnm(tn)
integer::rn1,rn2,yidx1(rn),yidx2(rn)
!double precision::mbin(mn,pedn,pedn),yv(rn,2),tyv(rn*2,1)
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1)

!double precision::qq(mn,rn,rn),qq2(mn,rn,rn)
!real::qq(mn,rn,rn),qq2(mn,rn,rn)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),rg,conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::blup_ebv(pedn,tn,mn),beta(tfn,1) !,cor
integer::i,j,m,k,l,io
double precision::sum_h2,sum_v(tn)


double precision::LC,LA,LG,ypy,x1,x2,x3,x10,x11,v1,v2,v10,v11
double precision::y1,y2,y3,y10,y11,z1,z2,z3,z10,z11
double precision::LKHP,MLKH,mva(mn),mva2(mn),mve,mve2,mcov(mn)
double precision::h2(mn),h2n(mn),tr1,tr2(mn),tr3

double precision::xm(trn,tfn),xmm(tn,pedn,1000)

double precision::xvmx(tfn,tfn),xvm(tfn,trn) 
double precision::xvm2(tfn,trn),xvm3(trn,tfn) 
double precision::pm(trn,trn),py(trn,1),xb(trn,1)

double precision::v_tmp(trn),v2_tmp(trn)

!double precision::aim(2+mn*3,2+mn*3),aiminv(2+mn*3,2+mn*3),up(2+mn*3,1),dldv(2+mn*3,1)
!double precision::sdm(2+mn*3,2+mn*3),sdmI(2+mn*3,2+mn*3)
double precision::aim(tn+mn*tn,tn+mn*tn)
double precision::sdmI(tn+mn*tn,tn+mn*tn)
double precision::up(tn+mn*tn,1),dldv(tn+mn*tn,1)

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
character(len=64)::fl4,fl5,fl11
logical::file_exist


xm=0;trn=0;tfn=0
do xi=1,tn
  tfn=tfn+fnm(xi)
  do i=1,rn
    !print*,xm1(i,:)
    if (yv(i,xi).ne.-99999) then
      trn=trn+1
      tyv(trn,1)=yv(i,xi)

      do j=1,fnm(xi)
        xm(trn,tfn-fnm(xi)+j)=xmm(xi,i,j)
      end do
      !print*,xm(rn1,1:fn1),xm1(i,:),rn1,i
      !pause
    end if
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

  ! Start iteration
  LKH=-1000000000
  MLKH=-1000000000


  ! V matrix    V= ZAZ'Va+Z2GZ2'Vpe+IVe*************************
  !do zi=1,20
  do zi=1,nit
    LKHP=LKH
  !LKH and AIREML (without sparce technique)
  ! V matrix    V = ZAZ'Va + IVe *************************

  call cpu_time (t2_cpu)

  !specify covariance making genetic correlation=rg 
  do xi=2,tn
    do yi=1,xi-1
      do k=1,mn
        vcva(k,xi,yi)=rg*sqrt(vcva(k,xi,xi))*sqrt(vcva(k,yi,yi))
      end do
    end do
  end do

  pm=0;trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    trny=0
    do yi=1,xi
      trny=trny+rnm(yi)

      do i=1,rnm(xi)
        do j=1,rnm(yi)
          do k=1,mn
            pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)+mbin(k,yidx(xi,i),yidx(yi,j))*vcva(k,xi,yi)
          end do
          if (xi.ne.yi) then
            pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end if
        end do
        if (xi==yi) then
          pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)
        end if

      end do
    end do

  end do
  !print*,trn,trnx,trny
  !print*,pm(1, 1:3)
  !print*,pm(trn, trn-3:trn)
  !pause





  call cpu_time (t2_cpu)
  call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA

    ! P = (VI)-(VI X (X'VI X)I' X' VI)******************************
    !xvm=MATMUL(transpose(xm),pm)
    do i=1,tfn
      do j=1,trn
        v1=0
        do k=1,trn
          v1=v1+xm(k,i)*pm(k,j)    !xm*pm
        end do
        xvm(i,j)=v1
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

    do i=1,tfn
      v2=0
      do j=1,trn
        v1=0
        do k=1,trn
          v1=v1+xvm2(i,k)*pm(k,j)    !*pm
        end do
        v2=v2+v1*tyv(j,1)           !*tyv
        xvm3(j,i)=v1
      end do
      beta(i,1)=v2
    end do
   !*********************************************
   !print*,'beta:',beta

    !vmx4=MATMUL(MATMUL(MATMUL(vmI,xm),xvmxI),MATMUL(transpose(xm),vmI))
    !pm=MATMUL(MATMUL(MATMUL(vmI,xm),xvmxI),MATMUL(transpose(xm),vmI))
    !do i=1,trn
    !  do j=1,trn
    !    !pm(i,j)=vmI(i,j)-vmx4(i,j)
    !    !pm(i,j)=vmI(i,j)-pm(i,j)
    !    pm(i,j)=vm(i,j)-pm(i,j)
    !  end do
    !end do
    !pm=pm-MATMUL(MATMUL(MATMUL(pm,xm),xvmx),MATMUL(transpose(xm),pm))
    do i=1,trn
      do j=1,trn
        v1=0
        do k=1,tfn
          v1=v1+xvm3(i,k)*xvm(k,j)
        end do
        pm(i,j)=pm(i,j)-v1
      end do
    end do

    !ypy estimation
    !ypyf=MATMUL(MATMUL(transpose(tyv),pm),tyv)
    v2=0
    do i=1,trn
      v1=0
      do j=1,trn
        v1=v1+tyv(j,1)*pm(j,i)
      end do
      v2=v2+v1*tyv(i,1)
    end do

    !LKH=-0.5*(LA+LG+ypyf(1,1))
    LKH=-0.5*(LA+LG+v2)
    !print*,LA,LG,v2  !ypyf

    call cpu_time (t2_cpu)
    !print*,'LKH:',t2_cpu-t1_cpu

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
    py=MATMUL(pm,tyv)       !for py

    !print*,'ai matrix'
    call cpu_time (t1_cpu)
    aim=0
    dldv=0

    !nb=tn+(tn**2-tn)/2  !no block for random effects (Vg, cov)
    nb=tn  !no block for random effects (Vg, cov)
    !Ve (diagonal)  **********************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      v2=0
      do i=1,rnm(xi)
        v1=0
        do j=1,rnm(xi)
          !v1=v1+py(j,1)*pm(j,i)
          v1=v1+py(trnx-rnm(xi)+j,1)*pm(trnx-rnm(xi)+j,trnx-rnm(xi)+i)
        end do
        v2=v2+v1*py(trnx-rnm(xi)+i,1)
      end do 
      aim(xi,xi)=v2
    end do

    !dldv Ve ***********************************************************
    trnx=0
    do xi=1,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      tr1=0
      do i=1,rnm(xi)
        tr1=tr1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)
      end do
      tr1=-0.5*tr1
      v2=0
      do i=1,rnm(xi)
        v2=v2+py(trnx-rnm(xi)+i,1)*py(trnx-rnm(xi)+i,1)
      end do 
      dldv(xi,1)=v2*(0.5)+tr1
    end do


    !Ve (off diagonal) aim(2,1) *************************************
    trnx=rnm(1)
    do xi=2,tn
      trnx=trnx+rnm(xi)   !total rn upto xi
      trny=0
      do yi=1,(xi-1)
        trny=trny+rnm(yi)  !total rn upto yi

        v2=0
        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(yi)
            v1=v1+py(trny-rnm(yi)+j,1)*pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)
          end do
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 
        aim(xi,yi)=v2
        !print*,v2
      end do
    end do
    !pause

    !print*,'Vg'
    !for Vg
    do k=1,mn !no. random effects

      !for Vg (doagonal) ******************************************************* 
      trnx=0
      do xi=1,tn
        trnx=trnx+rnm(xi)
        v2=0
        do i=1,rnm(xi)
          v10=0
          do j=1,rnm(xi)
            v10=v10+py(trnx-rnm(xi)+j,1)*mbin(k,yidx(xi,j),yidx(xi,i))     !py*pig
          end do
          !when rg=1, covariance part is not zero **************************
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            if (xi.ne.yi) then
              do j=1,rnm(yi)
                v10=v10+py(trny-rnm(yi)+j,1)*mbin(k,yidx(yi,j),yidx(xi,i))*0.5*(sqrt(vcva(k,yi,yi))/sqrt(vcva(k,xi,xi)))
              end do
            end if
          end do
          !*****************************************************************
          v_tmp(trnx-rnm(xi)+i)=v10
        end do

        !when rg=1, covariance part is not zero **************************
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          if (xi.ne.yi) then
            do j=1,rnm(yi)
              v10=0
              do i=1,rnm(xi)
                v10=v10+py(trnx-rnm(xi)+i,1)*mbin(k,yidx(xi,i),yidx(yi,j))*0.5*(sqrt(vcva(k,yi,yi)/vcva(k,xi,xi)))
              end do
              v_tmp(trny-rnm(yi)+j)=v10
            end do
          end if
        end do
        !******************************************************************

        do i=1,trn
          v10=0
          !do j=1,rnm(xi)
          do j=1,trn                  !*** rg=1 ***
            !v10=v10+v_tmp(trnx-rnm(xi)+j)*pm(trnx-rnm(xi)+j,i)      !*pm
            v10=v10+v_tmp(j)*pm(j,i)       !*pm   *** for rg=1 ***
          end do
          v2_tmp(i)=v10
        end do



        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(xi)
            v1=v1+v2_tmp(trnx-rnm(xi)+j)*mbin(k,yidx(xi,j),yidx(xi,i))     !*pig
          end do
          !when rg=1, covariance part is not zero **************************
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            if (xi.ne.yi) then
              do j=1,rnm(yi)
                v1=v1+v2_tmp(trny-rnm(yi)+j)*mbin(k,yidx(yi,j),yidx(xi,i))*0.5*(sqrt(vcva(k,yi,yi)/vcva(k,xi,xi)))
              end do
            end if
          end do
          !*****************************************************************
          v2=v2+v1*py(trnx-rnm(xi)+i,1)                                !*py
        end do

        !when rg=1, covariance part is not zero **************************
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          if (xi.ne.yi) then
            do j=1,rnm(yi)
              v1=0
              do i=1,rnm(xi)
                v1=v1+v2_tmp(trnx-rnm(xi)+i)*mbin(k,yidx(xi,i),yidx(yi,j))*0.5*(sqrt(vcva(k,yi,yi)/vcva(k,xi,xi)))
              end do
              v2=v2+v1*py(trny-rnm(yi)+j,1)
            end do
          end if
        end do
        !******************************************************************

        aim(tn+nb*k-nb+xi,tn+nb*k-nb+xi)=v2

        !for Vg x Ve (off diagonal) !******************************************
        trnv=0
        do vi=1,tn
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            v2=v2+v2_tmp(trnv-rnm(vi)+i)*py(trnv-rnm(vi)+i,1)                                !*py
          end do
          !aim(2+k*3-3+xi,vi)=v2             !with Ve1
          aim(tn+nb*k-nb+xi,vi)=v2
        end do

        !for Vg(1~[i-1]) x Vg(i)
        trnv=0
        do vi=1,xi-1
          trnv=trnv+rnm(vi)
          v2=0
          do i=1,rnm(vi)
            v1=0
            do j=1,rnm(vi)
              v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(k,yidx(vi,j),yidx(vi,i))   !*pig_tmp
            end do

            !when rg=1, covariance part is not zero **************************
            trnu=0
            do ui=1,tn
              trnu=trnu+rnm(ui)
              if (ui.ne.vi) then
                do j=1,rnm(ui)
                  v1=v1+v2_tmp(trnu-rnm(ui)+j)*mbin(k,yidx(ui,j),yidx(vi,i))*0.5*(sqrt(vcva(k,ui,ui)/vcva(k,vi,vi)))
                end do
              end if
            end do
            !*****************************************************************
            v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
          end do

          !when rg=1, covariance part is not zero **************************
          trnu=0
          do ui=1,tn
            trnu=trnu+rnm(ui)
            if (ui.ne.vi) then
              do j=1,rnm(ui)
                v1=0
                do i=1,rnm(vi)
                  v1=v1+v2_tmp(trnv-rnm(vi)+i)*mbin(k,yidx(vi,i),yidx(ui,j))*0.5*(sqrt(vcva(k,ui,ui)/vcva(k,vi,vi)))
                end do
                v2=v2+v1*py(trnu-rnm(ui)+j,1)
              end do
            end if
          end do
          !******************************************************************

          aim(tn+nb*k-nb+xi,tn+nb*k-nb+vi)=v2
        end do



        !check 
        !off diagonal *********************************************************
        do wi=1,k-1
          trnv=0
          do vi=1,tn
            trnv=trnv+rnm(vi)
            v2=0
            do i=1,rnm(vi)
              v1=0
              do j=1,rnm(vi)
                v1=v1+v2_tmp(trnv-rnm(vi)+j)*mbin(wi,yidx(vi,j),yidx(vi,i))     !*pig_tmp
              end do
              v2=v2+v1*py(trnv-rnm(vi)+i,1)                                !*py
            end do
            !aim(2+k*3-3+xi,2+wi*3-3+vi)=v2   !with the xi th Vg for first trait
            aim(tn+nb*k-nb+xi,tn+nb*wi-nb+vi)=v2   !with the xi th Vg for first trait
            !print*,v2
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
          v1=0
          do j=1,rnm(xi)
            v1=v1+pm(trnx-rnm(xi)+i,trnx-rnm(xi)+j)*mbin(k,yidx(xi,j),yidx(xi,i))
          end do
          !when rg=1, covariance part is not zero **************************
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            if (xi.ne.yi) then
              do j=1,rnm(yi)
                v1=v1+pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)*mbin(k,yidx(yi,j),yidx(xi,i))*0.5*(sqrt(vcva(k,yi,yi)/vcva(k,xi,xi)))
              end do
            end if
          end do
          !*****************************************************************
          tr1=tr1+v1
        end do

        !when rg=1, covariance part is not zero **************************
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          if (xi.ne.yi) then
            do j=1,rnm(yi)
              v1=0
              do i=1,rnm(xi)
                v1=v1+pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)*mbin(k,yidx(xi,i),yidx(yi,j))*0.5*(sqrt(vcva(k,yi,yi)/vcva(k,xi,xi)))
              end do
              tr1=tr1+v1
            end do
          end if
        end do
        !******************************************************************
        tr1=-0.5*tr1

        v2=0
        do i=1,rnm(xi)
          v1=0
          do j=1,rnm(xi)
            v1=v1+py(trnx-rnm(xi)+j,1)*mbin(k,yidx(xi,j),yidx(xi,i))
          end do

          !when rg=1, covariance part is not zero **************************
          trny=0
          do yi=1,tn
            trny=trny+rnm(yi)
            if (xi.ne.yi) then
              do j=1,rnm(yi)
                v1=v1+py(trny-rnm(yi)+j,1)*mbin(k,yidx(yi,j),yidx(xi,i))*0.5*(sqrt(vcva(k,yi,yi))/sqrt(vcva(k,xi,xi)))
              end do
            end if
          end do
          !*****************************************************************
          v2=v2+v1*py(trnx-rnm(xi)+i,1)
        end do 

        !when rg=1, covariance part is not zero **************************
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          if (xi.ne.yi) then
            do j=1,rnm(yi)
              v1=0
              do i=1,rnm(xi)
                v1=v1+py(trnx-rnm(xi)+i,1)*mbin(k,yidx(xi,i),yidx(yi,j))*0.5*(sqrt(vcva(k,yi,yi)/vcva(k,xi,xi)))
              end do
              v2=v2+v1*py(trny-rnm(yi)+j,1)
            end do
          end if
        end do
        !******************************************************************
        dldv(tn+nb*k-nb+xi,1)=v2*0.5+tr1 !************************************
      end do

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
    !print*,''
    !do i=1,tn+nb*mn
    !  print*,real(dldv(i,1))
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

    do xi=1,tn
      vcve(xi,xi)=vcve(xi,xi)+up(xi,1)
    end do
    do wi=1,mn
      do xi=1,tn
        vcva(wi,xi,xi)=vcva(wi,xi,xi)+up(tn+nb*wi-nb+xi,1)
        if (vcva(wi,xi,xi).le.0) vcva(wi,xi,xi)=0.001
      end do
      !nbc=0
      !do xi=2,tn
      !  do yi=1,xi-1
      !    nbc=nbc+1
      !    vcva(wi,xi,yi)=vcva(wi,xi,yi)+up(tn+nb*wi-nb+tn+nbc,1)
      !  end do
      !end do
    end do

    !print*,ve,ve2,va,va2,cov
    !print*,vcve,vcva

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
    !nbc=0
    !do xi=2,tn
    !  do yi=1,xi-1
    !    nbc=nbc+1
    !    write(41,'(a7,100f12.4)')'cov',vcva(wi,xi,yi),sqrt(sdmI(tn+wi*nb-nb+tn+nbc,tn+wi*nb-nb+tn+nbc))
    !  end do
    !end do
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
      x1=vcva(wi,xi,xi)/sum_v(xi)
      x11=sdmI(tn+nb*wi-nb+xi,xi) !for the row ************

      do wj=1,mn
        x11=x11+sdmI(tn+nb*wi-nb+xi,tn+nb*wj-nb+xi)
      end do
      !print*,x11

      !SE **********************************
      x2=x1**2*(sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/vcva(wi,xi,xi)**2+vvc(xi)/sum_v(xi)**2-2*x11/(vcva(wi,xi,xi)*sum_v(xi)))
      write(41,'(a7,100f12.4)')'h2',x1,sqrt(x2)

    end do
    write(41,*)''

    !nbc=0
    !do xi=2,tn
    !  do yi=1,xi-1
    !    nbc=nbc+1
    !    !z1=cov(i)/sqrt(va(i)*va2(i))
    !    z1=vcva(wi,xi,yi)/sqrt(vcva(wi,xi,xi)*vcva(wi,yi,yi))
    !    !var(cor) *************************************
    !    !z2=sdmI(2+i*3-2,2+i*3-2)/(4*va(i)**2)
    !    z2=sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+xi)/(4*vcva(wi,xi,xi)**2)
    !    !z2=z2+sdmI(2+i*3-1,2+i*3-1)/(4*va2(i)**2)
    !    z2=z2+sdmI(tn+nb*wi-nb+yi,tn+nb*wi-nb+yi)/(4*vcva(wi,yi,yi)**2)
    !    !z2=z2+sdmI(2+i*3,2+i*3)/(cov(i)**2)
    !    z2=z2+sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+tn+nbc)/(vcva(wi,xi,yi)**2)
    !    !z2=z2+2*sdmI(2+i*3-2,2+i*3-1)/(4*va(i)*va2(i))
    !    z2=z2+2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+yi)/(4*vcva(wi,xi,xi)*vcva(wi,yi,yi))
    !    !z2=z2-2*sdmI(2+i*3-2,2+i*3)/(2*va(i)*cov(i))
    !    z2=z2-2*sdmI(tn+nb*wi-nb+xi,tn+nb*wi-nb+tn+nbc)/(2*vcva(wi,xi,xi)*vcva(wi,xi,yi))
    !    !z2=z2-2*sdmI(2+i*3,2+i*3-1)/(2*cov(i)*va2(i))
    !    z2=z2-2*sdmI(tn+nb*wi-nb+tn+nbc,tn+nb*wi-nb+yi)/(2*vcva(wi,xi,yi)*vcva(wi,yi,yi))
    !    z2=z2*(z1**2)
    !    write(41,'(a7,100f12.4)')'cor',z1,sqrt(z2)
    !  end do
    !end do
    write(41,*)''
  end do

  write(41,*)''
  write(41,'(a7,100f12.4)')'LKH',LKH

  write(41,*)''
  do i=1,(tn+nb*mn)
    !write(41,*)real(sdmI(i,1:i))
    write(41,'(100g18.8)')(sdmI(i,1:i))
  end do

  close(41)

if (fl11.ne.'null') then

  print*,'for BLUP solution after convergence *******************'
  pm=0;trnx=0
  do xi=1,tn
    trnx=trnx+rnm(xi)
    trny=0
    do yi=1,xi
      trny=trny+rnm(yi)

      do i=1,rnm(xi)
        do j=1,rnm(yi)
          do k=1,mn
            pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)+mbin(k,yidx(xi,i),yidx(yi,j))*vcva(k,xi,yi)
            !print*,pm(trnx-rnm(xi)+i,trny-rnm(yi)+j),trnx-rnm(xi)+i,trny-rnm(yi)+j,mbin(k,yidx(xi,i),yidx(yi,j))
            !pause
          end do
          if (xi.ne.yi) then
            pm(trny-rnm(yi)+j,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trny-rnm(yi)+j)
          end if
        end do
        if (xi==yi) then
          pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)
        end if

      end do
    end do

  end do

  call cpu_time (t2_cpu)
  call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu) !,LA


  xb=matmul(pm,(tyv-matmul(xm,beta)))

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
        v1=0
        trny=0
        do yi=1,tn
          trny=trny+rnm(yi)
          do j=1,rnm(yi)
            v1=v1+mbin(k,yidx(yi,j),zi)*vcva(k,xi,yi)*xb(trny-rnm(yi)+j,1)
          end do
        end do
        blup_ebv(zi,xi,k)=v1
      end do !k
    end do !zi
  end do !xi

  open (unit=43,file=fl11,status='unknown')
  write(43,*)beta
  do xi=1,tn
    do i=1,pedn
      write(43,*)xi,blup_ebv(i,xi,:)
    end do
  end do
  close(43)

end if



end subroutine

!end module
