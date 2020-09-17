
subroutine aireml_m_gwas2 (mbin,yidx,yv,rn,rnm,trn,pedn,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res,filnam)        

!***********************************************************************
!aireml_b: multivariate analysis especially for case control 
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
!double precision::mbin(mn,pedn,pedn),yv(rn,2),tyv(rn*2,1)
real::mbin(mn,pedn,pedn),yv(rn,tn),tyv(trn,1),sv_fact
double precision::resw(trn),resw2(pedn,tn)

!double precision::qq(mn,rn,rn),qq2(mn,rn,rn)
!real::qq(mn,rn,rn),qq2(mn,rn,rn)

double precision::va(mn),va2(mn),cov(mn),ve,ve2,LKH,ypyf(1,1),tmpm(1,1),conv
double precision::vcva(mn,tn,tn),vcve(tn,tn),vvc(tn)
double precision::fix_vcva(mn,tn,tn),fix_vcve(tn,tn)
double precision::vcva_rsv(mn,tn,tn),vcve_rsv(tn,tn)
double precision::blup_ebv(pedn,tn,mn),blup_py(pedn,tn,mn),beta(tfn,1) !,cor
double precision::blup_r(pedn,tn,mn)

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
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

character(len=7)::cdum
character(len=64)::fl4,fl24,fl5,fl11,fl11r,wt_res,cdum_t(tn),filnam
logical::file_exist

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) !,igenom(:,:)
integer*1 b1,igen(0:1,0:1),w1,w2

integer::nid,nsnp,im,gwasi
integer,allocatable::snpv(:),nmiss(:)   
character(len=20),allocatable::mnam(:),markt(:),cmk(:)
real,allocatable::mpos(:),raf(:)

character(len=1),allocatable::mafreq1(:,:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which


!progressing report **************************
character*3 bsp
!character*1 creturn



nid=0
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
do
  nid=nid+1
  read(35,*,iostat=io)cdum4
  if (io.ne.0) exit
  !print*,trim(pedor(i,1)),trim(pedor(i,2))
end do
close(35)
nid=nid-1
print*,'no. ID    : ',nid


nsnp=0
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do
  nsnp=nsnp+1
  read(36,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(36)
nsnp=nsnp-1
print*,'no. marker: ',nsnp

ALLOCATE(mpos(nsnp),cno(nsnp),mnam(nsnp),mafreq1(nsnp,2))
!allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nsnp),sex(nid))
allocate(markt(nsnp),snpv(nid),raf(nsnp),nmiss(nsnp))


open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file


do i=1,nsnp
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
  !print*,mafreq1(i,:)
end do

!progressing report***************************************************
bsp(1:1)=char(8); bsp(2:2)=char(8); bsp(3:3)=char(8)

!!for fam file (ped) ***************************************************
!do i=1,nid
!  read(35,*)pedor(i,:),sex(i),cdum4
!end do !***************************************************************

! for pedigree and genotype (PLINK .bed, * files) ****************************
read(38)b1           !plink magic number 1

if (.not.btest(b1,0).and..not.btest(b1,1).and.btest(b1,2).and.btest(b1,3).and.&
& .not.btest(b1,4).and.btest(b1,5).and.btest(b1,6).and..not.btest(b1,7)) then
  write(*,'(a7)',advance='no')' 1 - ok'
else
  print*,'this may not be PLINK .bed file - check >>>'
end if

read(38)b1           !plink magic number 2
if (btest(b1,0).and.btest(b1,1).and..not.btest(b1,2).and.btest(b1,3).and.&
& btest(b1,4).and..not.btest(b1,5).and..not.btest(b1,6).and..not.btest(b1,7)) then
  write(*,'(a7)',advance='no')' 2 - ok'
else
  print*,'this may not be PLINK .bed file - check >>>'
end if

read(38)b1           !mode 

if (btest(b1,0)) then
  write(*,'(a37)')'  SNP-major mode for PLINK .bed file'
else
  print*,'should not be individual mode - check >>>'
end if

do im=1,nsnp
  if (nsnp<100 .or. mod(im,int(nsnp/100))==0) then
    write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nsnp)*100)
    !print*,im
!20  format (A3,I8)
  end if

  do i=1,int((nid-1)/4)+1    !no. ID / 4
    read(38)sqgenom(i,im)
  end do
  !print*,igenom(i,1:20)
  !pause

end do
write (*,'(a6)')'% done'    !closing progress report

close(35)
close(38)
close(36)

!!check duplications in ped file **********************************
!do i=1,nid
!  do j=i+1,nid
!    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
!      print*,'ID duplications in .ped: ',trim(pedor(i,1)),trim(pedor(i,2))
!      pause
!    end if
!  end do
!end do !**********************************************************


!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************




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
tfn=tfn+tn
print*,tfn

print*,'*** number of records used ***'
do i=1,tn
  print*,'trait',i,':',rnm(i)
end do
print*,''

if (wt_res.ne."null") then
  open (unit=50,file=wt_res,status='old')
  do i=1,pedn
    read(50,*)cdum,cdum,cdum_t(:)
    do j=1,tn
      if (cdum_t(j).ne.'NA') then
        read(cdum_t(j),*)resw2(i,j)
      else
        resw2(i,j)=-99999   !missing
      end if
      !print*,resw2(i,:)
    end do
  end do
  close(50)
  !pause

  k=0
  do j=1,tn
    do i=1,pedn
      if (resw2(i,j).ne.-99999) then
        k=k+1
        resw(k)=1/resw2(i,j)
        !print*,resw(k)
      end if
    end do
    if (k .ne. sum(rnm(1:j))) then
      print*,'# residual weights is not matched with # records'
      pause
    end if
  end do

else      !if null
  resw=1.0D0
end if




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

  ! y = Xb + Zu + e ---- fitting effects without QTL
  ! V matrix    V = ZAZ'Va + IVe *************************

  call cpu_time (t2_cpu)

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
          !pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)
          pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)=pm(trnx-rnm(xi)+i,trnx-rnm(xi)+i)+vcve(xi,xi)*resw(trnx-rnm(xi)+i)
        end if

      end do
    end do

  end do

  call cpu_time (t2_cpu)
  call cholesky_inv (pm,trn,LA)
  call cpu_time (t1_cpu)
  print*,'V inverse done - time:',real(t1_cpu-t2_cpu)!,LA
  print*,''


open (unit=47,file=trim(fl5)//".gwas",status='unknown')
write(47,'(a3,a17,a3,a8,2a14,a12,a14)')'CHR','SNP','A1','N','BETA','SE','CHI','P'
do gwasi=1,nsnp
  !print'(a4,i8,5x,f5.2,a12)','SNP #',gwasi,real(gwasi)/real(nsnp)*100,"% progressed"
  print'(a4,i8,5x,f5.2,a12)','SNP #',gwasi,real(gwasi)/real(nsnp)*100,"% progressed ***"
  !creturn=achar(13)
  !write (*,'(a,a4,i8,5x,f5.2,a12)',advance='no')creturn,'SNP #',gwasi,real(gwasi)/real(nsnp)*100,"% progressed ***"
  x1=0;x2=0
  do i=1,nid
    yi=int((i-1)/4)+1
    xi=mod((i-1),4)
    w1=igen(ibits(sqgenom(yi,gwasi),xi*2,1),ibits(sqgenom(yi,gwasi),xi*2+1,1))
    snpv(i)=w1
    if (w1.ne.3) then
      x1=x1+w1; x2=x2+1
    end if
  end do
  raf(gwasi)=(x1/x2)*0.5
  nmiss(gwasi)=x2
  do i=1,nid
    if (snpv(i)==3) snpv(i)=x1/x2
  end do

trn=0
do xi=1,tn
  do i=1,rn
    if (yv(i,xi).ne.-99999) then
      trn=trn+1
      xm(trn,tfn-tn+xi)=snpv(yidx(xi,trn))
      !print*,xm(rn1,1:fn1),xm1(i,:),rn1,i
      !pause
    end if
  end do
end do


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


!CDF functions *****************************************************
      !p = 1D+00 - k !0.011135D+00
      !q = 1.0D+00 - p
      p = huge ( p )
      q = huge ( q )
      x = beta(tfn,1)**2/xvmx(tfn,tfn)
      !x=-1.28155 !2 3.0D+00
      !mean = 0.0D+00 ! 5.0D+00
      sd = 1.0D+00 !0.875D+00    !degree of greedom for cdfchi
which=1  !Calculate P and Q from X, MEAN and SD;
!which=2  !Calculate X from P, Q, MEAN and SD;
!which=3  !Calculate MEAN from P, Q, X and SD;
!which=4  !Calculate SD from P, Q, X and MEAN;

call cdfchi ( which, p, q, x, sd, status, bound )


  !do xi=1,tfn
  !  write(47,'(i3,100f24.16)')xi,beta(xi,1),sqrt(xvmx(xi,xi))
  !end do
  write(47,'(a3,a17,a3,i8,2f14.5,f12.3,e14.5)')cno(gwasi),trim(markt(gwasi)),mafreq1(gwasi,1),nmiss(gwasi),real(beta(tfn,1)),real(sqrt(xvmx(tfn,tfn))),real(beta(tfn,1)**2/xvmx(tfn,tfn)),real(q)




end do !gwasi


close(47)


end subroutine

!end module
