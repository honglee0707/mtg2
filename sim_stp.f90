
subroutine sim_stp (filnam,filnam2)
 
implicit none
integer::nm,nid,fn,dn,i,j,k,k2,im,xi,yi,nl,prn,wn,nch,nmal,nfem,ich
integer,allocatable::idxmal(:),idxfem(:),sv1(:,:),nselmal(:),nselfem(:)
integer,allocatable::srtv(:,:),matpn(:,:),matcn(:,:),recmatp(:,:,:),recmatc(:,:,:)
integer,allocatable::mutpn(:,:),mutcn(:,:),mutmatp(:,:,:),mutmatc(:,:,:)
parameter(prn=15)
real::rnd
integer::ids,idd,ipar,jpar,kpar,xr,xr2,recn,iv1,recid(100),idx,k1,mutn
integer::mutid(100),trl,nmch,nmw,r1p2

double precision::mr,binpr(0:prn),mutpr(0:prn),pr,v1,v10,v2,v3
double precision,allocatable::chl(:)
integer,allocatable::nmch2(:)

character(len=100),allocatable::pedor(:,:)
character(len=20),allocatable::markt(:)
real,allocatable::mpos(:),yv(:),fst(:),lst(:)
real,allocatable::idxv(:,:),tvm(:)

integer,allocatable::ped(:,:)

character(len=1),allocatable::mafreq1(:,:)
integer*1,allocatable::genom(:,:),contom(:,:),genom1(:,:)
character(len=3)::cdum4,kchr(100)
character(len=3),allocatable::cno(:)

!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) 
integer*1 b1,igen(0:1,0:1)!,w1,w2

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,filnam2

integer,allocatable::seed(:)
integer::rseed,seed2

open (unit=37,file=trim(filnam2),status='old')             !snp list
  read(37,*)fn !# sires
  read(37,*)dn !# dams
  read(37,*)wn !# progeny
  !read(37,*)nch !# chr
  read(37,*)mr !# chr
  read(37,*)r1p2 !# random selection or phenotype selection
  read(37,*)seed2 !# random seed
  !read(37,*)ge !# generations
close(37)


! random number generator
call RANDOM_SEED(SIZE=rseed)
ALLOCATE(seed(rseed))
seed=seed2
call RANDOM_SEED(put=seed)
!dseed=seed2


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

nm=0;nch=0
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do
  nm=nm+1
  read(36,*,iostat=io)cdum4
  if (io.ne.0) exit

  do i=1,nch
    if (cdum4 == kchr(i)) goto 77
  end do
  nch=nch+1
  kchr(nch)=cdum4
77 continue
end do
close(36)
nm=nm-1
print*,'no. marker: ',nm!,kchr(1:nch)


ALLOCATE(mpos(nm),cno(nm),mafreq1(nm,2),fst(nch),lst(nch)) 
allocate(pedor(nid,4),yv(nid),sv1(nid,2),sex(nid),ped(fn*dn*wn,3))
allocate(sqgenom(int(nid/4)+1,nm),genom(nid*2,nm),genom1(2,nm),contom(fn*dn*wn*2,nm))
allocate(markt(nm),nselmal(nid),nselfem(nid))
!allocate(matpn(nch,nid*2),matcn(nch,nid*2),recmatp(nch,nid*4,1000))
!allocate(mutpn(nch,nid*2),mutcn(nch,nid*2),mutmatp(nch,nid*2,1000))
allocate(idxmal(nid),idxfem(nid),recmatc(nch,nid*4,1000),mutmatc(nch,nid*2,1000))
allocate(srtv(0,0),chl(nch),nmch2(nch))

!binomial probability work out - binpr(), mutpr() **************************
! (trl! / i! (trl-i)!) pr^i (1-pr)^(trl-i)  - binomial prob.



open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file



do i=1,nm
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
  !print*,mafreq1(i,:)
end do

  
fst=9999999999;lst=-9999999999;nmch2=0
do i=1,nm
  do j=1,nch
    if (cno(i)==kchr(j)) then
      nmch2(j)=nmch2(j)+1
      if (mpos(i) < fst(j)) fst(j)=mpos(i)
      if (mpos(i) > lst(j)) lst(j)=mpos(i)
    end if
  end do
end do

do j=1,nch
  chl(j)=(lst(j)-fst(j))*0.00000001       !100 Mbp as 1cM
end do    
!print*,chl
!print*,nmch2
!print*,fst
!print*,lst

!progressing report***************************************************
bsp(1:1)=char(8); bsp(2:2)=char(8); bsp(3:3)=char(8)

!for fam file (ped) ***************************************************
do i=1,nid
  read(35,*)pedor(i,:),sex(i),yv(i)
end do !***************************************************************


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


do im=1,nm
  if (nm<100 .or. mod(im,int(nm/100))==0) then
    write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nm)*100)
    !print*,im
!20  format (A3,I8)
  end if

  do i=1,int((nid-1)/4)+1    !no. ID / 4
    read(38)sqgenom(i,im)
  end do
  !pause

end do
write (*,'(a6)')'% done'    !closing progress report

close(35)
close(38)
close(36)


!check duplications in ped file **********************************
do i=1,nid
  do j=i+1,nid
    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
      print*,'ID duplications in .ped: ',trim(pedor(i,1)),trim(pedor(i,2))
      pause
    end if
  end do
end do !**********************************************************


!!indexed genotype coefficients into igen********************************
!igen(0,0)=2
!igen(0,1)=1
!igen(1,1)=0
!igen(1,0)=3           !missing *****************************************

nmal=0;nfem=0
do i=1,nid
  if (sex(i)==1) then
    nmal=nmal+1
    idxmal(nmal)=i
  elseif (sex(i)==2) then
    nfem=nfem+1
    idxfem(nfem)=i
  else
    print*,'no sex information, check >>>',pedor(i,1:2)
    pause
  end if
end do


PRINT*,'*********************************************************************'
PRINT*,'Structured population simualtion (e.g. livestock), given real genotype data'
print*,'# sire                  :',fn
print*,'# dam per sire          :',dn
print*,'# progeny per dam       :',wn
print*,''

  genom=0
  do i=1,nid
    if (nid<100 .or. mod(im,int(nid/100))==0) then
      write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nid)*100)
    end if

    do im=1,nm
      yi=int((i-1)/4)+1
      xi=mod((i-1),4)

      if (ibits(sqgenom(yi,im),xi*2,1)==0) then
        genom(i*2-1,im)=2
        genom(i*2,im)=2
        if (ibits(sqgenom(yi,im),xi*2+1,1)==1) genom(i*2,im)=1
      elseif (ibits(sqgenom(yi,im),xi*2,1)==1) then
        if (ibits(sqgenom(yi,im),xi*2+1,1)==1) then
          genom(i*2-1,im)=1
          genom(i*2,im)=1
        end if
      end if

    end do
  end do
  write (*,'(a6)')'% done'    !closing progress report

  deallocate(sqgenom)


  !create next generation
  do i=1,nid
    call random_number(rnd)
    sv1(i,1)=i
    if (r1p2==1) then
      sv1(i,2)=-int(rnd*10000000)
    elseif (r1p2==2) then
      sv1(i,2)=-int(yv(i)*10000000)
    end if

    !print*,rnd,int(rnd*1000000),sv1(i,2)
  end do
  call sort(sv1,nid)
  !do i=1,10
  !  print*,sv1(i,1),sv1(i,2),sex(sv1(i,1))
  !end do

  k=0
  do i=1,nid
    if (sex(sv1(i,1))==1) then
      k=k+1
      nselmal(k)=sv1(i,1)
      !print*,k,nselmal(k),sv1(i,1)
    end if
    if (k >= fn) exit
  end do
  !print*,''
  k=0
  do i=1,nid
    if (sex(sv1(i,1))==2) then
      k=k+1
      nselfem(k)=sv1(i,1)
      !print*,k,nselfem(k)
    end if
    if (k >= fn*dn) exit
  end do

  !create progeny
  idx=0
  do ipar=1,fn
    do jpar=1,dn
      do kpar=1,wn
        ids=nselmal(ipar)
        idd=nselfem(dn*ipar-dn+jpar)
        !print*,ipar,dn*ipar-dn+jpar
        idx=idx+1
        ped(idx,1)=ids
        ped(idx,2)=idd
        call random_number(rnd)
        ped(idx,3)=INT(2*rnd)+1

        !From dad**********************************************************
        do ich=1,nch

!nmch=nm/nch  !check
nmch=nmch2(ich)
trl=nmch         !trial - nmarkers/nch
pr=chl(ich)/real(nmch) !prob for occuring incidence (assuming 1M long - approx.)
!print*,pr,chl,nmch
!pause

v1=0
do i=1,trl
  v10=i
  v1=v1+log(v10)  ! trl!
end do
!print*,v1

do i=0,prn                  !incidence 0 ~ prn
  v2=0
  do j=1,i
    v10=j
    v2=v2+log(v10) ! i!
  end do
  !print*,v2

  v3=0
  do j=1,trl-i
    v10=j
    v3=v3+log(v10)  ! (trl-i)!
  end do
  !print*,v3

  binpr(i)=v1-(v2+v3)
  mutpr(i)=binpr(i)+i*log(mr)+(trl-i)*log(1-mr)
  binpr(i)=binpr(i)+i*log(pr)+(trl-i)*log(1-pr)
  !print*,binpr(i)
  !pause
end do
binpr=exp(binpr)
mutpr=exp(mutpr)
!print*,binpr(:)
!print*,mutpr(:)
!print*,mr,pr
!pause

do i=1,prn
  binpr(i)=binpr(i)+binpr(i-1)     !cumulative
  mutpr(i)=mutpr(i)+mutpr(i-1)     !cumulative
end do


          !print*,ich
          call random_number(rnd)
          xr=INT(2*rnd)+1         !random determination

          !determine how manay times rec. occur - recn *********************
          call random_number (rnd)
          do recn=0,prn    !10
            if (rnd < binpr(recn)) exit
          end do

          do i=1,recn
145         call random_number (rnd)
            iv1=int((nmch-1)*rnd)+2   !avoiding 1 
            do j=1,i-1
              if (recid(j)==iv1) goto 145    !recombined ID
            end do
            recid(i)=iv1
          end do
          !sorting************************'
          deallocate(srtv)
          allocate(srtv(recn,2))
          do i=1,recn
            srtv(i,2)=recid(i)
          end do
          call sort (srtv,recn)
          do i=1,recn
            recid(i)=srtv(i,2)
          end do !**********************************************************

          j=1     !indication for recombined marker ID
          do im=1,nmch       ! markers + QTL no.
            if (recn==0.or.j>recn) goto 501
            if (recid(j)==im) then
              xr=3-xr
              j=j+1
            end if
501         continue
            nmw=0
            do k=1,ich-1
              nmw=nmw+nmch2(k)
            end do
            genom1(1,nmw+im)=genom(ids*2-(2-xr),nmw+im)
            !print*,genom1(1,(ich-1)*nmch+im),ich,nmch,im
            !pause
          end do

          !determine how manay times mutation occur - mutn *****************
          call random_number (rnd)
          do mutn=0,prn
            if (rnd < mutpr(mutn)) exit
          end do

          do i=1,mutn       !no. mutations per chr.
143         call random_number(rnd)
            iv1=int(nmch*rnd)+1     !nmch - no. markers in one chr.
            do j=1,i-1
              if (mutid(j)==iv1) goto 143    !mutated ID
            end do
            mutid(i)=iv1
          end do

          do j=1,mutn
            k1=mutid(j)
            nmw=0
            do k=1,ich-1
              nmw=nmw+nmch2(k)
            end do
            if (genom1(1,nmw+k1)==1) then
              genom1(1,nmw+k1)=2
            elseif (genom1(1,nmw+k1)==2) then
              genom1(1,nmw+k1)=1
            end if
          end do

        end do  !ich

        !From mam**********************************************************'
        do ich=1,nch

nmch=nmch2(ich)
trl=nmch         !trial - nmarkers/nch
pr=chl(ich)/real(nmch) !prob for occuring incidence (assuming 1M long - approx.)
!print*,pr,chl,nmch
!pause

v1=0
do i=1,trl
  v10=i
  v1=v1+log(v10)  ! trl!
end do
!print*,v1

do i=0,prn                  !incidence 0 ~ prn
  v2=0
  do j=1,i
    v10=j
    v2=v2+log(v10) ! i!
  end do
  !print*,v2

  v3=0
  do j=1,trl-i
    v10=j
    v3=v3+log(v10)  ! (trl-i)!
  end do
  !print*,v3

  binpr(i)=v1-(v2+v3)
  mutpr(i)=binpr(i)+i*log(mr)+(trl-i)*log(1-mr)
  binpr(i)=binpr(i)+i*log(pr)+(trl-i)*log(1-pr)
  !print*,binpr(i)
  !pause
end do
binpr=exp(binpr)
mutpr=exp(mutpr)
!print*,binpr(:)
!print*,mutpr(:)
!print*,mr,pr
!pause

do i=1,prn
  binpr(i)=binpr(i)+binpr(i-1)     !cumulative
  mutpr(i)=mutpr(i)+mutpr(i-1)     !cumulative
end do



          !print*,'mum',ich
          call random_number(rnd)
          xr=INT(2*rnd)+1         !random determination

          !determine how manay times rec. occur - recn *********************
          call random_number (rnd)
          do recn=0,prn     !10
            if (rnd < binpr(recn)) exit
          end do

          do i=1,recn
149         call random_number (rnd)
            iv1=int((nmch-1)*rnd)+2   !avoiding 1 
            do j=1,i-1
              if (recid(j)==iv1) goto 149    !recombined ID
            end do
            recid(i)=iv1
          end do
          !sorting************************
          deallocate(srtv)
          allocate(srtv(recn,2))
          do i=1,recn
            srtv(i,2)=recid(i)
          end do
          call sort (srtv,recn)
          do i=1,recn
            recid(i)=srtv(i,2)
          end do !**********************************************************

          j=1     !indication for recombined marker ID
          do im=1,nmch       ! markers + QTL no.
            if (recn==0.or.j>recn) goto 401
            if (recid(j)==im) then
              xr=3-xr
              j=j+1
            end if
401         continue
            nmw=0
            do k=1,ich-1
              nmw=nmw+nmch2(k)
            end do
            genom1(2,nmw+im)=genom(idd*2-(2-xr),nmw+im)
          end do

          !determine how manay times mutation occur - mutn *****************
          call random_number (rnd)
          do mutn=0,prn
            if (rnd < mutpr(mutn)) exit
          end do

          do i=1,mutn       !no. mutations per chr.
147         call random_number(rnd)
            iv1=int(nmch*rnd)+1     !nmch - no. markers in one chr.
            do j=1,i-1
              if (mutid(j)==iv1) goto 147    !mutated ID
            end do
            mutid(i)=iv1
          end do

          do j=1,mutn
            k1=mutid(j)
            nmw=0
            do k=1,ich-1
              nmw=nmw+nmch2(k)
            end do
            if (genom1(2,nmw+k1)==1) then
              genom1(2,nmw+k1)=2
            elseif (genom1(2,nmw+k1)==2) then
              genom1(2,nmw+k1)=1
            end if
          end do
        end do !ich
        !print*,genom1(2,:)

        contom(idx*2-1,:)=genom1(1,:)   !paternal genotype
        contom(idx*2,:)=genom1(2,:)     !maternal genotype

      end do
    end do
  end do


  open (unit=101,file=trim(filnam2)//'.ped',status='unknown')
  do k=1,idx
    write(101,'(4I10,2I3,100000000I2)')k,k,ped(k,1),ped(k,2),ped(k,3),-9 &
&   ,(contom(k*2-1,j),contom(k*2,j),j=1,nm)

  end do
  close(101)

  open (unit=102,file=trim(filnam2)//'.map',status='unknown')
  do i=1,nm
    !write(102,'(*(G0,2X))')cno(i),markt(i),0,mpos(i)*100/100
    write(102,'(a3,a20,2f14.0)')cno(i),markt(i),0,mpos(i)
  end do
  close(102)


end subroutine

