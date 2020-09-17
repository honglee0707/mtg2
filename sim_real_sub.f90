
subroutine sim_real_sub (filnam,filnam2)
 
implicit none
integer::nm,nid,sn,dn,i,j,k,k2,zi,zj,im,ma1,idum,nm2,nm3,xi,yi,xi2,yi2,nl
integer,allocatable::alf(:,:),icmk(:)
integer*1,allocatable::igenom(:,:)
!real,allocatable::igenom(:,:)

character(len=100),allocatable::pedor(:,:)
character(len=20),allocatable::mnam(:),markt(:),cmk(:)
real,allocatable::mpos(:)
real,allocatable::idxv(:,:),tvm(:)
real,allocatable::bv(:),alpha(:)

character(len=1),allocatable::genom(:,:),mafreq1(:,:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

integer::ip,jp,icno,optn

real::rdum,x1,x2,x3,rij2,wij,wij2,fra,frb,frc,frd,sac,sab,sad,sbc,sbd,scd
real::rij,ssm,ssm2,v1


!sequential storing
integer*1,allocatable::sqgenom(:,:),sex(:) !,igenom(:,:)
integer*1 b1,igen(0:1,0:1),w1,w2
integer::iv1

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,filnam2

!narg=command_argument_count()
!if (narg.ne.2) then
!  print*,'name for plink bed format, snp list required ...'
!  stop
!end if

!  call get_command_argument(1,filnam)
!  print*,'plink file name:',trim(filnam)
!  call get_command_argument(2,filnam2)
!  print*,'snp list       :',trim(filnam2)

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

nm=0
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
do
  nm=nm+1
  read(36,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(36)
nm=nm-1
print*,'no. marker: ',nm


nl=0
open (unit=37,file=trim(filnam2),status='old')      !fam file
do
  nl=nl+1
  read(37,*,iostat=io)cdum4
  if (io.ne.0) exit
end do
close(37)
nl=nl-1
print*,'no. causal marker: ',nl


ALLOCATE(mpos(nm),cno(nm),mnam(nm),mafreq1(nm,2)) 
allocate(pedor(nid,4))
allocate(sqgenom(int(nid/4)+1,nm),sex(nid))
allocate(bv(nid),alpha(nl))
allocate(markt(nm))
allocate(icmk(nl),cmk(nl))


open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file
open (unit=37,file=trim(filnam2),status='old')             !snp list




do i=1,nm
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
  !print*,mafreq1(i,:)
end do

do i=1,nl
  read(37,*)icmk(i),cmk(i),alpha(i)  !,mafreq1(i,1),mafreq1(i,2)
  !print*,mafreq1(i,:)
end do

!progressing report***************************************************
bsp(1:1)=char(8); bsp(2:2)=char(8); bsp(3:3)=char(8)

!for fam file (ped) ***************************************************
do i=1,nid
  read(35,*)pedor(i,:),sex(i),cdum4
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
  !print*,igenom(i,1:20)
  !pause

end do
write (*,'(a6)')'% done'    !closing progress report

close(35)
close(38)
close(36)
close(37)


!check duplications in ped file **********************************
do i=1,nid
  do j=i+1,nid
    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
      print*,'ID duplications in .ped: ',trim(pedor(i,1)),trim(pedor(i,2))
      pause
    end if
  end do
end do !**********************************************************


!indexed genotype coefficients into igen********************************
igen(0,0)=2
igen(0,1)=1
igen(1,1)=0
igen(1,0)=3           !missing *****************************************


PRINT*,'*********************************************************************'
PRINT*,'breeding values'

  bv=0;x1=0;x2=0
  do i=1,nid

    if (nid<100 .or. mod(im,int(nid/100))==0) then
      write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nid)*100)
    end if

    do im=1,nl

      yi=int((i-1)/4)+1
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom(yi,icmk(im)),xi*2,1),ibits(sqgenom(yi,icmk(im)),xi*2+1,1))
      if (w1.ne.3) then
        bv(i)=bv(i)+w1*alpha(im)
      end if
    end do

    x1=x1+bv(i)
    x2=x2+bv(i)**2
  end do
  write (*,'(a6)')'% done'    !closing progress report

  open (UNIT=41,FILE=trim(filnam)//'.bv',STATUS='unknown')
  do i=1,nid
    !write(41,*)bv(i)
    write(41,*)(bv(i)-x1/nid)/sqrt((x2-(x1**2/nid))/(nid-1))
  end do
  close(41)


end subroutine

