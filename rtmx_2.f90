
subroutine rtmx_ (filnam,sclf,wt_par,f41)
 
implicit none
real::sclf     ! scale factor
integer::nm,nid,i,j,k,k1,k2,zi,zj,im,xi,yi
integer*1,allocatable::igenom(:,:)

character(len=100),allocatable::pedor(:,:)
character(len=50),allocatable::markt(:)
character(len=16)::cfst,csnd
real,allocatable::mafreq2(:),mpos(:),htz(:)
real,allocatable::wt(:,:)
!real,allocatable::kvm(:,:),kvm2(:,:)
double precision,allocatable::kvm(:,:),kvm2(:,:),rtmx(:,:),rtmx2(:,:)

character(len=1),allocatable::genom(:,:),mafreq1(:,:)
character(len=3)::cdum4
character(len=3),allocatable::cno(:)

integer::nsd,fst,snd,sdv,traitn,trait_no(10)    !check 10 is alway enough
real::x1,x2,x3,rij2,wij,wij2
integer,allocatable::sdt(:,:),idx(:)


!sequential storing
integer*1,allocatable::sqgenom(:,:),sqgenom_1(:,:),sqgenom_2(:,:),sex(:) 
integer*1 b1,igen(0:1,0:1),w1,w2

logical::wt_logic

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,cdum,f41,wt_par,wt_file





!narg=command_argument_count()
!if (narg.ne.4) then
!  print*,'1. name for plink bed format, 2. no. sections, 3. first section, 4. second section  required ...'
!  stop
!end if

!call get_command_argument(1,filnam)
!print*,'plink file name     :',trim(filnam)
!call get_command_argument(2,cdum)
!read(cdum,*)nsd
!print*,'number of section :',nsd
!call get_command_argument(3,cdum)
!read(cdum,*)fst
!print*,'first section     :',fst
!call get_command_argument(4,cdum)
!read(cdum,*)snd
!print*,'first section     :',snd

!if (snd>fst) then
!  print*,'1st block should larger than 2nd block'
!  pause
!end if

nsd=1;fst=1;snd=1            !check, later to modify

nid=0
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
do
  nid=nid+1
  read(35,*,iostat=io)cdum4
  if (io.ne.0) exit
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


ALLOCATE(mpos(nm),cno(nm),mafreq1(nm,2),mafreq2(nm)) 
allocate(pedor(nid,4))
allocate(sqgenom(1,nm),sex(nid))
allocate(htz(nm))
allocate(markt(nm),sdt(nsd,2),idx(nid))


wt_file=trim(filnam)//".wt"
traitn=1
idx=1
if (wt_par.ne."null") then
  open (unit=51,file=wt_par,status='old')
  read(51,*)wt_file
  read(51,*)traitn
  k1=0
  do i=1,traitn
    read(51,*)trait_no(i)
    idx(k1+1:trait_no(i))=i
    k1=trait_no(i)
  end do
end if

allocate(wt(nm,traitn))


!do i=1,nid
!  print*,i,idx(i)
!end do

!pop. division ***************************************************
sdv=int(nid/real(nsd))    !value for dividing pop.

sdt(1,1)=1
sdt(1,2)=sdt(1,1)+sdv
do i=2,nsd
  sdt(i,1)=sdt(i-1,2)+1
  sdt(i,2)=sdt(i,1)+sdv
  if (sdt(i,1).ge.nid) then
    print*,'nsd should be less than ',i-1
    pause
  end if
end do
if (sdt(nsd,2).ge.nid) then
  sdt(nsd,2)=nid
else
  print*,'check',nsd,sdt(nsd,2),nid
end if
!*****************************************************************

allocate(kvm(sdv+1,nm),kvm2(nm,sdv+1))
allocate(rtmx(sdv+1,sdv+1),rtmx2(sdv+1,sdv+1))

allocate(sqgenom_1(int(sdv/4)+10,nm),sqgenom_2(int(sdv/4)+10,nm))   ! +10 make sure enough memory


open (unit=38,file=trim(filnam)//'.bed',status='old',access='stream',form='unformatted')
open (unit=35,file=trim(filnam)//".fam",status='old')      !fam file
open (unit=36,file=trim(filnam)//".bim",status='old')      !fam file


do i=1,nm
  read(36,*)cno(i),markt(i),cdum4,mpos(i),mafreq1(i,1),mafreq1(i,2)
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

  k1=0;k2=0
  do i=1,int((nid-1)/4)+1    !no. ID / 4
    read(38)sqgenom(1,im)

    if (i>=int((sdt(fst,1)-1)/4)+1 .and. i<=int((sdt(fst,2)-1)/4)+1 ) then
      k1=k1+1
      sqgenom_1(k1,im)=sqgenom(1,im)
    end if
    if (i>=int((sdt(snd,1)-1)/4)+1 .and. i<=int((sdt(snd,2)-1)/4)+1 ) then
      k2=k2+1
      sqgenom_2(k2,im)=sqgenom(1,im)
    end if

  end do
  !print*,igenom(i,1:20)
  !pause

end do
write (*,'(a6)')'% done'    !closing progress report

close(35)
close(38)
close(36)

!print*,k1,k2



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
PRINT*,'allele frequency read from ',trim(filnam)//".freq"

open (unit=37,file=trim(filnam)//".freq",status='old')      !fam file
  do im=1,nm

    if (nm<100 .or. mod(im,int(nm/100))==0) then
      write (*,'(a3,i3)',advance='no')bsp,int((.1*im)/(.1*nm)*100)
    end if

    read(37,*)cdum4,cdum4,mafreq2(im),htz(im)
  end do
  write (*,'(a6)')'% done'    !closing progress report
close(37)

wt=1
!inquire(file=trim(filnam)//".wt",exist=wt_logic)
inquire(file=wt_file,exist=wt_logic)
if (wt_logic) then
  open (unit=39,file=wt_file,status='old')      !weighting file
  do im=1,nm
    read(39,*)cdum4,cdum4,(wt(im,k1),k1=1,traitn)
  end do
  print*,'wt file:',wt_file
end if


deallocate(sqgenom)   

PRINT*,'*********************************************************************'
print*,'relatedness estimation - scale factor:',sclf

deallocate(kvm,kvm2)
allocate(kvm(sdv+1,nm),kvm2(nm,sdv+1))

!sum of heterozygosity *****************************************************
!$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    do i=sdt(fst,1),sdt(fst,2)        !1st pop.
      yi=(int((i-1)/4)+1) - (int((sdt(fst,1)-1)/4)+1) +1 
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom_1(yi,im),xi*2,1),ibits(sqgenom_1(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        !kvm(i-sdt(fst,1)+1,im)=1   !for counting (e.g. 1.6)  **********
        !kvm(i-sdt(fst,1)+1,im)=wt(im)   !for counting (e.g. 1.6)  **********
        kvm(i-sdt(fst,1)+1,im)=wt(im,idx(i))  
      else
        kvm(i-sdt(fst,1)+1,im)=0
      end if
    end do
  end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    do i=sdt(snd,1),sdt(snd,2)        !2nd pop 
      yi=(int((i-1)/4)+1) - (int((sdt(snd,1)-1)/4)+1) +1 
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom_2(yi,im),xi*2,1),ibits(sqgenom_2(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        !kvm2(im,i-sdt(snd,1)+1)=1    !for counting (1,6) **********
        !kvm2(im,i-sdt(snd,1)+1)=wt(im)    !for counting (1,6) **********
        kvm2(im,i-sdt(snd,1)+1)=wt(im,idx(i)) 
      else
        kvm2(im,i-sdt(snd,1)+1)=0
      end if
    end do
  end do
!$OMP END PARALLEL DO

  !print*,'dgemm starts'
  !rtmx2=matmul(kvm,kvm2)
  call dgemm ('N','N',(sdv+1),(sdv+1),nm,1.0D0,kvm,(sdv+1),kvm2,nm,0.0D0,rtmx2,(sdv+1))

print*,'50% done' !*********************************************************

!sum of product ************************************************************
!$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    do i=sdt(fst,1),sdt(fst,2)        !1st pop.
      yi=(int((i-1)/4)+1) - (int((sdt(fst,1)-1)/4)+1) +1 
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom_1(yi,im),xi*2,1),ibits(sqgenom_1(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        !kvm(i-sdt(fst,1)+1,im)=(w1-mafreq2(im))/htz(im)**0.5
        !kvm(i-sdt(fst,1)+1,im)=(w1-mafreq2(im))*(htz(im)**0.5)**sclf
        !kvm(i-sdt(fst,1)+1,im)=wt(im)*(w1-mafreq2(im))*(htz(im)**0.5)**sclf
        kvm(i-sdt(fst,1)+1,im)=wt(im,idx(i))*(w1-mafreq2(im))*(htz(im)**0.5)**sclf
        !print*,i,sdt(fst,1),i-sdt(fst,1)+1,idx(i),wt(im,idx(i)),im
      else
        kvm(i-sdt(fst,1)+1,im)=0
      end if
    end do
  end do
  !pause
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE (yi,xi,w1) 
  do im=1,nm
    do i=sdt(snd,1),sdt(snd,2)        !2nd pop 
      yi=(int((i-1)/4)+1) - (int((sdt(snd,1)-1)/4)+1) +1 
      xi=mod((i-1),4)
      w1=igen(ibits(sqgenom_2(yi,im),xi*2,1),ibits(sqgenom_2(yi,im),xi*2+1,1))
      if (w1.ne.3) then
        !kvm2(im,i-sdt(snd,1)+1)=(w1-mafreq2(im))/htz(im)**0.5
        !kvm2(im,i-sdt(snd,1)+1)=(w1-mafreq2(im))*(htz(im)**0.5)**sclf
        !kvm2(im,i-sdt(snd,1)+1)=wt(im)*(w1-mafreq2(im))*(htz(im)**0.5)**sclf
        kvm2(im,i-sdt(snd,1)+1)=wt(im,idx(i))*(w1-mafreq2(im))*(htz(im)**0.5)**sclf
        !print*,i,sdt(snd,1),i-sdt(snd,1)+1,idx(i),wt(im,idx(i)),im
      else
        kvm2(im,i-sdt(snd,1)+1)=0
      end if
    end do
  end do
  !pause
!$OMP END PARALLEL DO

  !print*,'dgemm starts'
  !rtmx=matmul(kvm,kvm2)
  call dgemm ('N','N',(sdv+1),(sdv+1),nm,1.0D0,kvm,(sdv+1),kvm2,nm,0.0D0,rtmx,(sdv+1))

  deallocate(kvm,kvm2)

print*,'100% done' !*********************************************************

!write(cfst,*)fst
!write(csnd,*)snd

cfst=adjustl(cfst)
csnd=adjustl(csnd)

if (f41=="ascm.out") then
  !f41=trim(filnam)//".rtmx_"//trim(cfst)//"_"//trim(csnd)
  f41=trim(filnam)//".rtmx"      !check, later to modify
end if

open (UNIT=41,FILE=trim(f41),STATUS='unknown')
do i=sdt(fst,1),sdt(fst,2)        !1st pop.
  do j=sdt(snd,1),sdt(snd,2)        !2nd pop 
    if (j.lt.i) then
      write(41,*)i,j,real(rtmx(i-sdt(fst,1)+1,j-sdt(snd,1)+1)/rtmx2(i-sdt(fst,1)+1,j-sdt(snd,1)+1))
    elseif (j==i) then
      write(41,*)i,j,real(rtmx(i-sdt(fst,1)+1,j-sdt(snd,1)+1)/rtmx2(i-sdt(fst,1)+1,j-sdt(snd,1)+1))
    end if
  end do
end do
close(41)



end subroutine

