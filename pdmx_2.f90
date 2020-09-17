
subroutine pdmx_ (nid,pedor,filnam,grm,mn,f41)
 
implicit none
integer::nid,i,j,k,idum,ncol,mn
real::grm(mn,nid,nid),v1
character(len=100)::pedor(nid,4),cdum11,cdum12
real,allocatable::rtmx(:,:),kvm(:,:),kvm2(:,:),mat(:,:)
integer,allocatable::rid(:)

character(len=3)::cdum4

!progressing report **************************
character*3 bsp

!reading command line *****************
integer::narg,io
character(len=64)::filnam,cdum,f41
character(len=7000000)::string    !check size


open (unit=35,file=trim(filnam),status='old')      !pdmx file
read(35,'(a)')string
close(35)
!print*,string

do i=1,7000000
  read(string,*,iostat=io)(cdum4,j=1,i)
  if (io.ne.0) exit
end do
if (i == 7000000) then
  print*,"# column exceed >>> check",j-1-2
  pause
end if

ncol=j-1-2 ! (2 for FID IID)
print*,"no. column: ",ncol !,mn

allocate(kvm(nid,ncol),kvm2(nid,ncol),mat(ncol,ncol))
allocate(rtmx(nid,nid),rid(ncol))

open (unit=35,file=trim(filnam),status='old')      !fam file
do i=1,nid
  read(35,*)cdum11,cdum12,(kvm(i,j),j=1,ncol)
  if (cdum11.ne.pedor(i,1) .or. cdum12.ne.pedor(i,2)) then
    print*,"FID IID in .fam should be matched with FID IID in ",trim(filnam)
    print*,trim(cdum11),pedor(i,1),trim(cdum12),pedor(i,2)
    pause
  end if
end do
close(35)


if (mn==0) then
  PRINT*,'*********************************************************************'
  print*,"product matrix (WW') to be fit as random effects (like GRM)"

  !rtmx=matmul(kvm,transpose(kvm))

  !$OMP PARALLEL DO PRIVATE(i,j,k,v1)
  do i=1,nid
    do j=1,nid
      v1=0
      do k=1,ncol
        v1=v1+kvm(i,k)*kvm(j,k)
      end do
      rtmx(i,j)=v1
    end do
  end do
  !$OMP END PARALLEL DO 

  deallocate(kvm)
  print*,'done' !*********************************************************

elseif (mn==1) then

  PRINT*,'*********************************************************************'
  print*,"product matrix (WGW') to be fit as random effects (like GRM)"
  print*,"G scaled to a dimension of",ncol

  open (unit=36,file=trim(filnam)//'.id',status='old')      !fam file
  do i=1,ncol
    read(36,*)rid(i)
  end do
  close(36)
  
  do i=1,ncol
    do j=1,ncol
      mat(i,j)=grm(1,rid(i),rid(j))
    end do
  end do
  !print*,mat


  !$OMP PARALLEL DO PRIVATE(i,j,k,v1)
  do i=1,nid
    do j=1,ncol
      v1=0
      do k=1,ncol
        v1=v1+kvm(i,k)*mat(j,k)
      end do
      kvm2(i,j)=v1
    end do
  end do
  !$OMP END PARALLEL DO 

  !$OMP PARALLEL DO PRIVATE(i,j,k,v1)
  do i=1,nid
    do j=1,nid
      v1=0
      do k=1,ncol
        v1=v1+kvm2(i,k)*kvm(j,k)
      end do
      rtmx(i,j)=v1
    end do
  end do
  !$OMP END PARALLEL DO 

  deallocate(kvm,kvm2)
  print*,'done' !*********************************************************

end if


open (UNIT=41,FILE=trim(f41),STATUS='unknown')
do i=1,nid
  do j=1,i 
    if (rtmx(i,j).ne.0.0D0) then
      write(41,*)i,j,rtmx(i,j)
    end if
  end do
end do
close(41)

deallocate(rtmx)

end subroutine

