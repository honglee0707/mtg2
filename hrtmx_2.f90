subroutine hrtmx_ (hrtmx,grm,pedn,n,mn,fam)

!*********************************************************************
!To get H matrix for SSGBLUP
!S. Hong Lee (2017)
!*********************************************************************

implicit none

integer::n,pedn,i,j,w,k,kk,k1,k2,fn,mn,zi
character(len=64)::hrtmx
character(len=100)::ped(pedn,4),idx_mem(pedn*5),fam(n,4)  !check
real::grm(mn,n,n),nrm(pedn,pedn)
integer::idx_ped(pedn*5,4),idx_pedn,idx_fam2(n) !x5 check
integer::idx(pedn),idx_fam(n),idx_01(pedn),idx_total(pedn),idx_total2(pedn)
real,allocatable::idx_nrm(:,:)

double precision::A11((pedn-n),(pedn-n)),A21(n,(pedn-n)),GA22(n,n)
double precision::A22(n,n)
integer,allocatable::idx_ped2(:,:)
double precision::x1

double precision::A22A21(n,(pedn-n)),A12A22((pedn-n),n)
double precision::GA22A21(n,(pedn-n)),G22(n,n)
double precision::H11_tmp((pedn-n),n),H11_tmp2((pedn-n),(pedn-n))
double precision::H11_tmp4((pedn-n),n),H11_tmp3((pedn-n),(pedn-n))

!print*,pedn,n,mn
open (unit=48,file=hrtmx,status='old') 
do i=1,pedn
  read(48,*)ped(i,:)
end do
close(48)

do zi=1,pedn
  do i=1,pedn
    if ((ped(zi,3).ne.'0' .and. ped(zi,3) == ped(i,4)) .or. (ped(zi,4).ne.'0' .and. ped(zi,4) == ped(i,3))) then
      print*,'pedigree wrong >>> check',zi,i
      pause
    end if
  end do
end do



idx_pedn=0
do zi=1,pedn

  if (ped(zi,3) == '0') goto 10
  do i=1,idx_pedn
    if (ped(zi,3) == idx_mem(i)) goto 10
  end do

  do i=1,pedn
    if (ped(zi,3) == ped(i,2)) then
      if (i >= zi) then
        print*,'check >>> offspring is before parents'
        pause
      end if
      goto 10 
    end if
  end do

  idx_pedn=idx_pedn+1
  idx_mem(idx_pedn)=ped(zi,3)
  idx_ped(idx_pedn,2)=idx_pedn
  idx_ped(idx_pedn,3)=0
  idx_ped(idx_pedn,4)=0

10 continue

  if (ped(zi,4) == '0') goto 20
  do i=1,idx_pedn
    if (ped(zi,4) == idx_mem(i)) goto 20
  end do

  do i=1,pedn
    if (ped(zi,4) == ped(i,2)) then 
      if (i >= zi) then
        print*,'check >>> offspring is before parents'
        pause
      end if
      goto 20 
    end if
  end do

  idx_pedn=idx_pedn+1
  idx_mem(idx_pedn)=ped(zi,4)
  idx_ped(idx_pedn,2)=idx_pedn
  idx_ped(idx_pedn,3)=0
  idx_ped(idx_pedn,4)=0

20 continue

end do !zi

do zi=1,pedn
  idx_pedn=idx_pedn+1
  idx_mem(idx_pedn)=ped(zi,2)
  idx_ped(idx_pedn,2)=idx_pedn
  if (ped(zi,3)=='0') then
    idx_ped(idx_pedn,3)=0
    goto 30
  end if
  do i=1,(idx_pedn-1)
    if (idx_mem(i)==ped(zi,3)) then
      idx_ped(idx_pedn,3)=i
      goto 30
    end if
  end do
  print*,'something wrong >>> check',zi
  pause 
30 continue

  if (ped(zi,4)=='0') then
    idx_ped(idx_pedn,4)=0
    goto 40
  end if
  do i=1,(idx_pedn-1)
    if (idx_mem(i)==ped(zi,4)) then
      idx_ped(idx_pedn,4)=i
      goto 40
    end if
  end do
  print*,'something wrong >>> check',zi
  pause 
40 continue

end do !zi


!indexing
do zi=1,pedn
  do i=1,idx_pedn
    if (idx_mem(i)==ped(zi,2)) then
      idx(zi)=i
      goto 50
    end if
  end do
50 continue
end do

!open (unit=51,file='idx_ped',status='unknown')
!do i=1,idx_pedn
!  !print*,i,trim(idx_mem(i)),idx_ped(i,2:4)
!  write(51,*)i,trim(idx_mem(i)),idx_ped(i,2:4)
!end do
!close(51)
!do i=1,pedn
!  print*,idx(i)
!end do

allocate(idx_nrm(idx_pedn,idx_pedn),idx_ped2(idx_pedn,4))
idx_ped2=idx_ped(1:idx_pedn,:)

!do i=1,idx_pedn
!  print*,idx_ped2(i,2:4)
!end do

call Amat(idx_nrm,idx_ped2,idx_pedn)

do zi=1,pedn
  do j=1,pedn
    nrm(zi,j)=idx_nrm(idx(zi),idx(j))
    !print*,nrm(zi,j)
  end do
end do

deallocate(idx_nrm,idx_ped2)

if (mn > 0) then  !matching genotype info
  if (mn > 1) then
    print*,'# GRM should be not more than 1'
    pause
  end if

  do zi=1,n
    !print*,zi
    do i=1,pedn
      !print*,i
      if (fam(zi,2) == ped(i,2)) goto 60
    end do
    print*,'no matching between pedigree and fam >>> check',zi
    pause
60  continue
    idx_fam(zi)=i
  end do
  
  idx_01=0
  do i=1,n
    !print*,grm(1,i,i)
    !print*,idx_fam(i)
    idx_01(idx_fam(i))=1
  end do
  !print*,idx_01

  k=0
  do i=1,pedn
    if (idx_01(i)==0) then
      k=k+1
      idx_total(k)=i
      idx_total2(i)=k        !for getting back to the original order
    end if
  end do
  k1=k
  do i=1,pedn
    if (idx_01(i)==1) then
      k=k+1
      idx_total(k)=i
      idx_total2(i)=k        !for getting back to the original order
    end if
  end do
  k2=k-k1
  !print*,idx_total
  !pause
  !print*,idx_total2

  if (k.ne.pedn .or. k1.ne.pedn-n .or. k2.ne.n) then
    print*,'number not mached',k,k1,k2,pedn,n
    pause
  end if

  !GRM order index for G22
  do i=1,n
    do j=1,n
      if (ped(idx_total(k1+i),2)==fam(j,2)) then
        idx_fam2(i)=j
        goto 70
      end if
    end do
    print*,'something wrong >>> check'
70 continue
  end do
  !pause
  !print*,idx_fam2

  !A11
  do i=1,k1
    do j=1,k1
      A11(i,j)=nrm(idx_total(i),idx_total(j))
    end do
  end do
  !A21
  do i=1,k2
    do j=1,k1
      A21(i,j)=nrm(idx_total(k1+i),idx_total(j))
    end do
  end do
  !A22
  do i=1,k2
    do j=1,k2
      A22(i,j)=nrm(idx_total(k1+i),idx_total(k1+j))
    end do
  end do
  
  !G22=grm(1,:,:)
  do i=1,k2
    do j=1,k2
      G22(i,j)=grm(1,idx_fam2(i),idx_fam2(j))
    end do
  end do

  !G-A22
  GA22=G22-A22
  !print*,GA22

  !A22-1
  call cholesky_inv (A22,k2,x1)
  !print*,A22

  call dgemm ('N','N',n,(pedn-n),n,1.0D0,A22,n,A21,n,0.0D0,A22A21,n)
  do i=1,(pedn-n)
    do j=1,n
      A12A22(i,j)=A22A21(j,i)
    end do
  end do
  !print*,A12A22
  !pause
  call dgemm ('N','N',(pedn-n),n,n,1.0D0,A12A22,(pedn-n),GA22,n,0.0D0,H11_tmp,(pedn-n))
  call dgemm ('N','N',(pedn-n),(pedn-n),n,1.0D0,H11_tmp,(pedn-n),A22A21,n,0.0D0,H11_tmp2,(pedn-n))

  call dgemm ('N','N',n,(pedn-n),n,1.0D0,G22,n,A22A21,n,0.0D0,GA22A21,n)
  !print*,H11_tmp 

  nrm(1:k1,1:k1)=A11+H11_tmp2
  !print*,nrm(3,:)
  do i=1,(pedn-n)
    do j=1,n
      nrm(k1+j,i)=GA22A21(j,i)
      nrm(i,k1+j)=GA22A21(j,i)
    end do
  end do
  !print*,nrm(3,:)
  nrm((k1+1):(k1+k2),(k1+1):(k1+k2))=G22

  open (unit=50,file=trim(hrtmx)//'.hrtmx',status='unknown')
  do zi=1,pedn
    do j=1,zi
      if (nrm(idx_total2(zi),idx_total2(j)) .ne. 0) then
        write(50,*)zi,j,nrm(idx_total2(zi),idx_total2(j))
      end if
    end do
  end do
  close(50)

else ! if (mn ==0)

  open (unit=50,file=trim(hrtmx)//'.hrtmx',status='unknown')
  do zi=1,pedn
  !do zi=1,n
    do j=1,zi
      if (nrm(zi,j) .ne. 0) then
        !write(50,*)zi,j,nrm(zi,j)
        write(50,*)zi,j,nrm(zi,j)
      end if
      !if (A22(zi,j) .ne. 0) then
      !  write(50,*)zi,j,A22(zi,j)
      !end if
    end do
  end do
  close(50)

end if !mn > 0



end subroutine
