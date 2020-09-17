
subroutine mtg_matmat (mbin,n,vn,matmat,grm_known)

!***********************************************************************
!cholesky analysis 
!Author: Sang Hong Lee (c) 2009,2015,2018
!***********************************************************************
 
implicit none
character(len=64)::grm_known,cdum
integer::n,i,j,k,k2,zi,zj,im,ma1,info,vn,thn,M,err
integer::lwork,liwork,IL,IU,matmat

real::mbin(2,n,n)          !assuming two GRM
real::mat(n,n),v1

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

if (matmat==1) then
  print*,'symmetric matrices multiplcation'

  call cpu_time (t2_cpu)
  !$OMP PARALLEL DO PRIVATE (i,j,v1) 
  do i=1,n
    do j=1,n
      v1=0
      do k=1,n
        v1=v1+mbin(1,i,k)*mbin(2,k,j)
      end do
      mat(i,j)=v1
    end do
  end do
  !$OMP END PARALLEL DO
  call cpu_time (t1_cpu)

  print*,'symmetric mat x mat product done - time:',real(t1_cpu-t2_cpu)

  open (UNIT=39,FILE=trim(grm_known)//'.matmat1',STATUS='unknown')
  do i=1,n
    do j=1,n
      write(39,*)i,j,mat(i,j)
    end do
  end do
  close(39)


elseif (matmat==2) then
  print*,'triangular matrices multiplcation (i.e. cholesky)'

  call cpu_time (t2_cpu)
  !$OMP PARALLEL DO PRIVATE (i,j) 
  do i=1,n
    do j=1,(i-1)
      mbin(1,j,i)=0
      mbin(2,i,j)=0
    end do
  end do
  !$OMP END PARALLEL DO

  !$OMP PARALLEL DO PRIVATE (i,j,v1) 
  do i=1,n
    do j=1,n
      v1=0
      do k=1,n
        v1=v1+mbin(1,i,k)*mbin(2,k,j)
      end do
      mat(i,j)=v1
    end do
  end do
  !$OMP END PARALLEL DO
  call cpu_time (t1_cpu)

  print*,'triangular mat x mat product done - time:',real(t1_cpu-t2_cpu)

  open (UNIT=39,FILE=trim(grm_known)//'.matmat2',STATUS='unknown')
  do i=1,n
    do j=1,i
      write(39,*)i,j,mat(i,j)+mat(j,i)
    end do
  end do
  close(39)


end if




end subroutine

