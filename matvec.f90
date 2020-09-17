
subroutine mtg_matvec (mbin,n,vn,matvec,grm_known)

!***********************************************************************
!cholesky analysis 
!Author: Sang Hong Lee (c) 2009,2015,2018
!***********************************************************************
 
implicit none
character(len=64)::grm_known,cdum,matvec
integer::n,i,j,k,k2,zi,zj,im,ma1,info,vn,thn,M,err
integer::lwork,liwork,IL,IU

real::mbin(1,n,n)          !assumin single GRM
real::vec(n),vec2(n),v1

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

PRINT*,'mat x vec product **************************************************'
print*,''

  open (UNIT=44,FILE=matvec,STATUS='old')
  do i=1,n
    read(44,*)vec(i)
  end do
  close(44)

call cpu_time (t2_cpu)
!$OMP PARALLEL DO PRIVATE (v1) 
do i=1,n
  v1=0
  do j=1,i
    v1=v1+mbin(1,i,j)*vec(j)
  end do
  vec2(i)=v1
end do
!$OMP END PARALLEL DO
call cpu_time (t1_cpu)

print*,'mat x vec product done - time:',real(t1_cpu-t2_cpu)

open (UNIT=39,FILE=trim(grm_known)//'.matvec',STATUS='unknown')
do i=1,n
  write(39,*)real(vec2(i))  !,GRM(i,j)
end do
close(39)


end subroutine

