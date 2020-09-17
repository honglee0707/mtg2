
subroutine mtg_inv (mbin,n,vn,grm_known)

!***********************************************************************
!PCA analysis 
!Author: Sang Hong Lee (c) 2009,2015
!***********************************************************************
 
implicit none
character(len=64)::grm_known,cdum
integer::n,i,j,k,k2,zi,zj,im,ma1,info,vn,thn,M
integer::lwork,liwork,IL,IU

real::mbin(1,n,n)          !assumin single GRM

double precision,allocatable::eivec(:,:),eival(:),Z(:,:)
double precision,allocatable::work(:)
integer,allocatable::iwork(:),ISUPPZ(:)
double precision::x1,VL,VU,v1

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

allocate(eival(n),eivec(n,n),Z(n,n)) 
allocate(ISUPPZ(2*n)) 




PRINT*,'inverting GRM *****************************************************'
print*,''

eivec=mbin(1,:,:)
call cpu_time (t2_cpu)
call cholesky_inv (eivec,n,v1)
call cpu_time (t1_cpu)
print*,'inverse done - time:',real(t1_cpu-t2_cpu)

open (UNIT=39,FILE=trim(grm_known)//'.inv',STATUS='unknown')
do i=1,n
  do j=1,i
    if (eivec(i,j).ne.0) then
      write(39,*)i,j,real(eivec(i,j))  !,GRM(i,j)
    end if
  end do
end do
close(39)


end subroutine

