
subroutine mtg_pca (mbin,n,vn,grm_known)

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
double precision::x1,VL,VU

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

allocate(eival(n),eivec(n,n),Z(n,n)) 
allocate(ISUPPZ(2*n)) 




PRINT*,'PCA analysis *****************************************************'
print*,''

eivec=mbin(1,:,:)
!do i=1,n
!  do j=1,i
!    print*,i,j,eivec(i,j)
!  end do
!end do

call cpu_time (t2_cpu)
!print*,"call dsyevr"

IL=1
IU=vn

lwork=-1
liwork=-1
allocate(work(1)) 
allocate(iwork(1)) 

!call dsyevd ('V','L',n,eivec,n,eival,work,lwork,iwork,liwork,info)
call dsyevr ('V','A','L',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)

!print*,work(1),iwork(1)
lwork=int(abs(work(1)))
deallocate(work)
allocate(work(lwork))

liwork=int(abs(iwork(1)))
deallocate(iwork)
allocate(iwork(liwork))
!print*,n,lwork,liwork

!call dsyev ('V','L',n,eivec,n,eival,work,lwork,info)
!call dsyevd ('V','L',n,eivec,n,eival,work,lwork,iwork,liwork,info)
call dsyevr ('V','A','L',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)

print*,"0 means no error >>> ",info
print*,''
call cpu_time (t1_cpu)
print*,'pca done - time:',real(t1_cpu-t2_cpu)

open (UNIT=39,FILE=trim(grm_known)//'.evec',STATUS='unknown')
do i=1,n
  !write(39,'(100f12.5)')(real(eivec(i,n+1-j)),j=1,vn)
  !write(39,'(100f12.5)')(real(Z(i,n+1-j)),j=1,vn)
  !write(39,*)(real(Z(i,n+1-j)),j=1,vn)
  write(39,'(*(G0,1X))')(real(Z(i,n+1-j)),j=1,vn)
end do
close(39)

open (UNIT=40,FILE=trim(grm_known)//'.eval',STATUS='unknown')
do i=1,n
  write(40,*)real(eival(n+1-i))
end do
close(40)


end subroutine

