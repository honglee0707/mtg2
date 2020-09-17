
subroutine mtg_bend (mbin,n,vn,grm_known,opt,thn)

!***********************************************************************
!PCA analysis 
!Author: Sang Hong Lee (c) 2009,2015
!***********************************************************************
 
implicit none
character(len=64)::grm_known,cdum
integer::n,i,j,k,k2,zi,zj,im,ma1,info,vn,thn,M,np,p,OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
integer::lwork,liwork,IL,IU,opt

real::mbin(1,n,n)          !assumin single GRM

double precision,allocatable::eivec(:,:),eival(:),Z(:,:),v2(:,:)
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


if (opt==1) then

PRINT*,'Bending GRM to be PD *********************************************'
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
!call dsyevr ('V','A','L',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)
!call dsyevr ('V','A','L',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)
call dsyevr ('V','A','U',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)

!print*,work(1),iwork(1)
lwork=int(abs(work(1)))
deallocate(work)
allocate(work(lwork))

liwork=int(abs(iwork(1)))
deallocate(iwork)
allocate(iwork(liwork))
!print*,n,lwork,liwork


101 continue

!call dsyev ('V','L',n,eivec,n,eival,work,lwork,info)
!call dsyevd ('V','L',n,eivec,n,eival,work,lwork,iwork,liwork,info)
!call dsyevr ('V','A','L',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)
call dsyevr ('V','A','U',n,eivec,n,VL,VU,IL,IU,-1.0,M,eival,Z,n,ISUPPZ,work,lwork,iwork,liwork,info)


print*,"0 means no error >>> ",info
print*,''

k=0
do i=1,n
  if (eival(i).le.0.0D0) then
    eival(i)=0.001
    k=k+1
  end if
end do
print*,k

if (k==0) goto 110

!do i=1,n
!  do j=1,i
!    v1=0
!    do k=1,n
!      v1=v1+Z(i,k)*eival(k)*Z(j,k)
!    end do
!    eivec(i,j)=v1
!  end do
!end do

allocate(v2(n,thn))
!$OMP PARALLEL DO PRIVATE(i, j)
do p=1,thn
  do i=p,n,thn
    v2(:,p)= Z(i,:) * eival(:)
    do j=1, i
       eivec(j,i) =  dot_product(v2(:,p), Z(j,:)) 
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
deallocate(v2)
!print*,'parallel finish'

print*,'bent and check again >>>'
goto 101

110 continue

!do i=1,n
!  do j=1,i
!    v1=0
!    do k=1,n
!      v1=v1+Z(i,k)*eival(k)*Z(j,k)
!    end do
!    eivec(i,j)=v1
!  end do
!end do

allocate(v2(n,thn))
!$OMP PARALLEL DO PRIVATE(i, j)
do p=1,thn
  do i=p,n,thn
    v2(:,p)= Z(i,:) * eival(:)
    do j=1, i
      eivec(j,i) =  dot_product(v2(:,p), Z(j,:)) 
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
deallocate(v2)





call cpu_time (t1_cpu)
print*,'bending done - time:',real(t1_cpu-t2_cpu)


elseif (opt==2) then

PRINT*,'Adding 0.01 to the diagonal to make PD GRM ************************'
print*,''

do i=1,n
  do j=1,i
    eivec(j,i)=mbin(1,j,i)
    if (i==j) then
      eivec(j,i)=eivec(j,i)+0.01    !add 0.01 to diagonal
    end if
  end do
end do

else
end if  

open (UNIT=39,FILE=trim(grm_known)//'.bend',STATUS='unknown')
do i=1,n
  do j=1,i
    if (eivec(j,i).ne.0) then
      !write(39,*)i,j,real(eivec(i,j))  !,GRM(i,j)
      write(39,*)i,j,real(eivec(j,i))  !,GRM(i,j)
    end if
  end do
end do
close(39)


end subroutine

