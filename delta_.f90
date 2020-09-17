
subroutine delta_ (delta)        

!***********************************************************************
!delta_: delta method to get ratio and its variance 
!S. Hong Lee
!***********************************************************************

implicit none

integer::i,j,tn,tn2,fln(100),idxn
double precision,allocatable::vc2(:),vvc2(:,:)
double precision::cov1,vvt1,vt1,vtest1,z1,z2
character(len=7)::cdum
character(len=64)::fl12,delta

!print*,delta

open (UNIT=46,FILE=delta,STATUS='old')
read(46,*)fl12
read(46,*)tn,tn2
read(46,*)fln(1:tn2)
read(46,*)idxn
close(46)

print*,fl12

allocate(vc2(tn),vvc2(tn,tn))

open (UNIT=47,FILE=fl12,STATUS='old')
do i=1,tn
  read(47,*)cdum,vc2(i)
  !print*,vc2(i)
end do
do i=1,tn
  read(47,*)(vvc2(i,j),j=1,i)
  do j=1,i
    vvc2(j,i)=vvc2(i,j)
  end do
end do
close(47)

!do i=1,tn
!  print*,vvc2(i,:)
!end do

vt1=0;vvt1=0
do i=1,tn2
  vt1=vt1+vc2(fln(i))
  do j=1,tn2
    vvt1=vvt1+vvc2(fln(i),fln(j))
  end do
end do

cov1=0
do i=1,tn2
  cov1=cov1+vvc2(fln(i),idxn)
end do

vtest1=((vc2(idxn)/vt1)**2)*(vvc2(idxn,idxn)/vc2(idxn)**2+vvt1/vt1**2-2*cov1/(vc2(idxn)*vt1))

print*,'ratio: ',vc2(idxn)/vt1
print*,'SE   : ',vtest1**0.5

if (tn2==3) then !correlation and SE
  z1=vc2(fln(3))/sqrt(vc2(fln(1))*vc2(fln(2)))
  z2=vvc2(fln(1),fln(1))/(4*vc2(fln(1))**2)
  z2=z2+vvc2(fln(2),fln(2))/(4*vc2(fln(2))**2)
  z2=z2+vvc2(fln(3),fln(3))/(vc2(fln(3))**2)
  z2=z2+2*vvc2(fln(1),fln(2))/(4*vc2(fln(1))*vc2(fln(2)))
  z2=z2-2*vvc2(fln(1),fln(3))/(2*vc2(fln(1))*vc2(fln(3)))
  z2=z2-2*vvc2(fln(3),fln(2))/(2*vc2(fln(3))*vc2(fln(2)))
  z2=z2*(z1**2)

print*,''
print*,'cor  : ',z1
print*,'SE   : ',z2**0.5
  
end if


end subroutine
