
subroutine delta_2 (delta)        

!***********************************************************************
!delta_: delta method to get ratio and its variance 
!S. Hong Lee
!***********************************************************************

implicit none

integer::i,j,ti,ttn,tn,tn2,fln(100),fln2(3),idxn,io
double precision,allocatable::vc2(:),vvc2(:,:)
double precision::cov1,vvt1,vt1,vtest1,z1,z2
character(len=7)::cdum
character(len=64)::fl12,delta
character(len=100)::string

!CDF
real ( kind = 8 ) bound,mean,p,q,sd,x
integer ( kind = 4 ) status,which

!print*,delta

ttn=0
open (UNIT=46,FILE=delta,STATUS='old')
do
  ttn=ttn+1
  read(46,*,iostat=io)cdum
  if (io.ne.0) exit
end do
close(46)
ttn=ttn-1
!print*,ttn


open (UNIT=46,FILE=delta,STATUS='old')
read(46,*)fl12
print*,'var-cov-info matrix file: ',trim(fl12)
read(46,*)tn

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

do ti=1,(ttn-2)
  fln=0
  !read(46,*,iostat=io)cdum,fln
  read(46,'(a1,a)')cdum,string
  !print*,cdum,string
  do i=1,100
    read(string,*,iostat=io)(fln(j),j=1,i)
    if (io.ne.0) then
      tn2=i-1
      exit
    end if
  end do
  !print*,cdum,tn2,fln(1:tn2)
  if (cdum == "R" .or. "r") then
    tn2=0
    do i=1,100
      if (fln(i).ne.0) tn2=tn2+1
    end do

    vt1=0;vvt1=0
    do i=1,tn2
      vt1=vt1+vc2(fln(i))
      do j=1,tn2
        vvt1=vvt1+vvc2(fln(i),fln(j))
      end do
    end do

    idxn=fln(1)
    cov1=0
    do i=1,tn2
      cov1=cov1+vvc2(fln(i),idxn)
    end do

    vtest1=((vc2(idxn)/vt1)**2)*(vvc2(idxn,idxn)/vc2(idxn)**2+vvt1/vt1**2-2*cov1/(vc2(idxn)*vt1))

!CDF functions *****************************************************
      p = huge ( p )
      q = huge ( q )
      x = (vc2(idxn)/vt1)**2/vtest1
      sd = 1.0D+00 !0.875D+00    !degree of greedom for cdfchi
which=1  !Calculate P and Q from X, MEAN and SD;
call cdfchi ( which, p, q, x, sd, status, bound )

    print*,'Ratio:',real(vc2(idxn)/vt1),'SE:',real(vtest1**0.5),' p-value:',real(q)

  else if (cdum=="C" .or. "c") then
    !print*,fln(1:3)
    fln2(3)=fln(1)
    fln2(2)=fln(2)
    fln2(1)=fln(3)
    z1=vc2(fln2(3))/sqrt(vc2(fln2(1))*vc2(fln2(2)))
    z2=vvc2(fln2(1),fln2(1))/(4*vc2(fln2(1))**2)
    z2=z2+vvc2(fln2(2),fln2(2))/(4*vc2(fln2(2))**2)
    z2=z2+vvc2(fln2(3),fln2(3))/(vc2(fln2(3))**2)
    z2=z2+2*vvc2(fln2(1),fln2(2))/(4*vc2(fln2(1))*vc2(fln2(2)))
    z2=z2-2*vvc2(fln2(1),fln2(3))/(2*vc2(fln2(1))*vc2(fln2(3)))
    z2=z2-2*vvc2(fln2(3),fln2(2))/(2*vc2(fln2(3))*vc2(fln2(2)))
    z2=z2*(z1**2)

!CDF functions *****************************************************
      p = huge ( p )
      q = huge ( q )
      x = (z1)**2/z2
      sd = 1.0D+00 !0.875D+00    !degree of greedom for cdfchi
which=1  !Calculate P and Q from X, MEAN and SD;
call cdfchi ( which, p, q, x, sd, status, bound )

    print*,'Cor  :',real(z1),'SE:',real(z2**0.5),' p-value:',real(q)
  end if
  
end do

deallocate(vc2,vvc2)
close(46)


end subroutine
