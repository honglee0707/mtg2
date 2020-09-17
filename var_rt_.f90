
subroutine var_rt_ (mbin,pedn,mn,fl_var_rt,fl5)        

!***********************************************************************
!aireml_b: multivariate analysis especially for case control 
!S. Hong Lee

!mbin : matrices bin
!obs  : observed ID
!yv   : phenotypes
!rn   : no. phenotypes
!pedn : no. pedigree (ID) should be = rn
!fn   : no. fixed effects
!mn   : no. matrices (random effects)
!***********************************************************************

implicit none

INTEGER::n,nped,ix,iz,tn,tn2,zi,i,j,nit
integer::rn,pedn,mn,trn,trnx,trny,tfn,trnv,trnu,itit
real::mbin(mn,pedn,pedn)  !,yv(rn,tn),tyv(trn,1)
real::x1,x2,x3
real::x11,x12,x13,x14,x15,x16,x17,x18

character(len=7)::cdum
character(len=64)::fl4,fl5,fl11,fl_var_rt,fl12,fl13
logical::file_exist
integer,allocatable:: disc(:), valid(:)

nit=1
open (UNIT=46,FILE=fl_var_rt,STATUS='old')
!read(46,*)nit        !# replicates
read(46,*)tn,tn2     !# discovery, # validation
read(46,*)fl12,fl13  !file for discovery list, file for validation list
close(46)

print*,'discovery sample list: ',trim(fl12)
print*,'target sample list   : ',trim(fl13)
print*,'output file          : ',trim(fl5)
print*,''

allocate(disc(tn),valid(tn2))

open (UNIT=49,FILE=fl5,STATUS='unknown')
open (UNIT=47,FILE=fl12,STATUS='old')
open (UNIT=48,FILE=fl13,STATUS='old')

do zi=1,nit
  read(47,*)disc(:)
  read(48,*)valid(:)

  x1=0;x2=0;x3=0
  do i=1,tn
    do j=1,tn2
      x3=x3+1
      x1=x1+mbin(1,disc(i),valid(j))
      x2=x2+mbin(1,disc(i),valid(j))**2
    end do
  end do
  !write(49,*)x1/x3,(x2-(x1**2/x3))/(x3-1)
  print*,'Overall mean     : ',x1/x3 !,(x2-(x1**2/x3))/(x3-1)
  print*,'Overall variance : ',(x2-(x1**2/x3))/(x3-1)

end do
!write(49,*)''

if (nit==1) then
  do j=1,tn2
    x11=0;x12=0;x13=0;x14=0;x15=0;x16=0;x17=0;x18=0
    do i=1,tn
      x13=x13+1
      x11=x11+mbin(1,disc(i),valid(j))
      x12=x12+mbin(1,disc(i),valid(j))**2
      x14=x14+(mbin(1,disc(i),valid(j))-x1/x3)**2
      x15=x15+(mbin(1,disc(i),valid(j))-0)**2
      if (mbin(1,disc(i),valid(j))>0) then
        x18=x18+1
        x16=x16+mbin(1,disc(i),valid(j))
        x17=x17+mbin(1,disc(i),valid(j))**2
      end if
    end do
    !write(49,*)j,valid(j),x11/x13,(x12-(x11**2/x13))/(x13-1),x1/x3,x14/(x13-1)
    !write(49,'(*(G0,1X))')j,valid(j),x11/x13,(x12-(x11**2/x13))/(x13-1),x14/(x13-1),x15/(x13-1), &
&   !   x16/x13,(x17-(x16**2/x13))/(x13-1),x16/x18,(x17-(x16**2/x18))/(x18-1)
    write(49,'(*(G0,1X))')j,valid(j),x11/x13,(x12-(x11**2/x13))/(x13-1)

  end do
end if

close(48)
close(47)
close(49)



end subroutine

