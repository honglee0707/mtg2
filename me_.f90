subroutine me (ne,length,chrn)

!*********************************************************************
!To get Me given Ne, L and chr # (also see me_sim3.r in /R/julius/)
!S. Hong Lee (2016)
!*********************************************************************

implicit none

integer::ne,chrn
real::length,fxdx,fxdx2

!fxdx=1/((log(length*4*ne+1)+length*4*ne*(log(length*4*ne+1)-1))/(length**2*8*ne**2))
!fxdx2=1/((log(length*2*ne+1)+length*2*ne*(log(length*2*ne+1)-1))/(length**2*4*ne**2))

fxdx=((log(length*4*ne+1)+length*4*ne*(log(length*4*ne+1)-1))/(length**2*8*ne**2))+(1/(3*real(ne)))*(real(chrn)-1)
fxdx=chrn/fxdx

fxdx2=((log(length*2*ne+1)+length*2*ne*(log(length*2*ne+1)-1))/(length**2*4*ne**2))+(1/(3*real(ne)))*(real(chrn)-1)
fxdx2=chrn/fxdx2


print*,'Eq. 10 :',fxdx
print*,'Eq. 11 :',fxdx2




end subroutine
