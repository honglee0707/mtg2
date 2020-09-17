
subroutine t2b_convert (mbin,n,grm_known,t2b)

!***********************************************************************
!converting GRM format 
!Author: Sang Hong Lee (c) 2019
!***********************************************************************

implicit none
character(len=64)::grm_known
integer::n,i,j,k,t2b

real::mbin(1,n,n)          !assumin single GRM


if (t2b == 1) then
  print*,"converting txt to bin format"
  !open (UNIT=39,FILE=trim(grm_known)//'.bin',STATUS='unknown',form='unformatted',access='direct',recl=4,position="append")
  open (UNIT=39,FILE=trim(grm_known)//'.bin',STATUS='unknown',form='unformatted',access='direct',recl=1)
  k=0
  do i=1,n
    do j=1,i
      k=k+1
      write(39,rec=k)mbin(1,i,j)
    end do
  end do

elseif (t2b==2) then
  print*,"converting bin to txt format"
  open (UNIT=39,FILE=trim(grm_known)//'.txt',STATUS='unknown')
    do i=1,n
      do j=1,i
        write(39,*)i,j,mbin(1,i,j)
      end do
    end do
end if


end subroutine


