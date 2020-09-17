subroutine sim_coal_sub2 (sim_coal_file)

implicit none

integer::nm,nqtl,pt,pne,nch,ncg,pt2,pne2,dec12
integer::i,j,k,zi,k2,seed2,spf,empsz,ncase,ncont,ncase2,ncont2
double precision::mr,chl
!real::thd,se,mr,chl
real,allocatable::yor(:)
character(len=1),allocatable::genom(:,:)
character(len=64)::sim_coal_file

!open (UNIT=37,FILE='simdis.in',STATUS='old')
open (UNIT=37,FILE=sim_coal_file,STATUS='old')

spf=0
READ(37,*)nm,nch,chl  !,nqtl    !some SNP will be selected as causal mutation
READ(37,*)pt,pne      ! generation, Ne
READ(37,*)mr          !mutation rate
READ(37,*)seed2       !random seed
READ(37,*)empsz      !emplify size
READ(37,*)pt2,pne2,dec12      !emplify size

!allocate(genom(pne*2,nm),yor(pne))

!print*,'mr:',mr
  !ncase=1000
  !ncont=1000
  !ncase2=0
  !ncont2=0
  !thd=0;se=0
  !ncg=nm/2
  !call genedrop(nm,pne,pt,nch,chl,seed2,mr,yor,genom,spf,ncg,empsz,se,ncase,ncont,ncase2,ncont2,thd)
  call genedrop2(nm,pne,pt,nch,chl,seed2,mr,spf,empsz,sim_coal_file,pne2,pt2,dec12)


end subroutine
