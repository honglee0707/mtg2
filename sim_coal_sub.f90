subroutine sim_coal_sub (sim_coal_file)

implicit none

integer::nm,nqtl,pt,pne,nch,ncg
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
!READ(37,*)spf         !simplify recombination
!READ(37,*)ncg
READ(37,*)empsz      !emplify size
!READ(37,*)ncase,ncont      !emplify size
!READ(37,*)ncase2,ncont2      !emplify size
!READ(37,*)se          !residual stadard error
!READ(37,*)thd          !threshold

!allocate(genom(pne*2,nm),yor(pne))

!print*,'mr:',mr
  !ncase=1000
  !ncont=1000
  !ncase2=0
  !ncont2=0
  !thd=0;se=0
  !ncg=nm/2
  !call genedrop(nm,pne,pt,nch,chl,seed2,mr,yor,genom,spf,ncg,empsz,se,ncase,ncont,ncase2,ncont2,thd)
  call genedrop(nm,pne,pt,nch,chl,seed2,mr,spf,empsz,sim_coal_file)


end subroutine
