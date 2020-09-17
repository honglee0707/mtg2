subroutine genedrop2(nmarkers,neff,ngen,nch,chl,seed2,mr,spf,empsz,sim_coal_file,neff2,ngen2,dec12)

!************************************************************************
! Coalescence genedroping simulaiton 
! initiated from Julius' genedrop spread sheet 
! S. Hong Lee (2005, 2016) 
!************************************************************************ 

implicit none

integer::rseed,nmarkers,neff,ngen,neff1,neff2,ngen2,seed2,spf,prn,empsz,dec12
parameter(prn=15)
character(len=64)::sim_coal_file

double precision::chl,mr

real::rnd,dmq,dmq2,recq,recq2,dm,thd
integer::i,j,m1,km,k,k2,k1,pedn,an
integer::zi,w,xr,xr2,idd,ids,idx,im,nsel,ipar1,igen,nk2,nk1
integer::ncase,ncont,ncase2,ncont2 
integer::midxn,midx(nmarkers)
integer::nmal2,nfem2,kkk,popmal,popfem,nfem,popsiz,nmal
double precision::dseed

character(len=1),allocatable::genom(:,:),genom1(:,:)
character(len=1)::ck
character(len=1)::mafreq1(nmarkers,2)
integer,allocatable::nselmal(:),nselfem(:),seed(:),ibd(:,:)

double precision::binpr(0:prn),mutpr(0:prn),pr,v1,v2,v3,v10
real::gbv(neff),rsd(neff),gbv2(neff*empsz),rsd2(neff*empsz)
real::mafreq2(nmarkers,2),x1,x2,x10,x11,mux,sdx,gasdev,mux2,sdx2

integer::trl,iv1,recn,recid(100),mutn,mutid(100),nch,nmch,disy(neff*empsz)
integer::ped(neff*empsz,2)

integer::contped(neff*empsz,3),ncont10
character(len=1)::contom(neff*empsz*2,nmarkers)

integer::mind(nmarkers),ich,ma1,mano(nmarkers)

integer,allocatable::srtv(:,:),matpn(:,:),matcn(:,:),recmatp(:,:,:),recmatc(:,:,:)
integer,allocatable::mutpn(:,:),mutcn(:,:),mutmatp(:,:,:),mutmatc(:,:,:)

!print*,nmarkers,neff,ngen,nch,chl,seed2,mr,spf,empsz,sim_coal_file

! random number generator
call RANDOM_SEED(SIZE=rseed)
ALLOCATE(seed(rseed))
seed=seed2
call RANDOM_SEED(put=seed)
dseed=seed2

!nch=20     !no. chromosomes
nmch=nmarkers/nch      !no. QTL per each chr.


allocate(nselmal(neff*2),nselfem(neff*2))
allocate(srtv(0,0))
allocate(matpn(nch,neff*2),matcn(nch,neff*2),recmatp(nch,neff*4,1000),recmatc(nch,neff*4,1000))
allocate(mutpn(nch,neff*2),mutcn(nch,neff*2),mutmatp(nch,neff*2,1000),mutmatc(nch,neff*2,1000))
allocate(genom1(neff*2,nmarkers),genom(neff*2,nmarkers))


!print*,'no. mutations per genom:',int(nmarkers*mr)
!print*,''


!binomial probability work out - binpr(), mutpr() **************************
! (trl! / i! (trl-i)!) pr^i (1-pr)^(trl-i)  - binomial prob.

trl=nmch         !trial - nmarkers/nch
!pr=chl*.1/(.1*nmch)  !prob for occuring incidence (assuming 1M long - approx.)
pr=chl/real(nmch) !prob for occuring incidence (assuming 1M long - approx.)
!print*,pr,chl,nmch
!pause

v1=0
do i=1,trl
  v10=i
  v1=v1+log(v10)  ! trl!
end do
!print*,v1

do i=0,prn                  !incidence 0 ~ prn
  v2=0
  do j=1,i
    v10=j
    v2=v2+log(v10) ! i!
  end do
  !print*,v2

  v3=0
  do j=1,trl-i
    v10=j
    v3=v3+log(v10)  ! (trl-i)!
  end do
  !print*,v3

  binpr(i)=v1-(v2+v3)
  mutpr(i)=binpr(i)+i*log(mr)+(trl-i)*log(1-mr)
  binpr(i)=binpr(i)+i*log(pr)+(trl-i)*log(1-pr)
  !print*,binpr(i)
  !pause
end do
binpr=exp(binpr)
mutpr=exp(mutpr)
!print*,binpr(:)
!print*,mutpr(:)
!print*,mr,pr
!pause

do i=1,prn
  binpr(i)=binpr(i)+binpr(i-1)     !cumulative
  mutpr(i)=mutpr(i)+mutpr(i-1)     !cumulative
end do


!binpr=1
!mutpr=1

!print*,'**** cumulative probability for number of recombination per chr. ****'
!print*,binpr(:)
!print*,''
!print*,'**** cumulative probability for number of mutaitons per chr. ****'
!print*,mutpr(:)
!print*,''
!print*,'NOTE: cumulative prob. should be 1 or very close to 1'
!print*,''
!pause

print*,'**** Stochastic coalescence simulation *****************'
print*,''
print*,'#generation #recombination #mutation'

  nmal=neff/2
  nfem=neff/2  

  !SNP simulation ******************************************************
  !print*,'base generation *************************'
  do i=1,neff
    do j=1,nmarkers
      ck='1'
      call random_number(rnd)
      if (rnd.lt.0.5) THEN
        ck='2'
      END if
      genom(i*2-1,j)=ck  
      !genom(i*2-1,j)='1'     !total homozygous  

      ck='1'
      call random_number(rnd)
      if (rnd.lt.0.5) THEN
        ck='2'
      END if
      genom(i*2,j)=ck  
      !genom(i*2,j)='1'        !total homozygous  
    end do
  end do !**************************************************************
  !print*,'base generation assigned *****************'
  !print*,''

  !recombination matrix
  do ich=1,nch
    matpn(ich,:)=1
  end do
  do i=1,neff
    do ich=1,nch
      recmatp(ich,i*4-4+1,1)=i
      recmatp(ich,i*4-4+2,1)=1
      recmatp(ich,i*4-4+3,1)=i
      recmatp(ich,i*4-4+4,1)=1
      recmatp(ich,i*4-4+1,2)=1
      recmatp(ich,i*4-4+2,2)=nmch
      recmatp(ich,i*4-4+3,2)=2
      recmatp(ich,i*4-4+4,2)=nmch
    end do
  end do
  !print*,recmatp(1,1,1:10)
  !pause

  !mutation matrix
  do ich=1,nch
    mutpn(ich,:)=0
  end do

  !create next generation
  popsiz=neff
  popmal=nmal
  popfem=nfem

  neff1=neff

  kkk=0
  do igen=1,ngen   !historical generation starts********************
    if (igen > ngen2) then 
      if (dec12==1) then          !linealry decreasing
        neff=int(neff-real(neff1-neff2)/real(ngen-ngen2))
      elseif (dec12==2) then
        neff=neff2
      else
        print*,'check >>> line 202'
      end if

      nmal=neff/2
      nfem=neff/2  

      popsiz=neff
      popmal=nmal
      popfem=nfem
    end if


     print'(*(G0,10X))','',igen,matcn(1,1),mutcn(1,1),neff
    if (matcn(1,1)>400) then
      print*,'check >>> genomic length per chr > 5 '  !check
      print*,'matcn should re-allocate with > 1000'  !check
      pause
    end if
    nk1=0;nk2=0
    do i=1,popsiz*2
      call random_number(rnd)
      nsel=int(rnd*nmal)+1     !random selection and mating 
      nselmal(i)=nsel
    end do
    do i=1,popsiz*2
      call random_number(rnd)
      !print*,i,rnd,popsiz 
      nsel=nmal+int(rnd*nfem)+1
      !print*,nsel,nmal2,nfem2
      !pause
      nselfem(i)=nsel
    end do

    !create progeny
    do ipar1=1,popsiz*2
      !print*,ipar1
      ids=nselmal(ipar1)
      idd=nselfem(ipar1)
      if (nk1.ne.popmal.or.nk2.ne.popfem) then
220     call random_number(rnd)
        if (rnd.lt.0.5.and.nk1.lt.popmal) then
          nk1=nk1+1
          idx=nk1
        ELSEIF (rnd.ge.0.5.and.nk2.lt.popfem) then
          nk2=nk2+1
          idx=popmal+nk2
        else
          GOTO 220
        END if

        !print*,idx,ids,idd
        !pause

        !From dad**********************************************************
        do ich=1,nch
          !print*,ich
          call random_number(rnd)
          xr=INT(2*rnd)+1         !random determination
          !print*,xr
          xr2=3-xr       !for mutation
          !print*,xr
          !pause

          !determine how manay times rec. occur - recn *********************
          call random_number (rnd)
          do recn=0,prn    !10
            if (rnd < binpr(recn)) exit
          end do
          !print*,recn

          !simplified - one recombination per chr **************************
          if (spf==1) recn=1
          !*****************************************************************

          do i=1,recn
145         call random_number (rnd)
            !print*,rnd
            iv1=int((nmch-1)*rnd)+2   !avoiding 1 
            do j=1,i-1
              if (recid(j)==iv1) goto 145    !recombined ID
            end do
            recid(i)=iv1
          end do  
          !sorting************************'
          deallocate(srtv)
          allocate(srtv(recn,2))
          do i=1,recn
            srtv(i,2)=recid(i)
          end do
          call sort (srtv,recn)
          do i=1,recn
            recid(i)=srtv(i,2)
          end do !**********************************************************
          !print*,recid
          !pause

          !print*,recid(1:recn),'xr=',xr
          !print*,'ids:',matpn(ich,ids*2-1),matpn(ich,ids*2)
          !print*,recmatp(ich,ids*4-4+1,1:50)
          !print*,recmatp(ich,ids*4-4+2,1:50)
          !print*,recmatp(ich,ids*4-4+3,1:50)
          !print*,recmatp(ich,ids*4-4+4,1:50)
          !print*,''


          !j=1     !indication for recombined marker ID
          !do im=1,nmch       ! markers + QTL no.
          !  if (recid(j)==im) then
          !    xr=3-xr
          !    j=j+1
          !  end if
          !  genom1(idx*2-1,(ich-1)*nmch+im)=genom(ids*2-(2-xr),(ich-1)*nmch+im)
          !end do

          if (recn==0) then
            matcn(ich,idx*2-1)=matpn(ich,ids*2-2+xr)
            recmatc(ich,idx*4-4+1,1:matcn(ich,idx*2-1)*2)=recmatp(ich,ids*4-4+xr*2-1,1:matcn(ich,idx*2-1)*2)
            recmatc(ich,idx*4-4+2,1:matcn(ich,idx*2-1)*2)=recmatp(ich,ids*4-4+xr*2,1:matcn(ich,idx*2-1)*2)

            !mutation inheritance ***********************************
            mutcn(ich,idx*2-1)=mutpn(ich,ids*2-2+xr)
            mutmatc(ich,idx*2-1,1:mutcn(ich,idx*2-1))=mutmatp(ich,ids*2-2+xr,1:mutcn(ich,idx*2-1))

            goto 1001
          end if

          !mutation inheritance *************************************
          iv1=1; mutcn(ich,idx*2-1)=0
          do i=1,recn
            xr2=3-xr2
            do j=1,mutpn(ich,ids*2-2+xr2)
              if (mutmatp(ich,ids*2-2+xr2,j).ge.iv1.and.   &
                &  mutmatp(ich,ids*2-2+xr2,j)<recid(i)) then
                mutcn(ich,idx*2-1)=mutcn(ich,idx*2-1)+1
                mutmatc(ich,idx*2-1,mutcn(ich,idx*2-1))=mutmatp(ich,ids*2-2+xr2,j)
              end if
            end do
            iv1=recid(i)
          end do !***************************************************
          !finalizing*******************************************
          xr2=3-xr2
          do j=1,mutpn(ich,ids*2-2+xr2)
            if (mutmatp(ich,ids*2-2+xr2,j).ge.iv1.and.   &
              &  mutmatp(ich,ids*2-2+xr2,j).le.nmch) then
              mutcn(ich,idx*2-1)=mutcn(ich,idx*2-1)+1
              mutmatc(ich,idx*2-1,mutcn(ich,idx*2-1))=mutmatp(ich,ids*2-2+xr2,j)
            end if
          end do !*****************************************************

          k2=0
          j=1    !1 ~ first rec
          do k1=1,matpn(ich,ids*2-2+xr)     !prev recmat no.
            if (recid(j)-1.ge.recmatp(ich,ids*4-4+xr*2,k1*2-1) .and.  &
              &  recid(j)-1.le.recmatp(ich,ids*4-4+xr*2,k1*2)) then
              k2=k2+1
              recmatc(ich,idx*4-4+2,k2*2-1)=recmatp(ich,ids*4-4+xr*2,k1*2-1)
              recmatc(ich,idx*4-4+2,k2*2)=recid(j)-1
              recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1)
              recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)
              exit   !done for the first rec
            elseif (recid(j)-1.gt.recmatp(ich,ids*4-4+xr*2,k1*2)) then
              k2=k2+1
              recmatc(ich,idx*4-4+2,k2*2-1)=recmatp(ich,ids*4-4+xr*2,k1*2-1) !copy
              recmatc(ich,idx*4-4+2,k2*2)=recmatp(ich,ids*4-4+xr*2,k1*2)    !copy
              recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1) !copy
              recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)    !copy
            elseif (recid(j).lt.recmatp(ich,ids*4-4+xr*2,k1*2-1)) then
              print*,'how can it be >>> check'
              print*,recid(j),'xr=',xr,recid(1:recn)
              print*,'ids:',matpn(ich,ids*2-1),matpn(ich,ids*2)
              print*,recmatp(ich,ids*4-4+1,1:50)
              print*,recmatp(ich,ids*4-4+2,1:50)
              print*,recmatp(ich,ids*4-4+3,1:50)
              print*,recmatp(ich,ids*4-4+4,1:50)
              pause
            end if
          end do

          do j=2,recn
            xr=3-xr    !rec occur
            do k1=1,matpn(ich,ids*2-2+xr)     !finalyging prev recmat no.
              if (recid(j-1).ge.recmatp(ich,ids*4-4+xr*2,k1*2-1) .and.  &
                &  recid(j-1).le.recmatp(ich,ids*4-4+xr*2,k1*2)) then
                k2=k2+1
                recmatc(ich,idx*4-4+2,k2*2-1)=recid(j-1)
                recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1)
                if (recid(j)-1.gt.recmatp(ich,ids*4-4+xr*2,k1*2)) then
                  recmatc(ich,idx*4-4+2,k2*2)=recmatp(ich,ids*4-4+xr*2,k1*2)
                  recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)
                else
                  recmatc(ich,idx*4-4+2,k2*2)=recid(j)-1
                  recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)
                  goto 211   !done for the rec
                end if
              end if    !finalizing finish

              if (recid(j)-1.ge.recmatp(ich,ids*4-4+xr*2,k1*2-1) .and.  &
                &  recid(j)-1.le.recmatp(ich,ids*4-4+xr*2,k1*2)) then
                k2=k2+1
                recmatc(ich,idx*4-4+2,k2*2-1)=recmatp(ich,ids*4-4+xr*2,k1*2-1)
                recmatc(ich,idx*4-4+2,k2*2)=recid(j)-1
                recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1)
                recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)
                goto 211   !done for the rec
              elseif (recid(j)-1.gt.recmatp(ich,ids*4-4+xr*2,k1*2).and. &
                & recmatc(ich,idx*4-4+2,k2*2).lt.recmatp(ich,ids*4-4+xr*2,k1*2-1)) then
                k2=k2+1
                recmatc(ich,idx*4-4+2,k2*2-1)=recmatp(ich,ids*4-4+xr*2,k1*2-1) !copy
                recmatc(ich,idx*4-4+2,k2*2)=recmatp(ich,ids*4-4+xr*2,k1*2)    !copy
                recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1) !copy
                recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)    !copy
              elseif (recid(j).lt.recmatp(ich,ids*4-4+xr*2,k1*2-1)) then
                print*,'how can it be >>> check - can be for this time?'
                !nothing to do
                pause
              end if
            end do
211         continue

          end do

          j=recn    !the last rec ~ last marker
          xr=3-xr    !rec occur
          do k1=1,matpn(ich,ids*2-2+xr)     !prev recmat no.
            if (recid(j).ge.recmatp(ich,ids*4-4+xr*2,k1*2-1) .and.  &
              &  recid(j).le.recmatp(ich,ids*4-4+xr*2,k1*2)) then
              k2=k2+1
              recmatc(ich,idx*4-4+2,k2*2-1)=recid(j)
              recmatc(ich,idx*4-4+2,k2*2)=recmatp(ich,ids*4-4+xr*2,k1*2)
              recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1)
              recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)
              if (recmatc(ich,idx*4-4+2,k2*2)==nmch) then
                exit       !finish
              end if
            end if
            if (recmatc(ich,idx*4-4+2,k2*2)+1==recmatp(ich,ids*4-4+xr*2,k1*2-1)) then
                k2=k2+1
                recmatc(ich,idx*4-4+2,k2*2-1)=recmatp(ich,ids*4-4+xr*2,k1*2-1) !copy
                recmatc(ich,idx*4-4+2,k2*2)=recmatp(ich,ids*4-4+xr*2,k1*2)    !copy
                recmatc(ich,idx*4-4+1,k2*2-1)=recmatp(ich,ids*4-4+xr*2-1,k1*2-1) !copy
                recmatc(ich,idx*4-4+1,k2*2)=recmatp(ich,ids*4-4+xr*2-1,k1*2)    !copy
              if (recmatc(ich,idx*4-4+2,k2*2)==nmch) then
                exit       !finish
              end if
            end if
          end do

          matcn(ich,idx*2-1)=k2
1001      continue

          !check errors ***************************************************
          do i=1,matcn(ich,idx*2-1)-1
            if (recmatc(ich,idx*4-4+2,i*2)+1.ne.recmatc(ich,idx*4-4+2,(i+1)*2-1)) then
              print*,'xr=',xr,'recid:',recid(1:recn)
              print*,'ids:',matpn(ich,ids*2-1),matpn(ich,ids*2)
              print*,recmatp(ich,ids*4-4+1,1:50)
              print*,recmatp(ich,ids*4-4+2,1:50)
              print*,recmatp(ich,ids*4-4+3,1:50)
              print*,recmatp(ich,ids*4-4+4,1:50)
              print*,''
              print*,'matn=',matcn(ich,idx*2-1)
              print*,recmatc(ich,idx*4-4+1,1:50)
              print*,recmatc(ich,idx*4-4+2,1:50)
              pause
            end if
          end do !*********************************************************
          if (recmatc(ich,idx*4-4+2,matcn(ich,idx*2-1)*2).ne.nmch) then
              print*,'xr=',xr,'recid:',recid(1:recn)
              print*,'ids:',matpn(ich,ids*2-1),matpn(ich,ids*2)
              print*,recmatp(ich,ids*4-4+1,1:50)
              print*,recmatp(ich,ids*4-4+2,1:50)
              print*,recmatp(ich,ids*4-4+3,1:50)
              print*,recmatp(ich,ids*4-4+4,1:50)
              print*,''
              print*,'matn=',matcn(ich,idx*2-1)
              print*,recmatc(ich,idx*4-4+1,1:50)
              print*,recmatc(ich,idx*4-4+2,1:50)
              pause
          end if


          !determine how manay times mutation occur - mutn *****************
          call random_number (rnd)
          do mutn=0,prn
            if (rnd < mutpr(mutn)) exit
          end do

          do i=1,mutn       !no. mutations per chr.
143         call random_number(rnd)
            iv1=int(nmch*rnd)+1     !nmch - no. markers in one chr.
            do j=1,i-1
              if (mutid(j)==iv1) goto 143    !mutated ID
            end do
            mutid(i)=iv1
            mutcn(ich,idx*2-1)=mutcn(ich,idx*2-1)+1
            mutmatc(ich,idx*2-1,mutcn(ich,idx*2-1))=iv1
          end do

          !if (igen==2) then
          !  print*,'xr=',xr,'recid:',recid(1:recn)
          !  print*,'mut:',mutid(1:mutn)
          !  print*,'ids mut no:',mutpn(ich,ids*2-1),mutpn(ich,ids*2)
          !  print*,mutmatp(ich,ids*2-2+1,1:50)
          !  print*,mutmatp(ich,ids*2-2+2,1:50)
          !  print*,'idx mut no:',mutcn(ich,idx*2-1)
          !  print*,mutmatc(ich,idx*2-1,1:50)
          !  pause
          !end if

        end do  !ich


        !From mam**********************************************************'
        do ich=1,nch
          !print*,'mum',ich
          call random_number(rnd)
          xr=INT(2*rnd)+1         !random determination
          xr2=3-xr            !for mutaitons

          !determine how manay times rec. occur - recn *********************
          call random_number (rnd)
          do recn=0,prn     !10
            if (rnd < binpr(recn)) exit
          end do
          !print*,recn

          !simplified - one recombination per chr **************************
          if (spf==1) recn=1
          !*****************************************************************

          do i=1,recn
149         call random_number (rnd)
            iv1=int((nmch-1)*rnd)+2   !avoiding 1 
            do j=1,i-1
              if (recid(j)==iv1) goto 149    !recombined ID
            end do
            recid(i)=iv1
          end do  
          !sorting************************
          deallocate(srtv)
          allocate(srtv(recn,2))
          do i=1,recn
            srtv(i,2)=recid(i)
          end do
          call sort (srtv,recn)
          do i=1,recn
            recid(i)=srtv(i,2)
          end do !**********************************************************

          !j=1     !indication for recombined marker ID
          !do im=1,nmch       ! markers + QTL no.
          !  if (recid(j)==im) then
          !    xr=3-xr
          !    j=j+1
          !  end if
          !  genom1(idx*2,(ich-1)*nmch+im)=genom(idd*2-(2-xr),(ich-1)*nmch+im)
          !end do

          !print*,recid(1:recn),'xr=',xr
          !print*,'idd:',matpn(ich,idd*2-1),matpn(ich,idd*2)
          !print*,recmatp(ich,idd*4-4+1,1:50)
          !print*,recmatp(ich,idd*4-4+2,1:50)
          !print*,recmatp(ich,idd*4-4+3,1:50)
          !print*,recmatp(ich,idd*4-4+4,1:50)
          !pause


          if (recn==0) then
            matcn(ich,idx*2)=matpn(ich,idd*2-2+xr)
            recmatc(ich,idx*4-4+3,1:matcn(ich,idx*2)*2)=recmatp(ich,idd*4-4+xr*2-1,1:matcn(ich,idx*2)*2)
            recmatc(ich,idx*4-4+4,1:matcn(ich,idx*2)*2)=recmatp(ich,idd*4-4+xr*2,1:matcn(ich,idx*2)*2)

            !mutation inheritance ***********************************
            mutcn(ich,idx*2)=mutpn(ich,idd*2-2+xr)
            mutmatc(ich,idx*2,1:mutcn(ich,idx*2))=mutmatp(ich,idd*2-2+xr,1:mutcn(ich,idx*2))

            goto 1002
          end if

          !mutation inheritance *************************************
          iv1=1; mutcn(ich,idx*2)=0
          do i=1,recn
            xr2=3-xr2
            do j=1,mutpn(ich,idd*2-2+xr2)
              if (mutmatp(ich,idd*2-2+xr2,j).ge.iv1.and.   &
                &  mutmatp(ich,idd*2-2+xr2,j)<recid(i)) then
                mutcn(ich,idx*2)=mutcn(ich,idx*2)+1
                mutmatc(ich,idx*2,mutcn(ich,idx*2))=mutmatp(ich,idd*2-2+xr2,j)
              end if
            end do
            iv1=recid(i)
          end do !***************************************************
          !finalizing*******************************************
          xr2=3-xr2
          do j=1,mutpn(ich,idd*2-2+xr2)
            if (mutmatp(ich,idd*2-2+xr2,j).ge.iv1.and.   &
              &  mutmatp(ich,idd*2-2+xr2,j).le.nmch) then
              mutcn(ich,idx*2)=mutcn(ich,idx*2)+1
              mutmatc(ich,idx*2,mutcn(ich,idx*2))=mutmatp(ich,idd*2-2+xr2,j)
            end if
          end do !*****************************************************

          k2=0
          j=1    !1 ~ first rec
          do k1=1,matpn(ich,idd*2-2+xr)     !prev recmat no.
            if (recid(j)-1.ge.recmatp(ich,idd*4-4+xr*2,k1*2-1) .and.  &
              &  recid(j)-1.le.recmatp(ich,idd*4-4+xr*2,k1*2)) then
              k2=k2+1
              recmatc(ich,idx*4-4+4,k2*2-1)=recmatp(ich,idd*4-4+xr*2,k1*2-1)
              recmatc(ich,idx*4-4+4,k2*2)=recid(j)-1
              recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1)
              recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)
              exit   !done for the first rec
            elseif (recid(j)-1.gt.recmatp(ich,idd*4-4+xr*2,k1*2)) then
              k2=k2+1
              recmatc(ich,idx*4-4+4,k2*2-1)=recmatp(ich,idd*4-4+xr*2,k1*2-1) !copy
              recmatc(ich,idx*4-4+4,k2*2)=recmatp(ich,idd*4-4+xr*2,k1*2)    !copy
              recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1) !copy
              recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)    !copy
            elseif (recid(j).lt.recmatp(ich,idd*4-4+xr*2,k1*2-1)) then
              print*,'how can it be >>> check'
              pause
            end if
          end do

          do j=2,recn
            xr=3-xr    !rec occur
            do k1=1,matpn(ich,idd*2-2+xr)     !prev recmat no.
              if (recid(j-1).ge.recmatp(ich,idd*4-4+xr*2,k1*2-1) .and.  &
                &  recid(j-1).le.recmatp(ich,idd*4-4+xr*2,k1*2)) then
                k2=k2+1
                recmatc(ich,idx*4-4+4,k2*2-1)=recid(j-1)
                recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1)
                if (recid(j)-1.gt.recmatp(ich,idd*4-4+xr*2,k1*2)) then
                  recmatc(ich,idx*4-4+4,k2*2)=recmatp(ich,idd*4-4+xr*2,k1*2)
                  recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)
                else
                  recmatc(ich,idx*4-4+4,k2*2)=recid(j)-1
                  recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)
                  goto 221   !done for the rec
                end if
              end if    !finalizing

              if (recid(j)-1.ge.recmatp(ich,idd*4-4+xr*2,k1*2-1) .and.  &
                &  recid(j)-1.le.recmatp(ich,idd*4-4+xr*2,k1*2)) then
                k2=k2+1
                recmatc(ich,idx*4-4+4,k2*2-1)=recmatp(ich,idd*4-4+xr*2,k1*2-1)
                recmatc(ich,idx*4-4+4,k2*2)=recid(j)-1
                recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1)
                recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)
                goto 221   !done for the rec
              elseif (recid(j)-1.gt.recmatp(ich,idd*4-4+xr*2,k1*2).and. &
                & recmatc(ich,idx*4-4+4,k2*2).lt.recmatp(ich,idd*4-4+xr*2,k1*2-1)) then
                k2=k2+1
                recmatc(ich,idx*4-4+4,k2*2-1)=recmatp(ich,idd*4-4+xr*2,k1*2-1) !copy
                recmatc(ich,idx*4-4+4,k2*2)=recmatp(ich,idd*4-4+xr*2,k1*2)    !copy
                recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1) !copy
                recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)    !copy
              elseif (recid(j).lt.recmatp(ich,idd*4-4+xr*2,k1*2-1)) then
                print*,'how can it be >>> check - can be for this time'
                !nothing to do
                pause
              end if
            end do
221         continue

          end do

          j=recn    !the last rec ~ last marker
          xr=3-xr    !rec occur
          do k1=1,matpn(ich,idd*2-2+xr)     !prev recmat no.
            if (recid(j).ge.recmatp(ich,idd*4-4+xr*2,k1*2-1) .and.  &
              &  recid(j).le.recmatp(ich,idd*4-4+xr*2,k1*2)) then
              k2=k2+1
              recmatc(ich,idx*4-4+4,k2*2-1)=recid(j)
              recmatc(ich,idx*4-4+4,k2*2)=recmatp(ich,idd*4-4+xr*2,k1*2)
              recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1)
              recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)
              if (recmatc(ich,idx*4-4+4,k2*2)==nmch) then
                exit       !finish
              end if
            end if
            if (recmatc(ich,idx*4-4+4,k2*2)+1==recmatp(ich,idd*4-4+xr*2,k1*2-1)) then
                k2=k2+1
                recmatc(ich,idx*4-4+4,k2*2-1)=recmatp(ich,idd*4-4+xr*2,k1*2-1) !copy
                recmatc(ich,idx*4-4+4,k2*2)=recmatp(ich,idd*4-4+xr*2,k1*2)    !copy
                recmatc(ich,idx*4-4+3,k2*2-1)=recmatp(ich,idd*4-4+xr*2-1,k1*2-1) !copy
                recmatc(ich,idx*4-4+3,k2*2)=recmatp(ich,idd*4-4+xr*2-1,k1*2)    !copy
              if (recmatc(ich,idx*4-4+4,k2*2)==nmch) then
                exit       !finish
              end if
            end if
          end do

          matcn(ich,idx*2)=k2

1002      continue

          !check errors ***************************************************
          do i=1,matcn(ich,idx*2)-1
            if (recmatc(ich,idx*4,i*2)+1.ne.recmatc(ich,idx*4,(i+1)*2-1)) then
              print*,'xr=',xr,'recid:',recid(1:recn)
              print*,'idd:',matpn(ich,idd*2-1),matpn(ich,idd*2)
              print*,recmatp(ich,idd*4-4+1,1:50)
              print*,recmatp(ich,idd*4-4+2,1:50)
              print*,recmatp(ich,idd*4-4+3,1:50)
              print*,recmatp(ich,idd*4-4+4,1:50)
              print*,''
              print*,'matn=',matcn(ich,idx*2)
              print*,recmatc(ich,idx*4-4+3,1:50)
              print*,recmatc(ich,idx*4-4+4,1:50)
              pause
            end if
          end do !*********************************************************
          if (recmatc(ich,idx*4,matcn(ich,idx*2)*2).ne.nmch) then
              print*,'xr=',xr,'recid:',recid(1:recn)
              print*,'idd:',matpn(ich,idd*2-1),matpn(ich,idd*2)
              print*,recmatp(ich,idd*4-4+1,1:50)
              print*,recmatp(ich,idd*4-4+2,1:50)
              print*,recmatp(ich,idd*4-4+3,1:50)
              print*,recmatp(ich,idd*4-4+4,1:50)
              print*,''
              print*,'matn=',matcn(ich,idx*2)
              print*,recmatc(ich,idx*4-4+3,1:50)
              print*,recmatc(ich,idx*4-4+4,1:50)
              pause
          end if

          !determine how manay times mutation occur - mutn *****************
          call random_number (rnd)
          do mutn=0,prn
            if (rnd < mutpr(mutn)) exit
          end do

          do i=1,mutn       !no. mutations per chr.
147         call random_number(rnd)
            iv1=int(nmch*rnd)+1     !nmch - no. markers in one chr.
            do j=1,i-1
              if (mutid(j)==iv1) goto 147    !mutated ID
            end do
            mutid(i)=iv1
            mutcn(ich,idx*2)=mutcn(ich,idx*2)+1
            mutmatc(ich,idx*2,mutcn(ich,idx*2))=iv1
          end do

          !if (igen==2) then
          !  print*,'xr=',xr,'recid:',recid(1:recn)
          !  print*,'mut:',mutid(1:mutn)
          !  print*,'idd mut no:',mutpn(ich,idd*2-1),mutpn(ich,idd*2)
          !  print*,mutmatp(ich,idd*2-2+1,1:50)
          !  print*,mutmatp(ich,idd*2-2+2,1:50)
          !  print*,'idx mut no:',mutcn(ich,idx*2)
          !  print*,mutmatc(ich,idx*2,1:50)
          !  pause
          !end if


        end do !ich

        !print*,'idx:',matcn(ich,idx*2-1),matcn(ich,idx*2)
        !print*,recmatc(ich,idx*4-4+1,1:50)
        !print*,recmatc(ich,idx*4-4+2,1:50)
        !print*,recmatc(ich,idx*4-4+3,1:50)
        !print*,recmatc(ich,idx*4-4+4,1:50)
        !pause


        !do m1=1,int(nmarkers*mr)
        !  if (genom1(idx*2,mutid(m1))=='1') then
        !    genom1(idx*2,mutid(m1))='2'
        !  elseif (genom1(idx*2,mutid(m1))=='2') then
        !    genom1(idx*2,mutid(m1))='1'
        !  end if
        !end do 



      end if
    end do

    !do i=1,20
    !  print*,recmatp(ich,i,1:50)
    !end do
    !print*,''
    !do i=1,20
    !  print*,recmatc(ich,i,1:50)
    !end do
    !print*,''
    !print*,matcn(ich,1:10)
    !pause

    !for recombinaiton ************************************
    matpn=matcn
    do i=1,(nmal+nfem) !neff
      do ich=1,nch
        recmatp(ich,i*4-4+1,1:matpn(ich,i*2-1)*2)=recmatc(ich,i*4-4+1,1:matpn(ich,i*2-1)*2)
        recmatp(ich,i*4-4+2,1:matpn(ich,i*2-1)*2)=recmatc(ich,i*4-4+2,1:matpn(ich,i*2-1)*2)
        recmatp(ich,i*4-4+3,1:matpn(ich,i*2)*2)=recmatc(ich,i*4-4+3,1:matpn(ich,i*2)*2)
        recmatp(ich,i*4-4+4,1:matpn(ich,i*2)*2)=recmatc(ich,i*4-4+4,1:matpn(ich,i*2)*2)
      end do
    end do

    !for mutations ****************************************
    mutpn=mutcn
    do i=1,(nmal+nfem) !neff
      do ich=1,nch
        mutmatp(ich,i*2-2+1,1:mutpn(ich,i*2-1))=mutmatc(ich,i*2-2+1,1:mutpn(ich,i*2-1))
        mutmatp(ich,i*2-2+2,1:mutpn(ich,i*2))=mutmatc(ich,i*2-2+2,1:mutpn(ich,i*2))
      end do
    end do

   !make a new base generation **********************************************
   if (mod(igen,50)==0) then       !every 100 generation
     do i=1,(nmal+nfem) !neff
       do ich=1,nch
         do j=1,matpn(ich,i*2-1)
           nk1=recmatp(ich,i*4-4+1,j*2-1)         !ID
           nk2=recmatp(ich,i*4-4+1,j*2)           !M or P
           k1=recmatp(ich,i*4-4+2,j*2-1)          !first marker ID
           k2=recmatp(ich,i*4-4+2,j*2)            !last marker ID
           genom1(i*2-1,(ich-1)*nmch+k1:(ich-1)*nmch+k2)=genom(nk1*2-2+nk2,(ich-1)*nmch+k1:(ich-1)*nmch+k2)
         end do

         !mutation ************************************************
         do j=1,mutpn(ich,i*2-1)
           k1=mutmatp(ich,i*2-1,j)
           if (genom1(i*2-1,(ich-1)*nmch+k1)=='1') then
             genom1(i*2-1,(ich-1)*nmch+k1)='2'
           elseif (genom1(i*2-1,(ich-1)*nmch+k1)=='2') then
             genom1(i*2-1,(ich-1)*nmch+k1)='1'
           end if
         end do

         do j=1,matpn(ich,i*2)
           nk1=recmatp(ich,i*4-4+3,j*2-1)         !ID
           nk2=recmatp(ich,i*4-4+3,j*2)           !M or P
           k1=recmatp(ich,i*4-4+4,j*2-1)          !first marker ID
           k2=recmatp(ich,i*4-4+4,j*2)            !last marker ID
           genom1(i*2,(ich-1)*nmch+k1:(ich-1)*nmch+k2)=genom(nk1*2-2+nk2,(ich-1)*nmch+k1:(ich-1)*nmch+k2)
         end do
         !mutation ************************************************
         do j=1,mutpn(ich,i*2)
           k1=mutmatp(ich,i*2,j)
           if (genom1(i*2,(ich-1)*nmch+k1)=='1') then
             genom1(i*2,(ich-1)*nmch+k1)='2'
           elseif (genom1(i*2,(ich-1)*nmch+k1)=='2') then
             genom1(i*2,(ich-1)*nmch+k1)='1'
           end if
         end do
       end do  !ich
     end do !i
     !genom(:,1:nmch)=genom1(:,1:nmch)
     genom=genom1

     !initializing for recmat and mutmat
     !recombination matrix
     do ich=1,nch
       matpn(ich,:)=1
     end do
     do i=1,(nmal+nfem) !neff
       do ich=1,nch
         recmatp(ich,i*4-4+1,1)=i
         recmatp(ich,i*4-4+2,1)=1
         recmatp(ich,i*4-4+3,1)=i
         recmatp(ich,i*4-4+4,1)=1
         recmatp(ich,i*4-4+1,2)=1
         recmatp(ich,i*4-4+2,2)=nmch
         recmatp(ich,i*4-4+3,2)=2
         recmatp(ich,i*4-4+4,2)=nmch
       end do
     end do

     !mutation matrix
     do ich=1,nch
       mutpn(ich,:)=0
     end do

   end if

  END do   !igen ~ ngen
  !PRINT*,'mutation - size'
  !pause

  mafreq1='0'
  !print*,'finding no. marker alleles ***********************'
  do im=1,nmarkers
    !print*,im
    ma1=0;x1=0;x2=0
    do i=1,(nmal+nfem)*2 !neff*2
      if (genom(i,im).ne.'0') then
        x2=x2+1
        do j=1,ma1
          if (genom(i,im)==mafreq1(im,j)) then
            if (j==1) x1=x1+1
            goto 1724 !then
          end if
        end do ! j
        ma1=ma1+1
        !print*,ma1
        mafreq1(im,ma1)=genom(i,im)
        if (ma1==1) x1=x1+1
        !PRINT*,im,ma1,genom2or(i,:)
        !pause
1724    continue
      end if
    end do ! i
    mano(im)=ma1
    if (ma1.ne.2) then
      !print*,'allele no. not 2 >>> check',ma1,mafreq1(im,1:ma1)
      !pause
    end if
    mafreq2(im,1)=x1/x2
    mafreq2(im,2)=1-x1/x2
    !write(75,*)'marker ',im,' - ',ma1,' alleles',mafreq1(im,0)
    !pause
    !*********************************
  end do !im

  open (unit=111,file=trim(sim_coal_file)//'.frq',status='unknown')
  do i=1,nmarkers
    write(111,*)i,'  ',mafreq1(i,1),' ',mafreq1(i,2),'  ',mafreq2(i,:)
  end do
  close(111)


  deallocate(nselmal,nselfem,genom1)

  allocate(nselmal(neff*2*empsz),nselfem(neff*2*empsz))
  allocate(genom1(2,nmarkers))

  !print*,'generating current population to be investigated, x ',empsz,' *****'
  ncont10=0
  nk1=0;nk2=0
  !k=0
  do i=1,popsiz*2*empsz
    call random_number(rnd)
    nsel=int(rnd*nmal)+1     !random selection and mating 
    nselmal(i)=nsel
    !k=k+1
    !nselmal(i)=k
    !if (k==nmal) k=0
  end do
  !k=0
  do i=1,popsiz*2*empsz
    call random_number(rnd)
    nsel=nmal+int(rnd*nfem)+1
    nselfem(i)=nsel
    !k=k+1
    !nselfem(i)=nmal+k
    !if (k==nfem) k=0
  end do

  !creat progeny
  do ipar1=1,popsiz*empsz
    ids=nselmal(ipar1)
    idd=nselfem(ipar1)
    if (nk1.ne.popmal*empsz.or.nk2.ne.popfem*empsz) then
225   call random_number(rnd)
      if (rnd.lt.0.5.and.nk1.lt.popmal*empsz) then
        nk1=nk1+1
        idx=nk1
      ELSEIF (rnd.ge.0.5.and.nk2.lt.popfem*empsz) then
        nk2=nk2+1
        idx=popmal*empsz+nk2
      else
        GOTO 225
      END if

      ! pedigree recording to calculate h2
      ped(idx,1)=ids
      ped(idx,2)=idd

        !From dad**********************************************************
        do ich=1,nch
          !print*,ich
          call random_number(rnd)
          xr=INT(2*rnd)+1         !random determination

          !determine how manay times rec. occur - recn *********************
          call random_number (rnd)
          do recn=0,prn    !10
            if (rnd < binpr(recn)) exit
          end do
          !print*,recn

          !simplified - one recombination per chr **************************
          if (spf==1) recn=1
          !*****************************************************************

          do i=1,recn
155         call random_number (rnd)
            !print*,rnd
            iv1=int((nmch-1)*rnd)+2   !avoiding 1 
            do j=1,i-1
              if (recid(j)==iv1) goto 155    !recombined ID
            end do
            recid(i)=iv1
          end do
          !sorting************************'
          deallocate(srtv)
          allocate(srtv(recn,2))
          do i=1,recn
            srtv(i,2)=recid(i)
          end do
          call sort (srtv,recn)
          do i=1,recn
            recid(i)=srtv(i,2)
          end do !**********************************************************

          j=1     !indication for recombined marker ID
          do im=1,nmch       ! markers + QTL no.
            if (recn==0.or.j>recn) goto 501
            if (recid(j)==im) then
              xr=3-xr
              j=j+1
            end if
501         continue
            genom1(1,(ich-1)*nmch+im)=genom(ids*2-(2-xr),(ich-1)*nmch+im)
          end do


          !determine how manay times mutation occur - mutn *****************
          call random_number (rnd)
          do mutn=0,prn
            if (rnd < mutpr(mutn)) exit
          end do

          do i=1,mutn       !no. mutations per chr.
153         call random_number(rnd)
            iv1=int(nmch*rnd)+1     !nmch - no. markers in one chr.
            do j=1,i-1
              if (mutid(j)==iv1) goto 153    !mutated ID
            end do
            mutid(i)=iv1
          end do

          do j=1,mutn
            k1=mutid(j)
            if (genom1(1,(ich-1)*nmch+k1)=='1') then
              genom1(1,(ich-1)*nmch+k1)='2'
            elseif (genom1(1,(ich-1)*nmch+k1)=='2') then
              genom1(1,(ich-1)*nmch+k1)='1'
            end if
          end do

        end do    !ich


        !print*,'!From mam********************************************'
        do ich=1,nch
          !print*,'mum',ich
          call random_number(rnd)
          xr=INT(2*rnd)+1         !random determination

          !determine how manay times rec. occur - recn *********************
          call random_number (rnd)
          do recn=0,prn     !10
            if (rnd < binpr(recn)) exit
          end do
          !print*,recn

          !simplified - one recombination per chr **************************
          if (spf==1) recn=1
          !*****************************************************************

          do i=1,recn
159         call random_number (rnd)
            iv1=int((nmch-1)*rnd)+2   !avoiding 1 
            !print*,i,nmch,iv1,recn,recid(1:i-1)
            do j=1,i-1
              if (recid(j)==iv1) goto 159    !recombined ID
            end do
            recid(i)=iv1
          end do
          !sorting************************
          deallocate(srtv)
          allocate(srtv(recn,2))
          do i=1,recn
            srtv(i,2)=recid(i)
          end do
          call sort (srtv,recn)
          do i=1,recn
            recid(i)=srtv(i,2)
          end do !**********************************************************

          j=1     !indication for recombined marker ID
          do im=1,nmch       ! markers + QTL no.
            if (recn==0.or.j>recn) goto 401
            if (recid(j)==im) then
              xr=3-xr
              j=j+1
            end if
401         continue
            genom1(2,(ich-1)*nmch+im)=genom(idd*2-(2-xr),(ich-1)*nmch+im)
          end do

          !determine how manay times mutation occur - mutn *****************
          call random_number (rnd)
          do mutn=0,prn
            if (rnd < mutpr(mutn)) exit
          end do

          do i=1,mutn       !no. mutations per chr.
157         call random_number(rnd)
            iv1=int(nmch*rnd)+1     !nmch - no. markers in one chr.
            do j=1,i-1
              if (mutid(j)==iv1) goto 157    !mutated ID
            end do
            mutid(i)=iv1
          end do

          do j=1,mutn
            k1=mutid(j)
            if (genom1(2,(ich-1)*nmch+k1)=='1') then
              genom1(2,(ich-1)*nmch+k1)='2'
            elseif (genom1(2,(ich-1)*nmch+k1)=='2') then
              genom1(2,(ich-1)*nmch+k1)='1'
            end if
          end do


        end do  !ich

      end if !(nk1.ne.popmal*empsz ~ ...)


      ncont10=ncont10+1
      contped(ncont10,1)=idx
      contped(ncont10,2)=ids
      contped(ncont10,3)=idd
      contom(ncont10*2-1,:)=genom1(1,:)   !paternal genotype
      contom(ncont10*2,:)=genom1(2,:)     !maternal genotype
      !print*,ncont10,ipar1,popsiz,popsiz*2*empsz,neff


  end do  !ipar1


  !recording without causal mutation and < freq. of 0.01 
  k=0
  do im=1,nmarkers
    !if (mafreq2(im,1)<0.99.and.mafreq2(im,1)>0.01) then
      k=k+1
      midx(k)=im
    !end if
  
1471   continue

  end do
  midxn=k
  !print*,midxn

  !open (unit=101,file='simdis.ped',status='unknown')
  open (unit=101,file=trim(sim_coal_file)//'.ped',status='unknown')
  do k=1,ncont10
    write(101,'(2I10,4I3,1000000a2)')contped(k,1),contped(k,1),0,0,2,-9,'  ',(contom(k*2-1,midx(j)),' ',contom(k*2,midx(j)),'  ',j=1,midxn)
    !write(101,'(4I10,2I3,1000000a2)')contped(k,1),contped(k,:),2,-9,'  ',(contom(k*2-1,midx(j)),' ',contom(k*2,midx(j)),'  ',j=1,midxn)

  end do
  close(101)

  !open (unit=102,file='simdis.map',status='unknown')
  open (unit=102,file=trim(sim_coal_file)//'.map',status='unknown')
  v1=0;k=0;iv1=1
  do i=1,nmarkers
    v1=v1+1
    if (mod(i-1,nmch)==0) then
      v1=1;k=k+1
    end if
    if (midx(iv1)==i) then
      iv1=iv1+1
      write(102,'(*(G0,2X))')k,i,v1*(10000000*chl/(nmch*0.1)),v1*(10000000*chl/(nmch*0.1))!,mafreq1(i,1)," ",mafreq1(i,2)
    end if
  end do
  close(102)


end subroutine



