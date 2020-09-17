program mtg2

!**********************************************************************
! MTG2 (MT GREML, MT GBLUP, RRM, RNM)
! S. Hong Lee  (2009 - 2019)
!************************************************************************

use fileio
implicit none
integer::n,pedn,rn,xi,i,j,k,k2,kk,jj,tn,zi,tfn,trn,thn,pca,bend,inv,chol,mrrm,mrnm
integer,allocatable::fnm(:),rnm(:)
character(len=64),allocatable::grm_known(:),tzb(:)
double precision::ve,va,gasdev,h2,bva,bve,bh2,ncase,ncont,rnd

character (len=1024) :: linbuf                !for gz grm
type (ioport) :: port                         !for gz grm

character(len=100),allocatable::pedor(:,:),obs(:,:)
character(len=100)::cdum1,cdum2
character(len=64),allocatable::cdum_t(:) !for each trait

!double precision,allocatable::GRM(:,:,:),eval(:,:),y(:,:),idx_y(:,:)
real,allocatable::GRM(:,:,:),eval(:,:),y(:,:),idx_y(:,:)
double precision::x1,v1,conv
real::r1,rtmx
!double precision::r1

integer::io,ioerr,mn,narg,mg,bmg,zmg,tmg,cc,qc,cc_c,qc_c,matmat,t2b
integer,allocatable::idx_ped(:,:),idx_obs(:,:),yidx(:,:)
integer,allocatable::idx_in(:),inlst(:),outlst(:)

double precision::irgt
integer::icove,rrme,nit,frq,aug,OMP_GET_NUM_THREADS,gwas
character(len=64)::fl1,fl2,filnam,fl3,fl4,fl24,fl5,fl6,fl7,hrtmx,wt_par,wt_res
character(len=64)::fl8,fl9,fl10,fl11,fl11r,fl12,fl_rrm,fl_spl,plink,vgpy,sbv
character(len=64)::fl_pred,fl_rnm,fl_mcv,matvec
character(len=64)::delta,delta2,var_rt,pdmx,simcoal,simreal,simcoal2,simstp
character(len=64),allocatable::ccm(:,:),ccn(:,:,:),tmp_qcm(:)
real,allocatable::qcm(:,:)

character(len=3)::cdum4

integer*1,allocatable::sex(:)

!random number generator****************************************
integer::seed2,rseed
integer,allocatable::seed(:)
double precision::dseed

character(len=10000)::string

double precision,allocatable::xmm(:,:,:)
integer,allocatable::cc_zi(:,:),ccdet(:),ccdet2(:),ccdet3(:),qcdet(:),cc_lvl(:)

!time measure **********************************************************
INTEGER::now(8)
CHARACTER*8 date
CHARACTER*20 time,zone
double precision:: timef, t1_real, t2_real, elapsed_real_secs
double precision:: t1_cpu, t2_cpu, elapsed_cpu_secs
!***********************************************************************

real::h2_a,Me_a,k_a,p_a,p2_a,ckv_a
real::h2_b,k_b,p_b,p2_b,se_b
character(len=64)::indc
integer::N_a,M_a
real::t1_h2,t1_k,t1_p,se_c,rg_c,t2_h2,t2_k,t2_p

real::len
integer::chrn,ne

logical::file_exist

print*,"******************************************************************"
print*,"MTG2 version 2.15 (June2019)"
print*,"******************************************************************"

!prediction accuracy ***************************************************
h2_a=-99999
call get_command_argument(1, filnam)
if (filnam=='-pred_acc') then     !pred_acc
  print*,''
  print*,"Expected prediction accuracy given the following parameters *****" 
  narg=command_argument_count()
  if (narg.ne.9) then
    print*,"h2, N, M, Me, k, p, p2 and ckv should be specified"
    stop
  end if
  call get_command_argument(2, filnam)
      read(filnam,*)h2_a
      print*,'h2  : ',h2_a
  call get_command_argument(3, filnam)
      read(filnam,*)N_a
      print'(*(G0,2X))',' N   : ',N_a
  call get_command_argument(4, filnam)
      read(filnam,*)M_a
      print'(*(G0,2X))',' M   : ',M_a
  call get_command_argument(5, filnam)
      read(filnam,*)Me_a
      print*,'Me  : ',Me_a
  call get_command_argument(6, filnam)
      read(filnam,*)k_a
      print*,'k   : ',k_a
  call get_command_argument(7, filnam)
      read(filnam,*)p_a
      print*,'p   : ',p_a
  call get_command_argument(8, filnam)
      read(filnam,*)p2_a
      print*,'p2  : ',p2_a
  call get_command_argument(9, filnam)
      read(filnam,*)ckv_a
      print*,'ckv : ',ckv_a
end if  !*****************************************************************

fl_pred='null'
if (filnam=='-pred_acc2') then     !spline random regression ***
  call get_command_argument(2, filnam)
    fl_pred=filnam
    print*,'pred par file : ',trim(fl_pred)
end if


ne=-99999
call get_command_argument(1, filnam)
if (filnam=='-Me') then     !pred_acc
  print*,''
  print*,'Effective number of chromosome segment given the following parameters *****'
  narg=command_argument_count()
  if (narg.ne.4) then
    print*,"Ne, length and number of chromosome should be specified"
    stop
  end if
  call get_command_argument(2, filnam)
      read(filnam,*)ne
      print'(*(G0,2X))',' Effective pop. size           : ',ne
  call get_command_argument(3, filnam)
      read(filnam,*)len
      print*,'Genomic length per each chr   : ',len 
  call get_command_argument(4, filnam)
      read(filnam,*)chrn
      print'(*(G0,2X))',' Number of chr                 : ',chrn
end if


!trf function ***************************************************
h2_b=-99999
call get_command_argument(1, filnam)
if (filnam=='-trf_h2') then     !trf
  print*,''
  print*,"Transformation between observed and liability scale ***************" 
  narg=command_argument_count()
  if (narg.ne.6) then
    print*,"h2, se, k, p, and scale (obs or liab) should be specified"
    stop
  end if
  call get_command_argument(2, filnam)
      read(filnam,*)h2_b
      print*,'h2    : ',h2_b 
  call get_command_argument(3, filnam)
      read(filnam,*)se_b
      print*,'se    : ',se_b 
  call get_command_argument(4, filnam)
      read(filnam,*)k_b
      print*,'k     : ',k_b
  call get_command_argument(5, filnam)
      read(filnam,*)p_b
      print*,'p     : ',p_b
  call get_command_argument(6, indc)
      print*,'scale : ',trim(indc)
end if


!trf_rg function ***************************************************
rg_c=-99999
call get_command_argument(1, filnam)
if (filnam=='-trf_rg') then     !trf_rg
  print*,''
  print*,"Transformation to liability scale for genetic covariance *******" 
  print*,"i.e. SNP-coheritability                                  *******" 
  narg=command_argument_count()
  if (narg.ne.9) then
    print*,"t1_h2, t1_k, t1_p, t2_h2, t2_k, t2_p, rg, se should be specified"
    stop
  end if
  call get_command_argument(2, filnam)
      read(filnam,*)t1_h2
      print*,'t1_h2    : ',t1_h2 
  call get_command_argument(3, filnam)
      read(filnam,*)t1_k
      print*,'t1_k    : ',t1_k 
  call get_command_argument(4, filnam)
      read(filnam,*)t1_p
      print*,'t1_p    : ',t1_p 
  call get_command_argument(5, filnam)
      read(filnam,*)t2_h2
      print*,'t2_h2    : ',t2_h2 
  call get_command_argument(6, filnam)
      read(filnam,*)t2_k
      print*,'t2_k    : ',t2_k 
  call get_command_argument(7, filnam)
      read(filnam,*)t2_p
      print*,'t2_p    : ',t2_p 
  call get_command_argument(8, filnam)
      read(filnam,*)rg_c
      print*,'rg    : ',rg_c 
  call get_command_argument(9, filnam)
      read(filnam,*)se_c
      print*,'se    : ',se_c 
end if



narg=command_argument_count()
if (narg==0.or.narg==1) then
  print*,'-p fam file -d dat file -g grm file ...'
  stop
end if

fl1='null'
fl2='null'
fl5='ascm.out'
fl3='null';mg=-1;bmg=-1;zmg=-1;tmg=-1
fl11='null'
fl11r='null'
fl12='null'
fl4='null'
fl24='null'
aug=0
cc=0;qc=0
tn=2         !defoult: 2 traits 
icove=0       !defoult: no cove
irgt=-9      !defoult: no rg test
thn=1        !default: no. threads 
fl_rrm='null'       !rrm par file
mrrm=0              !# tratis for mrrm
mrnm=0              !# tratis for mrnm
fl_spl='null'       !rrm par file
rrme=0       !default: no residual random regression 
fl_rnm='null'       !rnm par file
fl_mcv='null'       !multiple covariates for rnm - par file
t2b=0        !converting GRM format between binary and txt
pca=0        !default: no PCs caculation
bend=0        !default: no bending matrix
inv=0        !default: no inverting matrix
chol=0        !default: no cholesky decomposition
matvec='null'        !default: no mat x vec
matmat=0        !default: no mat x mat 
nit=200      !default: maximum no iterations
conv=0.001   !default: convergence criteria
plink='null'
gwas=0
frq=0        !allele frequency
rtmx=-999    !relationship matrix scale factor
wt_par='null'    !weighting factor parameters
wt_res='null'    !weighting factor for residual structure
hrtmx='null'    !H matrix 
sbv='null'   !snp_blup option, a or b
vgpy='null'
delta='null'
delta2='null'
simcoal='null'
simcoal2='null'
simreal='null'
simstp='null'
var_rt='null'
pdmx='null'  !production matrix to fit random effects

do i=1, (narg/2)
  call get_command_argument(i*2-1, filnam)
  if (filnam=='-p'.or.filnam=='-fam') then
    call get_command_argument(i*2, filnam)
      fl1=filnam
      print*,'ID file   : ',trim(fl1)
  end if
  if (filnam=='-d'.or.filnam=='-pheno') then
    call get_command_argument(i*2, filnam)
      fl2=filnam
      print*,'dat file  : ',trim(fl2)
  end if
  if (filnam=='-g'.or.filnam=='-grm') then
    mg=0;tmg=0
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'grm file  : ',trim(fl3)
  end if
  if (filnam=='-bg') then
    bmg=0;tmg=0
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'bgrm file : ',trim(fl3)
  end if
  if (filnam=='-zg') then
    zmg=0;tmg=0
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'zgrm file  : ',trim(fl3)
  end if

  if (filnam=='-mg') then
    mg=1
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'grm file  : ',trim(fl3)
  end if
  if (filnam=='-mbg') then
    bmg=1
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'bgrm file  : ',trim(fl3)
  end if
  if (filnam=='-mzg') then
    zmg=1
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'zgrm file  : ',trim(fl3)
  end if
  if (filnam=='-tmg') then
    tmg=1
    call get_command_argument(i*2, filnam)
      fl3=filnam
      print*,'grm (mixed) file  : ',trim(fl3)
  end if
  if (filnam=='-eig') then
    call get_command_argument(i*2, filnam)
      fl12=filnam
      print*,'eig prefix: ',trim(fl12)
  end if

  if (filnam=='-sv') then
    call get_command_argument(i*2, filnam)
      fl4=filnam
      print*,'sv file   : ',trim(fl4)
  end if
  if (filnam=='-fix') then
    call get_command_argument(i*2, filnam)
      fl24=filnam
      print*,'fix file  : ',trim(fl24)
  end if
  if (filnam=='-out') then
    call get_command_argument(i*2, filnam)
      fl5=filnam
      print*,'out file  : ',trim(fl5)
  end if
  if (filnam=='-cc') then
    cc=1
    call get_command_argument(i*2, filnam)
      fl6=filnam
      print*,'cc file   : ',trim(fl6)
  end if
  if (filnam=='-qc') then
    qc=1
    call get_command_argument(i*2, filnam)
      fl7=filnam
      print*,'qc file   : ',trim(fl7)
  end if
  if (filnam=='-mod') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      fl8=filnam
      print*,'nr mode   : ',trim(fl8)
      read(fl8,*)tn
  end if
  if (filnam=='-cove') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      fl9=filnam
      print*,'cove      : ',trim(fl9)
      read(fl9,*)icove
  end if
  if (filnam=='-rg') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      fl10=filnam
      print*,'rg test   : ',trim(fl10)
      read(fl10,*)irgt
  end if
  if (filnam=='-bv') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      fl11=filnam
      print*,'ebv output: ',trim(fl11)
  end if
  if (filnam=='-bvr') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      fl11r=filnam
      print*,'ebv output: ',trim(fl11r)
  end if
  if (filnam=='-thread') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)thn
      print*,'threads   : ',thn
  end if
  if (filnam=='-mrrm') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)mrrm
      print*,'mrrm trait: ',mrrm
  end if
  if (filnam=='-mrnm') then     !number of traits ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)mrnm
      print*,'mrnm trait: ',mrnm
  end if
  if (filnam=='-rrm') then     !random regression ***
    call get_command_argument(i*2, filnam)
      fl_rrm=filnam
      print*,'rrm file  : ',trim(fl_rrm)
  end if
  if (filnam=='-spl') then     !spline random regression ***
    call get_command_argument(i*2, filnam)
      fl_spl=filnam
      print*,'spl file  : ',trim(fl_spl)
  end if
  if (filnam=='-rrme') then     !random regression with error term too  ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)rrme
      print*,'rrm for e : ',trim(filnam)
  end if
  if (filnam=='-rnm') then     !reaction norm for continous variable ***
    call get_command_argument(i*2, filnam)
      fl_rnm=filnam
      print*,'rnm file  : ',trim(fl_rnm)
  end if
  if (filnam=='-mcv') then     !reaction norm for continous variable ***
    call get_command_argument(i*2, filnam)
      fl_mcv=filnam
      print*,'mcv file  : ',trim(fl_mcv)
  end if

  if (filnam=='-t2b') then     
    call get_command_argument(i*2, filnam)
      read(filnam,*)t2b
      print*,'t2b       : ',trim(filnam)
  end if
  if (filnam=='-pca') then     !PCs ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)pca
      print*,'pca no.   : ',trim(filnam)
  end if
  if (filnam=='-bend') then     !PCs ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)bend
      print*,'bending   : ',trim(filnam)
  end if
  if (filnam=='-inv') then     !PCs ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)inv
      print*,'inverting : ',trim(filnam)
  end if
  if (filnam=='-chol') then     !PCs ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)chol
      print*,'cholesky  : ',trim(filnam)
  end if
  if (filnam=='-matvec') then     !PCs ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)matvec
      print*,'matvec    : ',trim(filnam)
  end if
  if (filnam=='-matmat') then     !PCs ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)matmat
      print*,'matmat    : ',trim(filnam)
  end if

  if (filnam=='-conv') then     !convergence criteria ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)conv
      print*,'conv      : ',trim(filnam)
  end if
  if (filnam=='-nit') then     !maximu no iterations ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)nit
      print*,'max nit   : ',trim(filnam)
  end if
  if (filnam=='-plink') then     !plink files prefix ***
    call get_command_argument(i*2, plink)
      print*,'plink     : ',trim(plink)
  end if
  if (filnam=='-gwas') then     !plink files prefix ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)gwas
      print*,'gwas      : ',trim(filnam)
  end if
  if (filnam=='-frq') then     !allele frequencies and var(x) ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)frq
      print*,'freq      : ',trim(filnam)
  end if
  if (filnam=='-rtmx') then     !scale parameter for grm ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)rtmx
      print*,'rtmx      : ',trim(filnam)
  end if
  if (filnam=='-wt') then     !scale parameter for grm ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)wt_par
      print*,'wt par    : ',trim(filnam)
  end if
  if (filnam=='-wtr') then     !scale parameter for grm ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)wt_res
      print*,'wt res    : ',trim(filnam)
  end if
  if (filnam=='-hrtmx') then     !ped info for h matrix (ssgblup) ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)hrtmx
      print*,'hrtmx      : ',trim(filnam)
  end if
  if (filnam=='-simreal') then     !allele frequencies and var(x) ***
    call get_command_argument(i*2, simreal)
      print*,'simreal      : ',trim(simreal)
  end if

  if (filnam=='-sbv') then     !snp_blup option ***
    call get_command_argument(i*2, sbv)
      print*,'snp_blup  : ',trim(sbv)
  end if
  if (filnam=='-vgpy') then     !snp_blup input ***
    call get_command_argument(i*2, vgpy)
      print*,'VgPy      : ',trim(vgpy)
  end if
  if (filnam=='-delta') then     !delta input ***
    call get_command_argument(i*2, delta)
      print*,'delta     : ',trim(delta)
  end if
  if (filnam=='-delta2') then     !delta input ***
    call get_command_argument(i*2, delta2)
      print*,'delta2    : ',trim(delta2)
  end if
  if (filnam=='-simcoal') then     !simulation input ***
    call get_command_argument(i*2, simcoal)
      print*,'simcoal   : ',trim(simcoal)
  end if
  if (filnam=='-simcoal2') then     !simulation input ***
    call get_command_argument(i*2, simcoal2)
      print*,'simcoal2  : ',trim(simcoal2)
  end if
  if (filnam=='-simstp') then     !simulation input ***
    call get_command_argument(i*2, simstp)
      print*,'simstp  : ',trim(simstp)
  end if
  if (filnam=='-var_rt') then     !var(relatioship) input ***
    call get_command_argument(i*2, var_rt)
      print*,'var_rt    : ',trim(var_rt)
  end if
  if (filnam=='-pdmx') then     !product matrix inpur ***
    call get_command_argument(i*2, pdmx)
      print*,'pdmx      : ',trim(pdmx)
  end if


  if (filnam=='-aug') then     !maximu no iterations ***
    call get_command_argument(i*2, filnam)
      read(filnam,*)aug
      print*,'augment y : ',trim(filnam)
  end if

end do
print*,''



!simple functions *****************************************************

!delta **********************************************************
if (delta.ne.'null') then
  call delta_ (delta)   !snp_blup
  goto 1100
end if
if (delta2.ne.'null') then
  call delta_2 (delta2)   !snp_blup
  goto 1100
end if

if (simcoal.ne.'null') then
  call sim_coal_sub (simcoal)   !snp_blup
  goto 1100
end if
if (simcoal2.ne.'null') then
  call sim_coal_sub2 (simcoal2)   !snp_blup
  goto 1100
end if

if (h2_a.ne.-99999) then
  call pred_acc (h2_a,N_a,M_a,Me_a,k_a,p_a,p2_a,ckv_a)   !snp_blup
  goto 1100
end if

if (fl_pred.ne.'null') then
  call pred_acc2 (fl_pred)   !snp_blup
  goto 1100
end if

if (ne.ne.-99999) then
  call Me (ne,len,chrn)   !snp_blup
  goto 1100
end if

if (h2_b.ne.-99999) then
  call trf_h2 (h2_b,se_b,k_b,p_b,indc)   !snp_blup
  goto 1100
end if

if (rg_c.ne.-99999) then
  call trf_rg (t1_h2,t1_k,t1_p,t2_h2,t2_k,t2_p,rg_c,se_c)   !snp_blup
  goto 1100
end if


!check
if (pca.ne.0 .and. (mg==1 .or. tmg==1 .or. bmg==1 .or. zmg==1)) then
  print*,'with -pca, -g, -bg or -zg should be specified'
  pause
end if 
if (bend.ne.0 .and. (mg==1 .or. tmg==1 .or. bmg==1 .or. zmg==1)) then
  print*,'with -bend, -g, -bg or -zg should be specified'
  pause
end if 
if (inv.ne.0 .and. (mg==1 .or. tmg==1 .or. bmg==1 .or. zmg==1)) then
  print*,'with -inv, -g, -bg or -zg should be specified'
  pause
end if 
if (chol.ne.0 .and. (mg==1 .or. tmg==1 .or. bmg==1 .or. zmg==1)) then
  print*,'with -chol, -g, -bg or -zg should be specified'
  pause
end if 
if (matvec.ne.'null' .and. (mg==1 .or. tmg==1 .or. bmg==1 .or. zmg==1)) then
  print*,'with -matvec, -g, -bg or -zg should be specified'
  pause
end if 
if (matmat.ne.0 .and. (mg.ne.1 .and. tmg.ne.1 .and. bmg.ne.1 .and. zmg.ne.1)) then
  print*,'with -matmat, -tmg, -mg, -mbg or -mzg should be specified'
  pause
end if 

if (rrme.ne.0 .and. fl12.eq.'null') then
  print*,'with -rrme, -eig should be specified'
  pause
end if 

if (fl12.ne.'null' .and. fl3.ne.'null') then
  print*,'either only -eig or -g should be specified'
  pause
end if

if (fl12.ne.'null' .and. fl_rnm.ne.'null') then
  print*,'with -rnm, -eig should not be used'
  pause
end if

if (fl12.ne.'null' .and. fl_mcv.ne.'null') then
  print*,'with -mcv, -eig should not be used'
  pause
end if

if (fl_rnm.ne.'null' .and. fl_mcv.ne.'null') then
  print*,'-rnm and -mcv cannot be used together'
  pause
end if

if ((fl12.ne.'null'.or.fl_rrm.ne.'null'.or.fl_spl.ne.'null'.or.fl_rnm.ne.'null') .and. wt_res.ne.'null') then
  print*,'with -eig, -rrm, -spl or -rnm, -wtr is not supported'
  pause
end if

if (fl_rrm.ne.'null' .and. fl_spl.ne.'null') then
  print*,'either only -rrm or -spl should be specified'
  pause
end if

if (frq==1 .and. plink.eq.'null') then
  print*,'with -frq, -plink should be specified'
  pause
end if

if (simreal.ne.'null' .and. plink.eq.'null') then
  print*,'with -simreal, -plink should be specified'
  pause
end if
if (simstp.ne.'null' .and. plink.eq.'null') then
  print*,'with -simstp, -plink should be specified'
  pause
end if

if (rtmx.ne.-999 .and. plink.eq.'null') then
  print*,'with -rtmx, -plink should be specified'
  pause
end if

if (sbv.ne.'null' .and. vgpy.eq.'null') then
  print*,'with -sbv, -VgPy file should be specified'
  pause
end if


if (fl1.eq.'null' .and. plink.ne.'null') then  !plink fam file provided
  fl1=trim(plink)//".fam"
end if

if (fl1.eq.'null' .and. plink.eq.'null') then
  print*,'-p should be specified'
  pause
end if

if (aug.ne.0 .and. fl12.eq.'null') then
  print*,'with -aug, -eig should be specified'
  pause
end if 


call OMP_SET_NUM_THREADS(thn)

!print*,OMP_GET_NUM_THREADS()


n=0
open (unit=38,file=fl1,status='old')      !fam file
do 
  n=n+1
  read(38,*,iostat=io)cdum4
  if (io.ne.0) exit
  !print*,trim(pedor(i,1)),trim(pedor(i,2))
end do
close(38)
n=n-1
print*,'no. ID: ',n
rn=n          !should be the same (missing = -99999)


if (fl12.ne."null") then  !if eigenvector and eigenvalues to be used
  if (fl_rrm.ne."null" .and. rrme.ne.0) then
    mn=2
  elseif (fl_spl.ne."null" .and. rrme.ne.0) then
    mn=2
  elseif (fl_rrm.ne."null" .and. rrme.eq.0) then
    mn=1
  elseif (fl_spl.ne."null" .and. rrme.eq.0) then
    mn=1
  elseif (fl_rrm.eq."null") then
    mn=1
  elseif (fl_spl.eq."null") then
    mn=1
  else
    print*,"someting wrong >>> check, -rrme 1 should be only with -eig and -rrm (-spl)"
  end if
else

  mn=0
  if (fl3.ne.'null') then
    if (mg==0.or.tmg==0.or.bmg==0.or.zmg==0) then
      mn=1
    else if (mg==1.or.tmg==1.or.bmg==1.or.zmg==1) then
      mn=0
      open (unit=39,file=fl3,status='old')      !fam file
      do
        mn=mn+1
        read(39,*,iostat=io)cdum4
        if (io.ne.0) exit
      end do
      close(39)
      mn=mn-1
    else
      print *,"check >>>"
    end if
    print*,'no. grm: ',mn
  end if
end if

if (hrtmx .ne. 'null') then
  pedn=0
  open (unit=48,file=hrtmx,status='old')      !fam file
  do 
    pedn=pedn+1
    read(48,*,iostat=io)cdum4
    if (io.ne.0) exit
    !print*,trim(pedor(i,1)),trim(pedor(i,2))
  end do
  close(48)
  pedn=pedn-1
  print*,'no. ID in ped (hrtmx): ',pedn
end if  


allocate(pedor(n,4),y(rn,tn))
allocate(GRM(mn,n,n),grm_known(mn),tzb(mn),obs(rn,2),eval(n,mn))
ALLOCATE(idx_ped(n,4),idx_obs(rn,2),idx_y(rn,tn),yidx(tn,rn))
ALLOCATE(idx_in(rn),sex(n))
allocate(ccm(n,1000),qcm(n,1000),tmp_qcm(1000),ccdet(1000),qcdet(10000))
allocate(cdum_t(tn),cc_lvl(tn))
allocate(xmm(tn,n,1000),fnm(tn),rnm(tn),ccdet2(tn),ccdet3(tn))

! for pedigree and haplotypes*********************************** 
open (unit=38,file=fl1,status='old')      !fam file
do i=1,n
  read(38,*)pedor(i,:),sex(i),cdum4
  !print*,trim(pedor(i,1)),trim(pedor(i,2))
end do
close(38)

!check duplications in ped file **********************************
do i=1,n
  do j=i+1,n
    if (pedor(i,1)==pedor(j,1).and.pedor(i,2)==pedor(j,2)) then
      print*,'ID duplications in .ped: ',trim(pedor(i,1)),' ',trim(pedor(i,2))
      pause
    end if
  end do
end do !**********************************************************

!print*,1,rn
! for data**************************************
if (fl2.ne.'null') then 
  open (unit=40,file=fl2,status='old')      !dat file
  do i=1,rn
    read(40,*)obs(i,:),cdum_t(:)
    !print*,obs(i,:),cdum_t(:)
    do j=1,tn
      if (cdum_t(j).ne.'NA') then
        read(cdum_t(j),*)y(i,j)
      else
        y(i,j)=-99999   !missing
      end if
    end do
    !print*,y(i,:)
    !pause
  end do
  close(40)

  !check duplications in dat file **********************************
  do i=1,rn
    do j=i+1,rn
      if (obs(i,1)==obs(j,1).and.obs(i,2)==obs(j,2)) then
        print*,'ID duplications in .dat: ',trim(obs(i,1)),' ',trim(obs(i,2))
        pause
      end if
    end do
  end do !**********************************************************

  !check members in ped match those in data file *******************
  do i=1,rn
    do j=1,n
      if (obs(i,1)==pedor(j,1).and.obs(i,2)==pedor(j,2)) then
        goto 70
      end if
    end do
    print*,'ID in .dat is not in .ped: ',trim(obs(i,1)),' ',trim(obs(i,2))
    pause
70  continue
  end do !***********************************************************

  ! indexd id, ped, observation  ***********************************
  do i=1,n
    idx_ped(i,2)=i        !order same as in .fam or ped
    idx_ped(i,3:4)=0      !because NRM or GRM is provided
  end do
  do i=1,rn               !order same as in .dat
    do j=1,n
      if (obs(i,1)==pedor(j,1).and.obs(i,2)==pedor(j,2)) goto 11
    end do
    print*,'no matched >>> check'
    pause

11  continue
    idx_obs(i,2)=j        !order same as in .dat
    idx_y(i,:)=y(i,:)
  end do

end if !fl2.ne.'null'


if (fl12.ne."null") then

else
  if (fl3.ne.'null') then
    if (mg==0.or.tmg==0.or.bmg==0.or.zmg==0) then
        grm_known(1)=fl3
    else if (mg==1.or.bmg==1.or.zmg==1) then
      open (unit=39,file=fl3,status='old')      !fam file
      do i=1,mn
        read(39,*)grm_known(i)
      end do
      close(39)
    else if (tmg==1) then
      open (unit=39,file=fl3,status='old')      !fam file
      do i=1,mn
        read(39,*)tzb(i),grm_known(i)
      end do
      close(39)
    end if
  end if

end if


!fixed effects ***********************************************
fnm=1          !at the moment, no. fixed effects=1
cc_c=0
if (cc==1) then
  open (unit=41,file=fl6,status='old')
  read(41,'(a)')string
  close(41)

  do i=1,1000
    read(string,*,iostat=io)(cdum4,j=1,i)
    if (io.ne.0) exit
  end do
  cc_c=i-1
  cc_c=cc_c-2
  !print*,'cc_c:',cc_c
  if (cc_c>1000) then
    print*,'no. column is over 1000 >>> check',cc_c
  end if

  inquire(file=trim(fl6)//'.det',exist=file_exist)
  if (file_exist) then
    open (UNIT=43,FILE=trim(fl6)//'.det',STATUS='old')
    read(43,*)ccdet(1:cc_c)
    close(43)

    open (unit=41,file=fl6,status='old')
    do i=1,n
      read(41,*)cdum1,cdum2,(ccm(i,j),j=1,cc_c)
      if (cdum1.ne.pedor(idx_obs(i,2),1).or.cdum2.ne.pedor(idx_obs(i,2),2)) then
        print*,'cc file order should match with dat file'
        pause
      end if
      do j=1,cc_c
        if (ccm(i,j)=='NA') then
          !idx_y(i,:)=-99999     !missing
          idx_y(i,ccdet(j))=-99999     !missing
        end if
      end do
    end do
    close(41)

  else 

    open (unit=41,file=fl6,status='old')
    do i=1,n
      read(41,*)cdum1,cdum2,(ccm(i,j),j=1,cc_c)
      if (cdum1.ne.pedor(idx_obs(i,2),1).or.cdum2.ne.pedor(idx_obs(i,2),2)) then
        print*,'cc file order should match with dat file'
        pause
      end if
      do j=1,cc_c
        if (ccm(i,j)=='NA') then
          idx_y(i,:)=-99999     !missing
        end if
      end do
    end do
    close(41)

    ccdet(1:cc_c)=1  !first trait
    do xi=2,tn
      do j=1,cc_c
        ccdet(cc_c*xi-cc_c+j)=xi
        ccm(:,cc_c*xi-cc_c+j)=ccm(:,j)
      end do
    end do
    cc_c=cc_c*tn
  end if 

end if

qc_c=0
if (qc==1) then
  open (unit=42,file=fl7,status='old')
  read(42,'(a)')string
  close(42)

  do i=1,1000
    read(string,*,iostat=io)(cdum4,j=1,i)
    if (io.ne.0) exit
  end do
  qc_c=i-1
  qc_c=qc_c-2
  !print*,'qc_c:',qc_c
  if (qc_c>1000) then
    print*,'no. column is over 1000 >>> check',qc_c
  end if

  inquire(file=trim(fl7)//'.det',exist=file_exist)
  if (file_exist) then
    open (UNIT=44,FILE=trim(fl7)//'.det',STATUS='old')
    read(44,*)qcdet(1:qc_c)
    close(44)

    open (unit=42,file=fl7,status='old')
    do i=1,n
      read(42,*)cdum1,cdum2,(tmp_qcm(j),j=1,qc_c)
      if (cdum1.ne.pedor(idx_obs(i,2),1).or.cdum2.ne.pedor(idx_obs(i,2),2)) then
        print*,'qc file order should match with dat file'
        pause
      end if

      do j=1,qc_c
        if (tmp_qcm(j)=='NA') then
          !idx_y(i,:)=-99999     !missing
          idx_y(i,qcdet(j))=-99999     !missing
        else
          read(tmp_qcm(j),*)qcm(i,j)
        end if
      end do
    end do
    close(42)

  else

    open (unit=42,file=fl7,status='old')
    do i=1,n
      read(42,*)cdum1,cdum2,(tmp_qcm(j),j=1,qc_c)
      if (cdum1.ne.pedor(idx_obs(i,2),1).or.cdum2.ne.pedor(idx_obs(i,2),2)) then
        print*,'qc file order should match with dat file'
        pause
      end if

      do j=1,qc_c
        if (tmp_qcm(j)=='NA') then
          idx_y(i,:)=-99999     !missing
        else
          read(tmp_qcm(j),*)qcm(i,j)
        end if
      end do
    end do
    close(42)

    qcdet(1:qc_c)=1  !first trait
    do xi=2,tn
      do j=1,qc_c
        qcdet(qc_c*xi-qc_c+j)=xi
        qcm(:,qc_c*xi-qc_c+j)=qcm(:,j)
      end do
    end do
    qc_c=qc_c*tn

  end if

end if


if (fl2.ne.'null') then 
  !find out no. record for the traits ********************
  rnm=0
  do xi=1,tn
    do i=1,rn
      if (idx_y(i,xi).ne.-99999) then
        rnm(xi)=rnm(xi)+1
        yidx(xi,rnm(xi))=idx_obs(i,2)      !indexing pedigree
      end if
    end do
  end do
end if !fl2.ne.'null'

!print*,ccdet(1:cc_c)
!print*,qcdet(1:qc_c)

allocate(cc_zi(tn,cc_c),ccn(tn,cc_c,1000))
cc_lvl=0;cc_zi=0;ccn='null';ccdet2=0
!making xm *******************************************
!do xi=1,tn
  do j=1,cc_c
    ccdet2(ccdet(j))=ccdet2(ccdet(j))+1
    !print*,'ccdet',j,ccdet(j),ccdet2(ccdet(j))
    do i=1,n
      do k=1,cc_zi(ccdet(j),ccdet2(ccdet(j)))
        if (ccm(i,j)==ccn(ccdet(j),ccdet2(ccdet(j)),k)) goto 50
      end do
      !if (idx_y(i,xi).ne.-99999) then
      if (idx_y(i,ccdet(j)).ne.-99999) then
        k=ccdet2(ccdet(j))
        !cc_lvl(ccdet(j))=cc_lvl(ccdet(j))+1
        cc_zi(ccdet(j),k)=cc_zi(ccdet(j),k)+1
        ccn(ccdet(j),k,cc_zi(ccdet(j),k))=ccm(i,j)
        !print*,cc_zi(ccdet(j),k)
        !print*,trim(ccn(ccdet(j),k,cc_zi(ccdet(j),k)))
      end if

50    continue
    end do
  end do

  ccdet2=0
  do j=1,cc_c
    ccdet2(ccdet(j))=ccdet2(ccdet(j))+1
    k=ccdet2(ccdet(j))
    !print*,'trait',ccdet(j),'cc level',cc_zi(ccdet(j),k),':',(trim(ccn(ccdet(j),k,k2)),' ',k2=1,cc_zi(ccdet(j),k))
    fnm(ccdet(j))=fnm(ccdet(j))+(cc_zi(ccdet(j),k)-1)
  end do

  xmm=0;ccdet2=0
  do j=1,cc_c

          jj=0;ccdet3=0
          do kk=1,(j-1)
            if (ccdet(kk)==ccdet(j)) then
              ccdet3(ccdet(kk))=ccdet3(ccdet(kk))+1
              jj=jj+cc_zi(ccdet(kk),ccdet3(ccdet(kk)))-1
            end if
          end do

    ccdet2(ccdet(j))=ccdet2(ccdet(j))+1
    k2=ccdet2(ccdet(j))
    do i=1,n
      do k=1,(cc_zi(ccdet(j),k2)-1)
        if (ccm(i,j)==ccn(ccdet(j),k2,k)) then
          xmm(ccdet(j),i,1+jj+k)=1                 !from cc
          goto 60
        end if
      end do
60    continue
    end do
  end do

  do j=1,qc_c
    jj=0
    do kk=1,j
      if(qcdet(kk)==qcdet(j)) then
        jj=jj+1
      end if
    end do
    do i=1,n
      xmm(qcdet(j),i,fnm(qcdet(j))+jj)=qcm(i,j)
    end do
  end do

  do j=1,qc_c
    fnm(qcdet(j))=fnm(qcdet(j))+1
  end do

jj=0
do xi=1,tn
  do i=1,n
    xmm(xi,i,1)=1          !mean
  end do

  jj=jj+1
  if (fl2.ne.'null') print'(i4,a7,i2,a7)',jj,'trait',xi,'mean'
  ccdet2=0
  do j=1,cc_c
    ccdet2(ccdet(j))=ccdet2(ccdet(j))+1
    k=ccdet2(ccdet(j))
    if (ccdet(j)==xi) then
      !print*,'trait',xi,'cc file #',j,'eff. name:',(trim(ccn(ccdet(j),k,k2)),' ',k2=1,cc_zi(ccdet(j),k)-1)
      do k2=1,cc_zi(ccdet(j),k)-1
        jj=jj+1

        inquire(file=trim(fl6)//'.det',exist=file_exist)
        if (file_exist) then
          print'(i4,a7,i2,i5,a14,a12,a10)',jj,'trait',xi,j,'th in cc file ',' factor name: ',trim(ccn(ccdet(j),k,k2))
        else
          print'(i4,a7,i2,i5,a14,a12,a10)',jj,'trait',xi,j-cc_c/tn*(xi-1),'th in cc file ',' factor name: ',trim(ccn(ccdet(j),k,k2))
        end if
      end do
    end if
  end do

  do j=1,qc_c
    if (qcdet(j)==xi) then
      jj=jj+1
      inquire(file=trim(fl7)//'.det',exist=file_exist)
      if (file_exist) then
        print'(i4,a7,i2,i5,a14)',jj,'trait',xi,j,'th in qc file '
      else
        print'(i4,a7,i2,i5,a14)',jj,'trait',xi,j-qc_c/tn*(xi-1),'th in qc file '
      end if

    end if
  end do


end do




PRINT*,'*********************************************************************'
PRINT*,'MTGREML, MTGBLUP, SNP BLUP, Random regression and many'
PRINT*,'The length (row) of ID and data should be the same'
PRINT*,'The order of GRM follows ID file'
PRINT*,'The order of covariate file should be the same as ID file'
PRINT*,'Cite "Maier et al (2015) AJHG 96: 283-294" or' 
PRINT*,'     "Lee and van der Werf (2016) Bioinformatics 32: 1420-1422' 
PRINT*,'*********************************************************************'
print*,''

if (fl12.ne."null") then  !if eigenvector and eigenvalues to be used

  open (unit=71,file=trim(fl12)//".evec",status='old')
  do i=1,n
    read (71,*)GRM(1,i,:)   !because order same as in .fam (check in RTMX)
  end do
  close(71)

  open (unit=72,file=trim(fl12)//".eval",status='old')
  do i=1,n
    read (72,*)eval(i,1)   !because order same as in .fam (check in RTMX)
  end do
  close(72)

  if (fl_rrm.ne."null" .and. rrme.ne.0) then
    eval(:,2)=1     !for rrm for residual(co)variance
  end if
  if (fl_spl.ne."null" .and. rrme.ne.0) then
    eval(:,2)=1     !for rrm for residual(co)variance
  end if

  print*,'grm (eig) reading done **************************'
  print*,''

elseif (tmg==1) then
  do zi=1,mn
    if (tzb(zi)=="txt" .or. tzb(zi)=="text" .or. tzb(zi)=="t") then
      open (unit=71,file=grm_known(zi),status='old')
      print*,'grm file:',grm_known(zi)
      do
        read (71,*,iostat=io)i,j,x1   !order same as in .fam (check in RTMX)
        if (io.ne.0) exit
        !print*,i,j,x1
        GRM(zi,i,j)=x1
        GRM(zi,j,i)=x1
      end do
      close(71)
    elseif (tzb(zi)=="gzip" .or. tzb(zi)=="gz" .or. tzb(zi)=="z" .or. tzb(zi)=="g") then
      call open_infile(grm_known(zi), port, ioerr)
      print*,'grm file:',grm_known(zi) !,mn
      do
        call readline(port, linbuf, advance='no', ios=ioerr)
        if (ioerr /= 0 .and. ioerr /= -2) then
          exit
        end if
        read(linbuf,*)i,j,cdum1,x1
        !print*,i,j,x1
        GRM(zi,i,j)=x1
        GRM(zi,j,i)=x1
      end do
    elseif (tzb(zi)=="bin" .or. tzb(zi)=="binary" .or. tzb(zi)=="b") then
      open (unit=71,file=grm_known(zi),form='unformatted',access='direct',recl=4)
      print*,'grm file:',grm_known(zi) !,mn
      k=0
      do i=1,n
        do j=1,i
          k=k+1
          !print*,k,io
          !read (71,rec=k,iostat=io)x1   !because order same as in .fam (check
          !in RTMX)
          read (71,rec=k)r1   !because order same as in .fam (check in RTMX)
          !if (io.ne.0) exit
          !print*,i,j,r1
          GRM(zi,i,j)=r1
          GRM(zi,j,i)=r1
        end do
      end do
      close(71)
    end if
  end do !zi
  print*,'grm reading done *****************************'
  print*,''

else

  if (bmg==-1 .and. zmg==-1) then
    do zi=1,mn
      open (unit=71,file=grm_known(zi),status='old')
      print*,'grm file:',grm_known(zi)
      do
        read (71,*,iostat=io)i,j,x1   !order same as in .fam (check in RTMX)
        if (io.ne.0) exit
        !print*,i,j,x1
        GRM(zi,i,j)=x1
        GRM(zi,j,i)=x1
      end do
      close(71)
    end do !zi
    print*,'grm reading done *****************************'
    print*,''

  elseif (mg==-1 .and. zmg==-1) then 
    do zi=1,mn
      open (unit=71,file=grm_known(zi),form='unformatted',access='direct',recl=4)
      print*,'grm file:',grm_known(zi) !,mn
      k=0
      do i=1,n
        do j=1,i
          k=k+1
          !print*,k,io
          !read (71,rec=k,iostat=io)x1   !because order same as in .fam (check in RTMX)
          read (71,rec=k)r1   !because order same as in .fam (check in RTMX)
          !if (io.ne.0) exit
          !print*,i,j,r1
          GRM(zi,i,j)=r1
          GRM(zi,j,i)=r1
        end do
      end do
      close(71)
    end do !zi
 
    print*,'grm reading done *****************************'
    print*,''

  elseif (mg==-1 .and. bmg==-1) then 
    do zi=1,mn
      call open_infile(grm_known(zi), port, ioerr)
      print*,'grm file:',grm_known(zi) !,mn
      do
        call readline(port, linbuf, advance='no', ios=ioerr)
        if (ioerr /= 0 .and. ioerr /= -2) then
          exit
        end if
        read(linbuf,*)i,j,cdum1,x1
        !print*,i,j,x1
        GRM(zi,i,j)=x1
        GRM(zi,j,i)=x1
      end do
    end do !zi

    print*,'grm reading done *****************************'
    print*,''
  end if


end if

!print*,'fn:',fnm !fn1,fn2
tfn=0;trn=0
do i=1,tn
  tfn=tfn+fnm(i)
  trn=trn+rnm(i)
end do
!print*,tfn

call cpu_time (t2_cpu)

!t2b ***************************************************************
if (t2b.ne.0 .and. (mg.eq.0 .or. bmg.eq.0 .or. zmg.eq.0)) then
  call t2b_convert(GRM,n,fl3,t2b)
  goto 1100
end if

!pca ***************************************************************
if (pca.ne.0 .and. (mg.eq.0 .or. bmg.eq.0 .or. zmg.eq.0)) then
  call mtg_pca (GRM,n,pca,fl3)
  goto 1100
end if

!bend ***************************************************************
if (bend.ne.0 .and. (mg.eq.0 .or. bmg.eq.0 .or. zmg.eq.0)) then
  call mtg_bend (GRM,n,n,fl3,bend,thn)
  goto 1100
end if

!inv ***************************************************************
if (inv.ne.0 .and. (mg.eq.0 .or. bmg.eq.0 .or. zmg.eq.0)) then
  call mtg_inv (GRM,n,n,fl3)
  goto 1100
end if

!chol ***************************************************************
if (chol.ne.0 .and. (mg.eq.0 .or. bmg.eq.0 .or. zmg.eq.0)) then
  call mtg_chol (GRM,n,n,fl3)
  goto 1100
end if

!matvec **************************************************************
!if (matvec.ne.'null' .and. (tmg.eq.0 .and. mg.and.0 .and. bmg.eq.0 .and. zmg.eq.0)) then
if (matvec.ne.'null') then
  call mtg_matvec (GRM,n,n,matvec,fl3)
  goto 1100
end if

!matmat **************************************************************
!if (matmat.ne.0 .and. (tmg.ne.0 .or. mg.ne.0 .or. bmg.ne.0 .or. zmg.ne.0)) then
if (matmat.ne.0) then
  if (mn.ne.2) then
    print*,'GRM # should be 2'
    pause
  end if
  call mtg_matmat (GRM,n,n,matmat,fl3)
  goto 1100
end if

!frequency *********************************************************
if (frq.ne.0) then
  call freq_vx (plink)
  goto 1100
end if

!constructing relationship matrix ************************************
if (rtmx.ne.-999) then
  call rtmx_ (plink,rtmx,wt_par,fl5)
  goto 1100
end if

!constructing H relationship matrix ************************************
if (hrtmx.ne.'null') then
  !print*,pedn,n,mn
  call hrtmx_ (hrtmx,GRM,pedn,n,mn,pedor)
  goto 1100
end if

!simulation ************************************
if (simreal.ne.'null') then
  call sim_real_sub (plink,simreal)
  goto 1100
end if
if (simstp.ne.'null') then
  call sim_stp (plink,simstp)   !snp_blup
  goto 1100
end if

!constructing product matrix ****************************************
if (pdmx.ne.'null') then
  call pdmx_ (n,pedor,pdmx,GRM,mn,fl5)   !product matrix
  goto 1100
end if

!snp_blup **********************************************************
if (sbv.ne.'null') then
  call sblup(plink,vgpy,sbv,fl5)   !snp_blup
  goto 1100
end if

!var_rt *************************************************************
if (var_rt.ne.'null') then
  call var_rt_ (GRM,n,mn,var_rt,fl5)   !snp_blup
  goto 1100
end if




!REML **************************************************************
if (fl3.ne.'null' .or. fl12.ne."null") then
  if (icove==0) then
    if (irgt.ne.-9.0D0) then
      if (irgt==0.0D0) then
        print*,'will be updated later'
        stop
      end if
      call aireml_b_h_rg (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,irgt,nit,conv)

    else
      if (fl12.ne."null") then    !eig
        if (fl_rrm.ne."null") then   !rrm
          if (mrrm==0) then     !uni rrm
            call aireml_rrm_eig (GRM,eval,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rrm,xmm,nit,conv)
          elseif (mrrm>=1) then  !multi rrm
            call aireml_rrm_eig_m (GRM,eval,yidx,idx_y,rn,rnm,trn,n,mrrm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rrm,xmm,nit,conv)
          end if
        elseif (fl_spl.ne."null") then
          call aireml_spl_eig (GRM,eval,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_spl,xmm,nit,conv)

        else
          if (aug==0) then
            if (gwas==2) then
              print*,'with -eig, approximate GWAS not available'
              pause
            end if
            if (gwas==1) then
              if (tn.ne.1) then
                print*,'-gwas should be used with -mod 1 (univariate GWAS only)'
                pause
              end if
              tfn=tfn+1
              call aireml_m_eig_gwas (GRM,eval,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,xmm,nit,conv,plink)

            else
              call aireml_m_eig (GRM,eval,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,xmm,nit,conv)
            end if 
          else  ! missing phenotype augmentation
            call aireml_m_eig_aug (GRM,eval,yidx,idx_y,rn,rnm,(n*tn),n,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,nit,conv)
          end if
        end if
      else    !if g, not eig
        if (fl_rrm.ne."null") then   !rrm
          if (mrrm==0) then     !uni rrm
            call aireml_rrm_h (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rrm,xmm,nit,conv)
          elseif (mrrm>=1) then  !multi rrm
            call aireml_rrm_m_h (GRM,yidx,idx_y,rn,rnm,trn,n,mrrm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rrm,xmm,nit,conv)
          end if

        elseif (fl_spl.ne."null") then
          call aireml_spl_h (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_spl,xmm,nit,conv)
        elseif (fl_rnm.ne."null") then
          if (mrnm==0) then   !uni rnm
            call aireml_rnm_h (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rnm,xmm,nit,conv)
          elseif (mrnm>=1) then  !multi rnm
            call aireml_rnm_m_h (GRM,yidx,idx_y,rn,rnm,trn,n,mrnm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl_rnm,fl24,xmm,nit,conv,thn)
          end if
        elseif (fl_mcv.ne."null") then
          if (mrnm==0) then   !uni rnm
            print*,'-mcv should be used with -mrnm 1 or higher value'
            pause
          elseif (mrnm>=1) then  !multi rnm
            call aireml_mcv_m_h (GRM,yidx,idx_y,rn,rnm,trn,n,mrnm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_mcv,fl24,xmm,nit,conv,thn)
          end if

        else
          if (gwas.ne.0) then
            if (tn .ne. 1) then 
              print*,'-gwas should be used with -mod 1 (univariate GWAS only)'
              pause
            end if
            if (gwas==1) then
              tfn=tfn+tn
              call aireml_m_gwas (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res,plink)
            elseif (gwas==2) then 
              if (fl4=='null') then
                print*,'-gwas 2 should be used with -sv (approximated GWAS)'
                pause
              end if
              tfn=tfn+tn
              call aireml_m_gwas2 (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res,plink)
            end if

          else
            call aireml_m_h (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res)
          end if
        end if
      end if
    end if

  elseif (icove==1) then
    if (irgt.ne.-9.0D0) then
      if (irgt==0.0D0) then
        print*,'will be updated later'
        stop
      end if
      call aireml_b_h2_rg (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,irgt,nit,conv)

    else
      if (fl12.ne."null") then    !eig
        if (fl_rrm.ne."null") then   !rrm
          if (mrrm==0) then     !uni rrm
            print*,'-cove 1 should be with -mrrm >>> check'
            pause
          elseif (mrrm>=1) then  !multi rrm
            call aireml_rrm_eig_m2 (GRM,eval,yidx,idx_y,rn,rnm,trn,n,mrrm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rrm,xmm,nit,conv)
          end if
        else
          call aireml_m_eig2 (GRM,eval,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,xmm,nit,conv)
        end if
      else     ! g
        if (fl_rrm.ne."null") then   !rrm
          if (mrrm==0) then     !uni rrm
            print*,'-cove 1 should be with -mrrm >>> check'
            pause
          elseif (mrrm>=1) then  !multi rrm
            call aireml_rrm_m_h2 (GRM,yidx,idx_y,rn,rnm,trn,n,mrrm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_rrm,xmm,nit,conv)
          end if
        elseif (fl_rnm.ne."null") then   !rnm
          if (mrnm==0) then     !uni rrm
            print*,'-cove 1 should be with -mrnm >>> check'
            pause
          elseif (mrnm>=1) then  !multi rrm
            call aireml_rnm_m_h2 (GRM,yidx,idx_y,rn,rnm,trn,n,mrnm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl_rnm,fl24,xmm,nit,conv,thn)
          end if
        elseif (fl_mcv.ne."null") then   !rnm
          if (mrnm==0) then     !uni rrm
            print*,'-cove 1 should be with -mrnm >>> check'
            pause
          elseif (mrnm>=1) then  !multi rrm
            call aireml_mcv_m_h2 (GRM,yidx,idx_y,rn,rnm,trn,n,mrnm,tn,fnm,tfn,mn,fl4,fl5,fl11,fl_mcv,fl24,xmm,nit,conv,thn)
          end if

        else  !multivariate LMM
          call aireml_m_h2 (GRM,yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,fl11r,fl24,xmm,nit,conv,wt_res)
        end if
      end if
    end if

  end if
else
  if (icove==0) then
    call null_m_h (yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,nit,conv)
  elseif (icove==1) then
    call null_m_h2 (yidx,idx_y,rn,rnm,trn,n,tn,fnm,tfn,mn,fl4,fl5,fl11,xmm,nit,conv)
  end if
end if


1100 continue

call cpu_time (t1_cpu)
!print*,'Estimation done - time:',real(t1_cpu-t2_cpu) !,LA



end program

