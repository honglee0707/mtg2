subroutine trf_rg (t1_h2,t1_k,t1_p,t2_h2,t2_k,t2_p,rg,se)

!**********************************************************************
! Expected prediction accuracy (also see ccgpa5.r in /R/Daetwyler/)
! S. Hong Lee  (2016)
!************************************************************************

implicit none

  real::t1_h2,se,t1_k,t1_p,rg
  real::t2_h2,t2_k,t2_p,cv,rdum
  character(len=64)::indc

   print*,"*************************************************************"
   print*,"t1_    : first trait"
   print*,"t2_    : second trait"
   print*,"h2     : estimated h2 on the observed scale using LMM"
   print*,"k      : population prevalence"
   print*,"p      : proportion of cases in training sample"
   print*,"rg     : genetic correlation"
   print*,"se     : standard error of rg"
   print*,"method : bivariate LMM (Bioinformatics 28: 2540_2542)"
   print*,"*************************************************************"
   print*,""

   rdum=0
   indc='obs'
   print*,'first trait'
   call trf_h2(t1_h2,rdum,t1_k,t1_p,indc)
   print*,'second trait'
   call trf_h2(t2_h2,rdum,t2_k,t2_p,indc)

   cv=sqrt(t1_h2*t2_h2)

  rg=rg*cv
  se=se*cv
  print*,'SNP-coheritability'
  print*,'SNP-coheritability (genet covar on the liab scale): ',rg
  print*,'SE of SNP-coheritability                          : ',se
  print*,""


end subroutine


