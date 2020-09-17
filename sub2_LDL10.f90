
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! making NRM 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine Amat(A,ped,nped)

REAL::A(nped,nped)
REAL::inbco
INTEGER::i,z,nped,ped(nped,4),dn
A=0
A(1, 1) = 1
 do z = 2, nped
   If (ped(z, 3) .ne. 0 .And. ped(z, 4).ne. 0) Then
     do i = 1, z - 1
       A(z, i) = 0.5 * (A(ped(z, 3), i) + A(ped(z, 4), i))
       A(i, z) = A(z, i)
       A(z, z) = 1 + 0.5 * A(ped(z, 3), ped(z, 4))
     end do !Next i
   ElseIf (ped(z, 3).ne. 0 .And. ped(z, 4) == 0) Then
     do i = 1, z - 1
       A(z, i) = 0.5 * (A(ped(z, 3), i))
       A(i, z) = A(z, i)
       A(z, z) = 1
     end do !Next i
   ElseIf (ped(z, 3) == 0 .And. ped(z, 4) .ne. 0) Then
     do i = 1, z - 1
       A(z, i) = 0.5 * (A(ped(z, 4), i))
       A(i, z) = A(z, i)
       A(z, z) = 1
     end do !Next i
   ElseIf (ped(z, 3) == 0 .And. ped(z, 4) == 0) Then
     do i = 1, z - 1
       A(z, i) = 0
       A(i, z) = A(z, i)
       A(z, z) = 1
     end do !Next i
   End If
 end do !Next z

end subroutine

!************************************************************************
! Program to produce standard normal distributed (mean 0, variance 1)
! Random deviates (gaussian distribution)
! See numerical recipes page 203
! from MTBLUP (Julius van der Werf)
!************************************************************************

        REAL FUNCTION GASDEV(DSEED)
        DOUBLE PRECISION DSEED
        COMMON /GSDV/GSET,ISET
        ISET=0
1           IF (ISET.EQ.0)THEN
                V1=2.*RAN1(DSEED)-1.
                V2=2.*RAN1(DSEED)-1.
                R=V1**2+V2**2
                IF(R.GE.1.)GO TO 1
                FAC=SQRT(-2.*LOG(R)/R)
                GSET=V1*FAC
                GASDEV=V2*FAC
                ISET=1
                ELSE
                GASDEV=GSET
             ISET=0
             END IF
        RETURN
        END FUNCTION

! PROGRAM TO PRODUCE UNIFORMLY DISRIBUTED (BETWEEN 0 AND 1)
! RANDOM DEVIATES
!   SAME AS  REAL FUNCTION GGUBFS(IMSL )
      REAL FUNCTION RAN1(DSEED)
      DOUBLE PRECISION DSEED
      DOUBLE PRECISION D2P31M,D2P31
      DATA D2P31M/2147483647.D0/
      DATA D2P31 /2147483648.D0/
      DSEED=DMOD(16807.D0*DSEED,D2P31M)
      RAN1=DSEED/D2P31
      RETURN
      END FUNCTION


!end module


  subroutine minv(mat,n,imat)

  implicit none
  integer::n
  real::mat(n,n),imat(n,n)
  real::lmat(n,n),diam(n),lmat2(n,n),invm(n, n),invm2(n,n)
  integer::i,j,k,k1,k2,zi
  real::v1,x1

  !cholesky decomposition (Golub and Van Loan, 1983)
  !i=1~n, j=i+1~n
  !l(i,i)=sqrt(mat(i,i)-sum(l(i,k)^2,k=1~i-1)
  !l(j,i)=(mat(j,i)-sum(l(j,k)*l(i,k),k=1~i-1)/l(i,i)

  !do i=1,n
  !  print*,mat(i,:)
  !end do
  !pause

  do i=1,n
    x1 = 0
    do k=1,i-1
      x1=x1+lmat(i,k)**2
    end do
    lmat(i,i)=sqrt(mat(i,i)-x1)
 
    do j=i+1,n
      x1 = 0
      do k=1,i-1
        x1 = x1 + lmat(j, k) * lmat(i, k)
      end do
      lmat(j,i)=(mat(j,i)-x1)/lmat(i,i)
    end do
    !print*,lmat(i,:)
  end do
  !pause

  !now inverse lower triangle matrix

  do i=1,n
    diam(i) = 1 / lmat(i, i)
  end do

  lmat2(1, 1) = 1
  do i=2,n
    lmat2(i, i) = 1
    do j=1,i-1
      lmat2(i,j) = diam(i)*lmat(i, j)
    end do
  end do


  !inverse from lmat2
  do i = 1, n
    invm(i, i) = 1
  end do

  k1 = 1
  do zi = 1, n - 1
    k1 = k1 + 1
    do i = k1, n
      v1 = 0
      do j = i - (k1 - 2), i
        v1 = v1 + invm(i, j) * lmat2(j, i - (k1 - 1))
      end do
      invm(i, i - (k1 - 1)) = -v1
    end do
  end do


  do i = 1, n
    do j = 1, i
      invm2(i, j) = invm(i, j) * diam(j)
    end do
  end do

  !imat=matmul(transpose(invm2),invm2)

  do zi=1,n
    do i=1,zi !n
      v1=0
      do j=1,n
        v1 = v1 + invm2(j, i) * invm2(j, zi)
      end do
      imat(zi,i)=v1
      imat(i,zi)=v1
    end do
    !print*,imat(zi,:)
  end do
  !pause


  end subroutine




  subroutine cholesky_inv (mat,n,det)

  implicit none
  integer::n
  double precision::mat(n,n)!,imat(n,n)
  !double precision::lmat(n,n),diam(n),lmat2(n*(n+1)/2),invm(n, n),invm2(n,n)
  integer::i,j,k,k1,k2,zi,err
  double precision::v1,x1,det
  double precision:: t1_cpu, t2_cpu

  !cholesky decomposition (Golub and Van Loan, 1983)
  !i=1~n, j=i+1~n
  !l(i,i)=sqrt(mat(i,i)-sum(l(i,k)^2,k=1~i-1)
  !l(j,i)=(mat(j,i)-sum(l(j,k)*l(i,k),k=1~i-1)/l(i,i)

  !if (n<40) then
  !do i=1,n
  !  print*,mat(i,:)
  !end do
  !pause
  !end if

  !call cpu_time (t1_cpu)
  call dpotrf ('L',n,mat,n,err)
  !call cpu_time (t2_cpu)
  if (err.ne.0) then
    !print*,'spptrf:',err!,t2_cpu-t1_cpu
    mat(1,1)=-1
  end if

  det=0
  do j=1,n
    det=det+log(mat(j,j))
  end do
  det=det*2   !because L'L
  !print*,det

  !now inverse lower triangle matrix using lapack subroutines

  !call cpu_time (t1_cpu)
  call dpotri ('L',n,mat,n,err)
  !call cpu_time (t2_cpu)
  !print*,'spptri:',err,t2_cpu-t1_cpu
  if (err.ne.0) then
    !print*,err
  end if

  do i=1,n
    do j=1,i
      !imat(i,j)=mat(i,j)
      !imat(j,i)=mat(i,j)
      mat(j,i)=mat(i,j)
    end do
  end do

  end subroutine



  !subroutine lu_minv (mat,n,imat,det)
  subroutine lu_minv (mat,n,det)

  implicit none
  integer::n
  double precision::mat(n,n)!,imat(n,n)
  double precision::work(n*n)
  integer::i,j,k,k1,k2,zi,err,ipiv(n)
  double precision::v1,x1,det
  double precision:: t1_cpu, t2_cpu

  call dgetrf (n,n,mat,n,ipiv,err)

  det=0
  do j=1,n
    det=det+log(mat(j,j))
  end do
  !print*,det

  call dgetri (n,mat,n,ipiv,work,n*n,err)
  print*,err

  !imat=mat


  end subroutine




  subroutine msol2(mat,n,bmat,fn,det)

  implicit none
  integer::n,fn
  real::mat(n,n),imat(n,n),bmat(n,fn)
  real::lmat(n,n),diam(n),lmat2(n*(n+1)/2),invm(n, n),invm2(n,n)
  integer::i,j,k,k1,k2,zi,err
  real::v1,x1,det
  double precision:: t1_cpu, t2_cpu

  !cholesky decomposition (Golub and Van Loan, 1983)
  !i=1~n, j=i+1~n
  !l(i,i)=sqrt(mat(i,i)-sum(l(i,k)^2,k=1~i-1)
  !l(j,i)=(mat(j,i)-sum(l(j,k)*l(i,k),k=1~i-1)/l(i,i)

  !do i=1,n
  !  print*,mat(i,:)
  !end do
  !pause

  k=1
  do j=1,n
    do i=j,n
      lmat2(k+i-j)=mat(i,j)
    end do
    k=k+n-j+1
  end do
  
  call cpu_time (t1_cpu)
  call spptrf ('L',n,lmat2,err)
  call cpu_time (t2_cpu)
  print*,'spptrf:',err,t2_cpu-t1_cpu

  k=1;det=0
  do j=1,n
    det=det+log(lmat2(k))
    k=k+n-j+1
  end do
  det=det*2   !because L'L

  !now solve linear system
  !X* = V-1 X, X* = (LL')-1 X, Ly = X, L'X* = y

  call cpu_time (t1_cpu)
  call spptrs ('L',n,fn,lmat2,bmat,n,err)
  call cpu_time (t2_cpu)
  print*,'spptrs:',err,t2_cpu-t1_cpu

  !imat=0
  !k=1
  !do j=1,n
  !  do i=j,n
  !    imat(i,j)=lmat2(k+i-j)
  !    imat(j,i)=lmat2(k+i-j)
  !  end do
  !  k=k+n-j+1
  !end do



  end subroutine


FUNCTION FLEGENDRE(N,X) 

  REAL*8 PI,PIM1,PIM2 
  REAL*8 FI 
  INTEGER I 
  ! CHECK FOR VALID VALUES OF N AND X HERE BEFORE PROCEEDING C C PROCEEDING WITH CALCULATIONS

  IF (N.EQ.0) THEN 
    FLEGENDRE=1 
  ELSEIF (N.EQ.1) THEN 
    FLEGENDRE=X 
  ELSE 
    ! SYNTHETIC CALCULATIONS 
    PIM1=1 
    PI=X 
    DO 20 I=2,N 
    FI=I 
    PIM2=PIM1 
    PIM1=PI
    PI=((I+I-1)*X*PIM1-(I-1)*PIM2)/FI 
    20 CONTINUE 
    FLEGENDRE=PI 
  ENDIF 


END

 subroutine spline (wk1,kn,wk2,tn,phi)
  implicit none

  integer::kn,tn
  double precision::wk1(kn),wk2(tn),phi(tn,kn),h(kn-1),Zs(kn,kn-2)
  double precision::delta(kn,(kn-2)),Gs(kn-2,kn-2),d2(kn-2,kn-2),Zv(kn,kn-2)
  integer::i,j,k,err,zi
  !real::det
  double precision::det
  real::XI(kn),FI(kn),P2(kn)
  REAL :: X, F, DX, ALPHA, BETA, GAMMA, ETA

  !print*,'spline'
  !Verbyla et al. (1999) Appl. Statist. p 278.
  do i=1,kn-1
    h(i)=wk1(i+1)-wk1(i)
  end do

  delta=0
  do i=1,kn-2
    delta(i,i)=1/h(i)
    delta(i+1,i)=-(1/h(i)+1/h(i+1))
    delta(i+2,i)=1/h(i+1)
  end do

  Gs=0
  do i=1,kn-2
    Gs(i,i)=(h(i)+h(i+1))/3
    if (i<kn-2) then
      Gs(i,i+1)=h(i+1)/6
      Gs(i+1,i)=h(i+1)/6
    end if
  end do

  d2=matmul(transpose(delta),delta)
  call cholesky_inv (d2,kn-2,det)
  !print*,det
  Zs=matmul(delta,d2)

  call dpotrf ('U',kn-2,Gs,kn-2,err)
  !print*,err
  do i=2,kn-2
    do j=1,i-1
      Gs(i,j)=0
    end do
  end do

  Zv=matmul(Zs,transpose(Gs))

  !Zv=Zv/1.171219 !Zv*0.854 !0.8538      ! check make it the same as asreml
  !do i=1,kn
  !  print*,Zv(i,:)
  !  print*,Zs(i,:)
  !end do
  !print*,''

  do zi=1,kn-2

    XI=wk1
    FI=Zv(:,zi)
    call cubic_spline(kn-1,XI,FI,P2) !from splin3.f90 in /randR/ (Tao Pang, 2006)

    X = XI(1)
    !DO I = 1, tn-1
    DO I = 1, tn

      if (i==tn .and. wk2(i)==wk1(kn)) then
        phi(tn,1)=1
        phi(tn,2)=wk2(i)
        phi(tn,2+zi)=Zv(kn,zi)  !check
        exit
      end if

      X = wk2(i)
      ! Find the interval that x resides
      K = 1
      DX = X-XI(1)
      DO WHILE (DX .GE. 0)
        K = K + 1
        DX = X-XI(K)
      END DO
      K = K - 1

      ! Find the value of function f(x)
      DX = XI(K+1) - XI(K)
      ALPHA = P2(K+1)/(6*DX)
      BETA = -P2(K)/(6*DX)
      GAMMA = FI(K+1)/DX - DX*P2(K+1)/6
      ETA = DX*P2(K)/6 - FI(K)/DX
      F = ALPHA*(X-XI(K))*(X-XI(K))*(X-XI(K)) &
         +BETA*(X-XI(K+1))*(X-XI(K+1))*(X-XI(K+1)) &
         +GAMMA*(X-XI(K))+ETA*(X-XI(K+1))
      !WRITE (6, *) X, F
      phi(i,2+zi)=F    !check
      !phi(i,2+zi)=F/1.171219    !check, make same as asreml
      !phi(i,2+zi)=F*3    !check, make same as asreml

      phi(i,1)=1
      phi(i,2)=wk2(i)
    END DO
    !phi(tn,1)=1
    !phi(tn,2)=wk2(i)
    !phi(tn,2+zi)=Zv(kn,zi)  !check
    !phi(tn,2+zi)=Zv(kn,zi)/1.171219  !check, same as asreml
    !phi(tn,2+zi)=Zv(kn,zi)*3  !check, same as asreml

  end do ! zi

  !phi=0
  !do i=1,tn
    !do j=1,(kn-1)
    !  if (wk2(i).ge.wk1(j) .and. wk2(i).le.wk1(j+1)) then
    !    do k=1,kn-2
    !      phi(i,k+2)=Zv(j,k)+(wk2(i)-wk1(j))*(Zv(j+1,k)-Zv(j,k))/(wk1(j+1)-wk1(j))
    !      phi(i,k+2)=phi(i,k+2)*0.8538
    !      !print*,phi(i,:)
    !    end do
    !  end if
    !end do

    !print*,phi(i,:)
    !phi(i,1)=1
    !phi(i,2)=wk2(i)
  !end do


  !print*,'spline finish'
  end subroutine

SUBROUTINE CUBIC_SPLINE (N, XI, FI, P2)
!
! Function to carry out the cubic-spline approximation
! with the second-order derivatives returned.
!
  INTEGER :: I
  INTEGER, INTENT (IN) :: N
  REAL, INTENT (IN), DIMENSION (N+1):: XI, FI
  REAL, INTENT (OUT), DIMENSION (N+1):: P2
  REAL, DIMENSION (N):: G, H
  REAL, DIMENSION (N-1):: D, B, C
!
! Assign the intervals and function differences
!
  DO I = 1, N
    H(I) = XI(I+1) - XI(I)
    G(I) = FI(I+1) - FI(I)
  END DO
!
! Evaluate the coefficient matrix elements
  DO I = 1, N-1
    D(I) = 2*(H(I+1)+H(I))
    B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
    C(I) = H(I+1)
  END DO
!
! Obtain the second-order derivatives
!
  CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)
  P2(1) = 0
  P2(N+1) = 0
  DO I = 2, N
    P2(I) = G(I-1)
  END DO
END SUBROUTINE CUBIC_SPLINE

SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
!
! Functione to solve the tridiagonal linear equation set.
!
  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  REAL, INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL, INTENT (OUT), DIMENSION (L):: Z
  REAL, DIMENSION (L):: Y, W
  !REAL, DIMENSION (L-1):: V, T
  REAL, DIMENSION (L):: V, T    !check if this is ok
!
! Evaluate the elements in the LU decomposition
!
  !print*,L
  !pause
  W(1) = D(1)
  V(1)  = C(1)
  T(1)  = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)
!
! Forward substitution to obtain y
!
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO
!
! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO
END SUBROUTINE TRIDIAGONAL_LINEAR_EQ


subroutine legendre (kn,wk,tn,phi)
  integer::kn,tn
  double precision::wk(tn),phi(tn,kn),xk(tn)
  double precision::am2(20,20),mm2(tn,20)
  double precision::am(kn,kn),mm(tn,kn)

  if (kn > 20) then
    print*,"# order is more than 20, which is the maximum at the moment"
    pause
  end if 

  !Legendre polynomial function ***********************
  am2=0

am2(1,1)=0.7071068

am2(2,2)=1.224745

am2(1,3)=-0.7905694
am2(3,3)=2.371708

am2(2,4)=-2.806243
am2(4,4)=4.677072

am2(1,5)=0.7954951
am2(3,5)=- 7.954951
am2(5,5)=9.280777

am2(2,6)=4.397265
am2(4,6)=- 20.52057
am2(6,6)=18.46851

am2(1,7)=-0.7967218
am2(3,7)=16.731162
am2(5,7)=- 50.19347
am2(7,7)= 36.80855

am2(2,8)=-5.990715
am2(4,8)= 53.91644
am2(6,8)=- 118.6162
am2(8,8)=+ 73.42906

am2(1,9)=0.7972005
am2(3,9)=- 28.69922
am2(5,9)=+ 157.8457
am2(7,9)=- 273.5992
am2(9,9)=+ 146.571

am2(2,10)=7.585119
am2(4,10)=- 111.2484
am2(6,10)=+ 433.8688
am2(8,10)=- 619.8126
am2(10,10)=+ 292.6893

am2(1,11)=-0.7974349
am2(3,11)=+ 43.85892
am2(5,11)=- 380.1106
am2(7,11)=+ 1140.332
am2(9,11)=- 1384.689
am2(11,11)=+584.6464


am2(2,12)=-9.17999
am2(4,12)=+ 198.8998
am2(6,12)=- 1193.399
am2(8,12)=+ 2898.254
am2(10,12)=- 3059.268
am2(12,12)=+ 1168.084

am2(1,13)=0.7975667
am2(3,13)=- 62.2102
am2(5,13)=+ 777.6276
am2(7,13)=- 3525.245
am2(9,13)=+ 7176.392
am2(11,13)=- 6697.965
am2(13,13)=+ 2334.139

am2(2,14)=10.77512
am2(4,14)=- 323.2537
am2(6,14)=+ 2747.657
am2(8,14)=- 9943.9
am2(10,14)=+ 17401.82
am2(12,14)=- 14554.25
am2(14,14)=+ 4664.825

am2(1,15)=-0.7976481
am2(3,15)=+ 83.75305
am2(5,15)=- 1423.802
am2(7,15)=+ 9017.412
am2(9,15)=- 27052.24
am2(11,15)=+ 41480.09
am2(13,15)=- 31424.31
am2(15,15)=+ 9323.698

am2(2,16)=-12.37042
am2(4,16)=+ 490.6933
am2(6,16)=- 5593.904
am2(8,16)=+ 27969.52
am2(10,16)=- 71477.66
am2(12,16)=+ 97469.54
am2(14,16)=- 67478.91
am2(16,16)=+ 18637.03

am2(1,17)=0.7977018
am2(3,17)=- 108.4874
am2(5,17)=+ 2404.805
am2(7,17)=- 20200.36
am2(9,17)=+ 82965.78
am2(11,17)=- 184368.4
am2(13,17)=+ 226270.3
am2(15,17)=- 144216.2
am2(17,17)=+ 37255.86


am2(2,18)=13.96582
am2(4,18)=- 707.6017
am2(6,18)=+ 10401.75
am2(8,18)=- 68354.33
am2(10,18)=+ 237341.4
am2(12,18)=- 466052.2
am2(14,18)=+ 519827.5
am2(16,18)=- 306945.8
am2(18,18)=+ 74479.49

am2(1,19)=-0.7977391
am2(3,19)=+ 136.4134
am2(5,19)=- 3819.575
am2(7,19)=+ 40996.77
am2(9,19)=- 219625.6
am2(11,19)=+ 658876.7
am2(13,19)=- 1158026
am2(15,19)=+ 1183477
am2(17,19)=- 650912.2
am2(19,19)=+ 148901.5

am2(2,20)=-15.5613
am2(4,20)=+ 980.362
am2(6,20)=- 18038.66
am2(8,20)=+ 150322.2
am2(10,20)=- 676449.8
am2(12,20)=+ 1783368
am2(14,20)=- 2835097
am2(16,20)=+ 2673092
am2(18,20)=- 1375856
am2(20,20)=+ 297699.8

  wk=wk-wk(1)
  xk=(wk-wk(tn)/2)/(wk(tn)/2)

  mm2=1
  mm2(:,2)=xk
  mm2(:,3)=xk**2
  mm2(:,4)=xk**3
  mm2(:,5)=xk**4
  mm2(:,6)=xk**5
  mm2(:,7)=xk**6
  mm2(:,8)=xk**7
  mm2(:,9)=xk**8
  mm2(:,10)=xk**9
  mm2(:,11)=xk**10
  mm2(:,12)=xk**11
  mm2(:,13)=xk**12
  mm2(:,14)=xk**13
  mm2(:,15)=xk**14
  mm2(:,16)=xk**15
  mm2(:,17)=xk**16
  mm2(:,18)=xk**17
  mm2(:,19)=xk**18
  mm2(:,20)=xk**19


    am=am2(1:kn,1:kn)
    mm=mm2(1:tn,1:kn)

    phi=matmul(mm,am)

end subroutine


subroutine pol (kn,wk,tn,phi)
  integer::kn,tn
  double precision::wk(tn),phi(tn,kn),xk(tn)
  double precision::mm2(tn,20)
  double precision::mm(tn,kn)

  if (kn > 20) then
    print*,"# order is more than 20, which is the maximum at the moment"
    pause
  end if

  xk=wk 

  mm2=1
  mm2(:,2)=xk
  mm2(:,3)=xk**2
  mm2(:,4)=xk**3
  mm2(:,5)=xk**4
  mm2(:,6)=xk**5
  mm2(:,7)=xk**6
  mm2(:,8)=xk**7
  mm2(:,9)=xk**8
  mm2(:,10)=xk**9
  mm2(:,11)=xk**10
  mm2(:,12)=xk**11
  mm2(:,13)=xk**12
  mm2(:,14)=xk**13
  mm2(:,15)=xk**14
  mm2(:,16)=xk**15
  mm2(:,17)=xk**16
  mm2(:,18)=xk**17
  mm2(:,19)=xk**18
  mm2(:,20)=xk**19

  phi=mm2(1:tn,1:kn)

end subroutine


subroutine sort (v1,dn)

integer::dn
!dimension::v1(dn,2),v2(dn,2)
integer::v1(dn,2),v2(dn,2)
integer::xi,k,i,j

    v2=v1
    do xi=1,dn-1
      do i=xi,dn
        A2=v1(xi,2)
        A1=v1(xi,1)
        !if (v1(xi,2)<v1(i,2)) then
        if (v1(xi,2)>v1(i,2)) then   !asscending
          v2(xi,2)=v1(i,2)
          v2(xi,1)=v1(i,1)
          do j=xi+1,i-1
            v2(j+1,2)=v1(j,2)
            v2(j+1,1)=v1(j,1)
          end do
          v2(xi+1,2)=A2
          v2(xi+1,1)=A1
          do k=1,i
            v1(k,2)=v2(k,2)
            v1(k,1)=v2(k,1)
          end do
        end if
      end do
    end do

end subroutine


subroutine sort2 (v1,v2,dn)

integer::dn
integer::v1(dn),v11(dn)
double precision::v2(dn),v12(dn)
integer::xi,k,i,j

    v11=v1;v12=v2
    do xi=1,dn-1
      do i=xi,dn
        A2=v2(xi)
        A1=v1(xi)
        !if (v1(xi,2)<v1(i,2)) then
        if (v2(xi)>v2(i)) then   !asscending
          v12(xi)=v2(i)
          v11(xi)=v1(i)
          do j=xi+1,i-1
            v12(j+1)=v2(j)
            v11(j+1)=v1(j)
          end do
          v12(xi+1)=A2
          v11(xi+1)=A1
          do k=1,i
            v2(k)=v12(k)
            v1(k)=v11(k)
          end do
        end if
      end do
    end do

end subroutine


