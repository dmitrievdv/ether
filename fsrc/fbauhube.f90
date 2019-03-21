!-------------------------- Module fbauhube.f90 ----------------------------
!* Ref.: "Numerical algorithms with C, by Gisela Engeln-Muellges and       *
!*        Frank Uhlig, Springer-Verlag, 1996" [BIBLI 11].                  *
!*                                                                         *
!*                                       F90 Release By J-P Moreau, Paris. *
!*                                              (www.jpmoreau.fr)          *
!---------------------------------------------------------------------------
MODULE FBAUHUBE

  private ! Modification added March 21 2019

  Integer, Parameter :: ITERMAX = 100         ! Maximal number of function
                                              ! evaluations per root
  Real*8, Parameter :: MACH_EPS = 1.D-15      ! Small real number
  Real*8, Parameter :: EPSROOT = 31.622D6     ! SQRT(MACH_EPS) 
  Real*8, Parameter :: EPS = 64.0 * MACH_EPS  ! Accuracy for functional value
  Real*8, Parameter :: BETA = 8.0 * EPS
  Real*8, Parameter :: QR = 0.1D0             ! Real and imaginary parts of the
  Real*8, Parameter :: QI = 0.9D0             ! spiralization constant

  Integer, Parameter :: NMAX = 25

  public :: bauhub, NMAX ! Modification added March 21 2019

CONTAINS

 Function comabs    &   ! Complex absolute value ....................
              (     &                                                 
               ar,  &             ! Real part .......................
               ai   &             ! Imaginary part ..................
              )
 !*====================================================================*
 !*                                                                    *
 !*  Complex absolute value of   a                                     *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      ar,ai    REAL*8  ar, ai;                                      *
 !*               Real, imaginary parts of  a                          *
 !*                                                                    *
 !*   Return value :                                                   *
 !*   =============                                                    *
 !*      Absolute value of a                                           *
 !*                                                                    *
 !*   Functions used :   DSQRT, DABS                                   *
 !*   ==============                                                   *
 !*                                                                    *
 !*====================================================================*
   IMPLICIT REAL*8 (A-H,O-Z) 
   if (ar == 0.d0.and.ai == 0.d0) then
     comabs=0.d0;
     return
   end if

   ar1 = DABS (ar)
   ai1 = DABS (ai)

   if (ai1 > ar1) then                             !Switch  ai1 and ar1
     t=ai1; ar1=ai1; ai1=t
   end if

   if (ai1==0.d0) then
     comabs=ar1
   else
     comabs=ar1 * DSQRT(1.d0 + (ai1/ar1)**2)
   end if
   return
 End Function comabs


  Function icomdiv &  ! Complex division .............................
           (             &
            ar,          &         ! Real part of numerator ..........
            ai,          &         ! Imaginary part of numerator .....
            br,          &         ! Real part of denominator ........
            bi,          &         ! Imaginary part of denominator ...
            cr,          &         ! Real part of quotient ...........
            ci           &         ! Imaginary part of quotient ......
           )
 !*====================================================================*
 !*                                                                    *
 !*  Complex division  c = a / b                                       *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      ar,ai    REAL*8  ar, ai;                                      *
 !*               Real, imaginary parts of numerator                   *
 !*      br,bi    REAL*8  br, bi;                                      *
 !*               Real, imaginary parts of denominator                 *
 !*                                                                    *
 !*   Output parameters:                                               *
 !*   ==================                                               *
 !*      cr,ci    REAL*8  cr, ci;                                      *
 !*               Real , imaginary parts of the quotient               *
 !*                                                                    *
 !*   Return value :                                                   *
 !*   =============                                                    *
 !*      = 0      ok                                                   *
 !*      = 1      division by 0                                        *
 !*                                                                    *
 !*====================================================================*
  IMPLICIT REAL*8(A-H,O-Z) 
  if (br == 0.d0.and.bi == 0.d0) then
    icomdiv=1
    return
  end if

  if (DABS (br) > DABS (bi)) then
    tmp  = bi / br
    br   = tmp * bi + br
    cr   = (ar + tmp * ai) / br
    ci   = (ai - tmp * ar) / br
  else
    tmp  = br / bi
    bi   = tmp * br + bi
    cr   = (tmp * ar + ai) / bi
    ci   = (tmp * ai - ar) / bi
  end if

  icomdiv=0
  return
End Function icomdiv



  Integer Function bauhub (  &  !Bauhuber's method for complex polynomials
            real0,      &          ! Are the coefficients real ? .....
            scale,      &          ! Scaling ? .......................
            n,          &          ! degree of polynomial ............
            ar,         &          ! Real parts of coefficients ......
            ai,         &          ! Imaginary parts, coefficients ...
            rootr,      &          ! Real parts of roots .............
            rooti,      &          ! Imaginary parts of roots ........
            absf        &          ! Absolute value of function values
           )
 !*====================================================================*
 !*                                                                    *
 !*  bauhub uses Bauhuber's Method to find all real or complex roots   *
 !*  of a polynomial of degree n :                                     *
 !*                                               n-1           n      *
 !*      P(x) = a(0) + a(1) * x + ... + a(n-1) * x    + a(n) * x ,     *
 !*                                                                    *
 !*  where a(i), i=0, ..., n, can be complex.                          *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Application:                                                     *
 !*   ===========                                                      *
 !*      Find roots of arbitrary polynomials with complex coefficients.*
 !*      If the polynomial roots are ill-condi=ditioned, i.e., if small*
 !*      changes in the coefficients lead to large changes in the roots*
 !*      the polynomial should not be scaled. Otherwise scaling helps  *
 !*      with stability and performance.                               *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      real0    integer;                                             *
 !*        = 0    Polynomial coefficients are complex                  *
 !*       <> 0    Polynomial coefficients are real                     *
 !*      scale    integer;                                             *
 !*        = 0    no scaling                                           *
 !*       <> 0    Scaling, see procedure scpoly()                      *
 !*      n        integer;                                             *
 !*               degree of the polynomial (must be >= 1)              *
 !*      ar, ai   REAL*8;                                              *
 !*               Real and imaginary parts of the polynomial           *
 !*               coefficients (ar(0), ..., ar(n))                     *
 !*                                                                    *
 !*   Output parameters:                                               *
 !*   =================                                                *
 !*      rootr    REAL*8;     (Vector of length  n+1 )                 *
 !*               rootr(0),..,rootr(n-1) are the real parts of the n   *
 !*               roots                                                *
 !*      rooti    REAL*8;     (Vector of length n+1 )                  *
 !*               rooti(0),..,rooti(n-1) are the imaginary parts       *
 !*               of the roots                                         *
 !*      absf     REAL*8;                                              *
 !*               absf(0),..,absf(n-1) are the magnitudes of the       *
 !*               polynomial values at the computed roots              *
 !*                                                                    *
 !*   Return value:                                                    *
 !*   ============                                                     *
 !*      = 0      all is ok                                            *
 !*      = 1      n < 1 or invalid input parameter                     *
 !*      = 2      ar(n) = 0.0 and ai(n) = 0.0                          *
 !*      = 3      Iteration maximum ITERMAX exceeded                   *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Functions or Procedures used:                                    *
 !*   ============================                                     *
 !*     ibauroot():  determines one root of the polynomial             *
 !*     scpoly():    Scales the polynomial                             *
 !*     chorner():   Evaluates the polynomial                          *
 !*     polydiv():   factors off a root                                *
 !*     comabs():    magnitude of a complex number                     *
 !*                                                                    *
 !*====================================================================*
  IMPLICIT REAL*8 (A-H,O-Z) 
  integer real0, scale
  Dimension ar(0:n),ai(0:n),rootr(0:n),rooti(0:n),absf(0:n) 
  Integer res

  scalefak = 1.d0 
  if (n < 1) then     !n is to small!
    bauhub=1
    return
  end if

  if (ar(n) == 0.d0.and.ai(n) == 0.d0) then  ! Leading coefficient must
    bauhub=2                                 ! differ from zero
    return
  end if                              

  do i = 0, n                        ! store given coefficients in root
    rootr(i) = ar(i)
    rooti(i) = ai(i)
    if (i < n)  absf(i) = 0.d0
  end do

  scalefak = 1.d0
  if (scale.ne.0) then                   ! Scale polynomial, if desired
    call scpoly (n, rootr, rooti, scalefak)
  end if

  x0r = 0.d0                                           ! Starting value
  x0i = 0.d0

  do i = 0, n-1
                                                 ! compute the ith root
    res = ibauroot (n, i, rootr, rooti, x0r, x0i)

    rootr(i) = scalefak * x0r                              ! store root
    rooti(i) = scalefak * x0i

    if (res.ne.0) then
      bauhub = 3                           ! Iteration maximum reached?
      return
    end if

    ! Polynomial value of input polynomial at (rootr(i), rooti(i))

    call chorner (n, 0, ar, ai, rootr(i), rooti(i),tempr, tempi, t1, t2, t3, t4, t5)

    absf(i) = comabs(tempr, tempi)                       ! store error

    call polydiv (n, i, rootr, rooti, x0r, x0i)        ! reduce degree

    if (real0.ne.0) then                          ! New starting value
      x0i = -x0i                               ! depending on real x..
    else
      x0r = 0.d0
      x0i = 0.d0
    end if
  end do

  bauhub = 0                                          ! normal exit
  return
End Function bauhub


  Subroutine scpoly (n,       &      ! length of vector .............
                    ar,       &      ! Real part of the vector ......
                    ai,       &      ! Imaginary part of the vector .
                    scal)            ! Scaling factor ...............
 !*====================================================================*
 !*                                                                    *
 !*  scalpoly scales the polynomial P :                                *
 !*                                               n-1           n      *
 !*      P(x) = a(0) + a(1) * x + ... + a(n-1) * x    + a(n) * x ,     *
 !*                                                                    *
 !*  where all a(i), i=0, ..., n, can be complex.                      *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Eingabeparameter:                                                *
 !*   ================                                                 *
 !*      n        integer;                                             *
 !*               degree of the polynomila (>=1)                       *
 !*      ar, ai   REAL*8;                                              *
 !*               Real and imaginary parts of the coefficients         *
 !*               a(0),..,a(n)                                         *
 !*                                                                    *
 !*   Output parameters:                                               *
 !*   =================                                                *
 !*      ar, ai   REAL*8;                                              *
 !*               Real and imaginary parts of the coefficients         *
 !*               a(0),..,a(n) of the scaled polynomial.               *
 !*      scal     REAL*8                                               *
 !*               Scaling factor                                       *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Functions used:                                                  *
 !*   ==============                                                   *
 !*      comabs:  Modulus of a complex number                          *
 !*      Max   :  Maximum of two real numbers                          *
 !*                                                                    *
 !*====================================================================*
  IMPLICIT REAL*8 (A-H,O-Z) 
  Dimension ar(0:n),ai(0:n)

  scal = 0.d0
                              ! scal =
  p = comabs (ar(n), ai(n))   !              a(i)  1/(n-i)
  do i = 0, n-1               !   max{ cabs( ---- )       , i=0..n-1
                              !              a(n)
    if (ar(i).ne.0.d0.or.ai(i).ne.0.d0) then

      ai(i) = ai(i) / p
      ar(i) = ar(i) / p

      pot = comabs(ar(i),ai(i))**(1.d0/(n-i))
      scal = Max(scal, pot)
    end if
  end do

  ar(n) = ar(n)/p                     ! Absolute value of a(n) = 1
  ai(n) = ai(n)/p

  if (scal == 0.d0)  scal = 1.d0

  p=1.d0
  do i = n-1, 0, -1
    p = p*scal                !                    n-i
    ar(i)=ar(i)/p             ! a(i) = a(i) / (scal    ), i=0..n-1
    ai(i)=ai(i)/p             !
  end do
  return
End Subroutine scpoly


  Subroutine chorner (n,      &      ! highest degree in polynomial .
                     iu,      &      ! lowest degree in polynomial ..
                     ar,      &      ! Real parts of coefficients ...
                     ai,      &      ! Imaginary parts, coefficients 
                     xr,      &      ! Real part of x ...............
                     xi,      &      ! Imaginary part of x ..........
                     pr,      &      ! Real part of function value ..
                     pi,      &      ! Imaginary part of function v. 
                     p1r,     &      ! Real part of first derivative 
                     p1i,     &      ! Imaginary part, first deriv. .
                     p2r,     &      ! Real part, second derivative .
                     p2i,     &      ! Imaginary part, second deriv. 
                     rf1)            ! Error estimate for 1st deriv. 
 !*====================================================================*
 !*                                                                    *
 !*  Horner scheme for polynomial with complex coefficients.           *
 !*  We compute :                                                      *
 !*    1. Polynomial value of the polynomial P (complex) of degree     *
 !*       n - iu,                                                      *
 !*    2. value of first derivative at x,                              *
 !*    3. value of 2nd derivative at x,                                *
 !*    4. an error estimate for the first derivative.                  *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      n        int n;                                               *
 !*               Maximal degree of the polynomial ( >= 1 )            *
 !*      ar, ai   REAL*8  ar(), ai();                                  *
 !*               Real and imaginary parts of the coefficients of the  *
 !*               polynomial with the coefficients a(iu), ..., a(n)    *
 !*      x0r,x0i  REAL*8  x0r, x0i;                                    *
 !*               Real and imaginary parts of the point of evaluation  *
 !*                                                                    *
 !*   Ausgabeparameter:                                                *
 !*   ================                                                 *
 !*      pr, pi   REAL*8  pr, pi;                                      *
 !*               Real and imaginary part of the polynomial            *
 !*      p1r, p1i REAL*8  p1r, p1i;                                    *
 !*               Real and imaginary parts of the 1st derivative there *
 !*      p2r, p2i REAL*8  p2r, p2i;                                    *
 !*               Real and imaginary parts of the 2nd derivative       *
 !*      rf1      REAL*8  rf1;                                         *
 !*               Error estimate for the first derivative              *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Functions used:                                                  *
 !*   ==============                                                   *
 !*     comabs():    modulus of a complex number                       *
 !*====================================================================*
 !*                                                                    *
 !*   Constant used:  EPS                                              *
 !*   =============                                                    *
 !*                                                                    *
 !*====================================================================*
  IMPLICIT REAL*8 (A-H,O-Z) 
  Dimension ar(0:n), ai(0:n)

  p2r = ar(n)
  p2i = ai(n)

  pr = p2r
  p1r = p2r
  pi = p2i
  p1i = p2i

  rf1 = comabs(pr, pi)
  i1 = n - iu

  j=n-iu
  do i = n - 1, iu, -1
    temp = pr                         ! Polynomial value (pr,pi)
    pr = pr * xr - pi * xi + ar(i)
    pi = pi * xr + temp * xi + ai(i)
    if (i == iu)  goto 20                          ! exit i loop

    temp = p1r                        ! 1st derivative (p1r,p1i)
    p1r = p1r * xr - p1i * xi
    p1i = p1i * xr + temp * xi

    temp = comabs(p1r, p1i)           ! Error estimate for the 1st
    p1r = p1r + pr                    ! derivative of P
    p1i = p1i + pi

    temp1 = comabs (pr, pi)
    temp = Max(temp, temp1)

    if (temp > rf1) then
      rf1 = temp
      i1 = j - 1
    end if

    if (i - iu <= 1) goto 10                    !go to next i, j

    temp = p2r                           ! 2nd derivative (p2r,p2i)
    p2r = p2r * xr - p2i * xi + p1r
    p2i = p2i * xr + temp * xi + p1i
10  j=j-1
  end do

20 temp = comabs(xr, xi)

  if (temp.ne.0.d0) then
    rf1 = rf1 * temp**i1 * (1.d0*i1 + 1.d0)
  else
    rf1 = comabs(p1r, p1i)
  end if

  rf1 = rf1 * EPS

  p2r = p2r + p2r
  p2i = p2i + p2i

  return
End Subroutine chorner

Subroutine quadsolv   &  ! Complex quadratic equation ................
             (        &
               ar,    &             !second degree coefficient .......
               ai,    &
               br,    &             !linear coefficient ..............
               bi,    &
               cr,    &             !polynomial constant .............
               ci,    &
               tr,    &             !solution ........................
               ti     &
             )
 !*====================================================================*
 !*                                                                    *
 !*  Compute the least magnitude solution of the quadratic equation    *
 !*  a * t**2 + b * t + c = 0. Here a, b, c and t are complex.         *
 !*                                         2                          *
 !*  Formeula used: t = 2c / (-b +/- sqrt (b  - 4ac)).                 *
 !*  This formula is valid for a=0 .                                   *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*  Input parameters:                                                 *
 !*  ================                                                  *
 !*      ar, ai   coefficient of t**2             REAL*8  ar, ai;      *
 !*      br, bi   coefficient of t                REAL*8  br, bi;      *
 !*      cr, ci   constant term                   REAL*8  cr, ci;      *
 !*                                                                    *
 !*  Output parameter:                                                 *
 !*  ================                                                  *
 !*      tr, ti   complex solution of minimal magnitude                *
 !*                                               REAL*8  tr, ti;      *
 !*                                                                    *
 !*====================================================================*
  IMPLICIT REAL*8 (A-H,O-Z) 

  pr = br * br - bi * bi
  pi = 2.d0 * br * bi                       !  p = b * b

  qr0 = ar * cr - ai * ci
  qi0 = ar * ci + ai * cr                   !  q = a * c

  pr = pr - 4.d0 * qr0       
  pi = pi - 4.d0 * qi0                      ! p = b * b - 4 * a * c

  h  = DSQRT (pr * pr + pi * pi)            ! q = sqrt (p)

  qr0 = h + pr
  if (qr0 > 0.d0) then
    qr0 = DSQRT (qr0 * 0.5d0)
  else
    qr0 = 0.d0
  end if

  qi0 = h - pr
  if (qi0 > 0.d0) then
    qi0 = DSQRT (qi0 * 0.5d0)
  else
    qi0 = 0.d0
  end if

  if (pi < 0.d0)  qi0 = -qi0

  h = qr0 * br + qi0 * bi      ! p = -b +/- q, choose sign for large}
                               ! magnitude  p
  if (h > 0.d0) then
    qr0 = -qr0
    qi0 = -qi0
  end if

  pr = qr0 - br
  pi = qi0 - bi
  h = pr * pr + pi * pi                       ! t = (2 * c) / p

  if (h == 0.d0) then
    tr = 0.d0
    ti = 0.d0
  else
    tr = 2.d0 * (cr * pr + ci * pi) / h
    ti = 2.d0 * (ci * pr - cr * pi) / h
  end if
  return
End Subroutine quadsolv


    Function ibauroot (n,      &    ! largest degree ...............
                       iu,     &    ! lowest degree ................
                       ar,     &    ! Real parts of the coefficients
                       ai,     &    ! Imaginary parts, coefficients 
                       x0r,    &    ! Real part of the root ........
                       x0i)         ! Imaginary part of the root ...
 !*====================================================================*
 !*                                                                    *
 !*  ibauroot computes one root of the polynomial P of degree n-iu:    *
 !*                                                 n-iu               *
 !*      P(x) = a(iu) + a(iu+1) * x + ... + a(n) * x                   *
 !*                                                                    *
 !*  with complex  coefficients a(i), i=iu, ..., n.                    *
 !*  This program uses Newton's method on the function P(x) / P'(x).   *
 !*  The iteration is stabilized using spiralization and extrapolation.*
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      n        int n;                                               *
 !*               Maximal degree of the polynomial ( >= 1 )            *
 !*      iu       int iu;                                              *
 !*               Index for the constant term of the polynomial,       *
 !*               n-iu is the degree of the polynomial with            *
 !*               coefficients a(iu), ..., a(n)                        *
 !*      ar, ai   REAL*8  ar(), ai();                                  *
 !*               Real and imaginary parts of the coefficients         *
 !*                                                                    *
 !*   Output parameter:                                                *
 !*   ================                                                 *
 !*      x0r,x0i  REAL*8  x0r, x0i;                                    *
 !*               Real and imaginary part of the computed root         *
 !*                                                                    *
 !*   Return value :                                                   *
 !*   =============                                                    *
 !*      = 0      all is ok                                            *
 !*      = 1      Division by zero                                     *
 !*      = 2      Iteration number ITERMAX exceeeded                   *
 !*      = 3      Improper input                                       *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Functions used:                                                  *
 !*   ==============                                                   *
 !*     chorner():   Complex Horner scheme                             *
 !*     comabs():    magnitude of a complex number                     *
 !*     icomdiv():    Complex division                                 *
 !*     quadsolv():  solves quadratic equations                        *
 !*                                                                    *
 !*   Constants used : ITERMAX,                                        *
 !*   ===============  QR, QI, MACH_EPS, EPS, EPSROOT, BETA            *
 !*                                                                    *
 !*====================================================================*
  IMPLICIT REAL*8 (A-H,O-Z)  
  Dimension ar(0:n), ai(0:n)

  Integer rc, result0, endit

  result0=2; endit=0; iter=0; i=0
  xoldr=0.d0; xoldi=0.d0; dxr=0.d0; dxi=0.d0; bdze=0.d0

  if (n < 1) then
    ibauroot=3
    return
  end if

  if (n - iu == 1) then                     ! Polynomial of degree 1
    call quadsolv (0.d0, 0.d0, ar(n), ai(n), ar(n-1), ai(n-1), x0r, x0i)
    ibauroot=0
    return
  end if

  if (n - iu == 2) then                     ! Polynomial of degree 2
    call quadsolv (ar(n),ai(n), ar(n-1),ai(n-1), ar(n-2),ai(n-2), x0r,x0i)
    ibauroot=0
    return
  end if

  xnewr = x0r;  xnewi = x0i
  endit = 0

  call chorner (n, iu, ar, ai, xnewr, xnewi,  &   ! Evaluate polynomial
                pr, pi, p1r, p1i, p2r, p2i, ss)

  iter=iter+1

  abs_pnew = comabs (pr, pi)
  if (abs_pnew < EPS) then
    ibauroot=0                              ! Starting value is a
    return                                  ! good approximation
  end if

  abs_pold = abs_pnew
  dzmin = BETA * (EPSROOT + comabs (xnewr, xnewi))

  do while (iter < ITERMAX)                 ! Bauhuber-Iteration

    abs_p1new = comabs (p1r, p1i)

    if (abs_pnew > abs_pold) then           ! Spiralization step

      i = 0                                 ! dx = dx * q
      iter=iter+1
      temp = dxr

      dxr = QR * dxr - QI * dxi
      dxi = QR * dxi + QI * temp
    
    else
      
      dzmax = 1 + comabs(xnewr, xnewi)
      h1 = p1r * p1r - p1i * p1i - pr * p2r + pi * p2i
      h2 = 2.d0 * p1r * p1i - pr * p2i - pi * p2r
      if (abs_p1new > 10.d0 * ss.and.comabs(h1, h2) > 100.d0 * ss * ss) then
        i=i+1
        if (i > 2) i = 2
        tempr = pr * p1r - pi * p1i
        tempi = pr * p1i + pi * p1r

        rc = icomdiv(-tempr, -tempi, h1, h2, dxr, dxi)
        
	if (rc.ne.0) then
          ibauroot = 1
          return
        end if

        if (comabs(dxr, dxi) > dzmax) then
          temp = dzmax / comabs(dxr, dxi)          ! Newton step
          dxr = dxr * temp
          dxi = dxi * temp
          i = 0
        end if

        if (i == 2.and.comabs(dxr, dxi) < dzmin / EPSROOT.and.comabs(dxr, dxi) > 0.d0) then
          i = 0                                    ! Extrapolation step
          rc = icomdiv(xnewr - xoldr, xnewi - xoldi, dxr, dxi, h3, h4)
          if (rc.ne.0) then
            ibauroot = 1
            return
          end if

          h3 = h3 + 1.d0
          h1 = h3 * h3 - h4 * h4
          h2 = 2.d0 * h3 * h4
          rc = icomdiv (dxr, dxi, h1, h2, h3, h4)
          if (rc.ne.0) then
            ibauroot = 1
            return
          end if

          if (comabs(h3, h4) < 50.d0 * dzmin) then
            dxr = dxr + h3
            dxi = dxi + h4
          end if
        end if

        xoldr = xnewr
        xoldi = xnewi
        abs_pold = abs_pnew
      else
        i = 0                               ! Close to a saddle point
        h = dzmax / abs_pnew
        dxr = h * pr
        dxi = h * pi

        xoldr = xnewr
        xoldi = xnewi
        abs_pold = abs_pnew

        do while(DABS(comabs(u,v)/abs_pnew-1.d0) < EPSROOT)

          call chorner (n, iu, ar, ai, xnewr+dxr, xnewi+dxi, u, v, h, h1, h2, h3, h4)
          iter=iter+1

          dxr = dxr + dxr
          dxi = dxi + dxi                   ! dx = dx * 2

        end do
      end if
    end if

    if (endit.ne.0) then
      if (comabs(dxr, dxi) < 0.1d0 * bdze) then
        xnewr = xnewr + dxr
        xnewi = xnewi + dxi
      end if

      result0 = 0
      goto 10                                 ! stop iteration
    else
      xnewr = xoldr + dxr
      xnewi = xoldi + dxi
      dzmin = BETA * (EPSROOT + comabs (xnewr, xnewi))
      call chorner (n, iu, ar, ai, xnewr, xnewi, pr, pi, p1r, p1i, p2r, p2i, ss)
      iter=iter+1
      abs_pnew = comabs(pr,pi)

      if (abs_pnew == 0.d0) then
        result0 = 0
        goto 10
      end if

      if (comabs(dxr, dxi) < dzmin.or.abs_pnew < EPS) then
        endit = 1
        bdze = comabs(dxr, dxi)
      end if
    end if
  end do  !End Bauhuber iteration

10 x0r = xnewr
   x0i = xnewi
  ibauroot=result0
  
  return
End Function ibauroot


  Subroutine polydiv (n,      &      ! maximal degree ...............
                     iu,      &      ! minimal degree ...............
                     ar,      &      ! Real parts of coefficients ...
                     ai,      &      ! Imaginary parts, coefficients 
                     x0r,     &      ! Real part of x ...............
                     x0i)            ! Imaginary part of x ..........
 !*====================================================================*
 !*                                                                    *
 !*  polydiv computes the coefficients of the polynomial Q, with       *
 !*  P(x) = Q(x) * (x - x0), where x0 is a computed root of P.         *
 !*  Both P and Q and x0 may be complex.                               *
 !*                                                                    *
 !*====================================================================*
 !*                                                                    *
 !*   Input parameters:                                                *
 !*   ================                                                 *
 !*      n        int n;                                               *
 !*               Highest degree of the polynomial ( >= 1 )            *
 !*      ar, ai   REAL*8  ar(), ai();                                  *
 !*               Real and imaginary parts of the coefficienten of P   *
 !*               of degree n-iu, with a(iu), ..., a(n)                *
 !*      x0r,x0i  REAL*8  x0r, x0i;                                    *
 !*               Real and imaginary parts of the root x0              *
 !*                                                                    *
 !*   Output parameter:                                                *
 !*   ================                                                 *
 !*      ar, ai   REAL*8  ar(), ai();                                  *
 !*               Real and imaginary parts of the coefficients         *
 !*               ar(iu+1),..,ar(n) of the remainder polynomial Q      *
 !*                                                                    *
 !*====================================================================*
   IMPLICIT REAL*8 (A-H,O-Z) 
   Dimension ar(0:n), ai(0:n)
   do i = n - 1, iu+1, -1
     temp = ar(i+1)
     ar(i) = ar(i) + temp * x0r - ai(i+1) * x0i
     ai(i) = ai(i) + ai(i+1) * x0r + temp * x0i
   end do
 End Subroutine polydiv 

END MODULE

! ------------------------- END fbauhube.f90  -------------------------
