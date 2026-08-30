!*****************************************************************************80
!*****************************************************************************80
!*****************************************************************************80
function algdiv ( a, b )

!*****************************************************************************80
!
!! ALGDIV computes ln ( Gamma ( B ) / Gamma ( A + B ) ) when 8 <= B.
!
!  Discussion:
!
!    In this algorithm, DEL(X) is the function defined by
!
!      ln ( Gamma(X) ) = ( X - 0.5 ) * ln ( X ) - X + 0.5 * ln ( 2 * PI ) 
!                      + DEL(X).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, define the arguments.
!
!    Output, real (real64) ALGDIV, the value of ln(Gamma(B)/Gamma(A+B)).
!
  use iso_fortran_env, only: real64
  implicit none
  

  real (real64) a
  real (real64) algdiv
  real (real64) alnrel
  real (real64) b
  real (real64) c
  real (real64), parameter :: c0 =  0.833333333333333D-01
  real (real64), parameter :: c1 = -0.277777777760991D-02
  real (real64), parameter :: c2 =  0.793650666825390D-03
  real (real64), parameter :: c3 = -0.595202931351870D-03
  real (real64), parameter :: c4 =  0.837308034031215D-03
  real (real64), parameter :: c5 = -0.165322962780713D-02
  real (real64) d
  real (real64) h
  real (real64) s11
  real (real64) s3
  real (real64) s5
  real (real64) s7
  real (real64) s9
  real (real64) t
  real (real64) u
  real (real64) v
  real (real64) w
  real (real64) x
  real (real64) x2

  if ( b < a ) then
    h = b / a
    c = 1.0D+00 / ( 1.0D+00 + h )
    x = h / ( 1.0D+00 + h )
    d = a + ( b - 0.5D+00 )
  else
    h = a / b
    c = h / ( 1.0D+00 + h )
    x = 1.0D+00 / ( 1.0D+00 + h )
    d = b + ( a - 0.5D+00 )
  end if
!
!  Set SN = (1 - X**N)/(1 - X).
!
  x2 = x * x
  s3 = 1.0D+00 + ( x + x2 )
  s5 = 1.0D+00 + ( x + x2 * s3 )
  s7 = 1.0D+00 + ( x + x2 * s5 )
  s9 = 1.0D+00 + ( x + x2 * s7 )
  s11 = 1.0D+00 + ( x + x2 * s9 )
!
!  Set W = DEL(B) - DEL(A + B).
!
  t = ( 1.0D+00 / b )**2
  w = (((( &
          c5 * s11  * t &
        + c4 * s9 ) * t &
        + c3 * s7 ) * t &
        + c2 * s5 ) * t &
        + c1 * s3 ) * t &
        + c0

  w = w * ( c / b )
!
!  Combine the results.
!
  u = d * alnrel ( a / b )
  v = a * ( log ( b ) - 1.0D+00 )

  if ( v < u ) then
    algdiv = ( w - v ) - u
  else
    algdiv = ( w - u ) - v
  end if

  return
end
function alnrel ( a )

!*****************************************************************************80
!
!! ALNREL evaluates the function ln ( 1 + A ).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, the argument.
!
!    Output, real (real64) ALNREL, the value of ln ( 1 + A ).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) alnrel
  real (real64), parameter :: p1 = -0.129418923021993D+01
  real (real64), parameter :: p2 =  0.405303492862024D+00
  real (real64), parameter :: p3 = -0.178874546012214D-01
  real (real64), parameter :: q1 = -0.162752256355323D+01
  real (real64), parameter :: q2 =  0.747811014037616D+00
  real (real64), parameter :: q3 = -0.845104217945565D-01
  real (real64) t
  real (real64) t2
  real (real64) w
  real (real64) x

  if ( abs ( a ) <= 0.375D+00 ) then

    t = a / ( a +  2.0D+00  )
    t2 = t * t

    w = ((( p3 * t2 + p2 ) * t2 + p1 ) * t2 + 1.0D+00 ) &
      / ((( q3 * t2 + q2 ) * t2 + q1 ) * t2 + 1.0D+00 )

    alnrel =  2.0D+00  * t * w

  else

    x = 1.0D+00 + real ( a, kind = real64 )
    alnrel = log ( x )

  end if

  return
end
function apser ( a, b, x, eps )

!*****************************************************************************80
!
!! APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
!
!  Discussion:
!
!    APSER is used only for cases where
!
!      A <= min ( EPS, EPS * B ), 
!      B * X <= 1, and 
!      X <= 0.5.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, X, the parameters of the
!    incomplete beta ratio.
!
!    Input, real (real64) EPS, a tolerance.
!
!    Output, real (real64) APSER, the computed value of the
!    incomplete beta ratio.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) aj
  real (real64) apser
  real (real64) b
  real (real64) bx
  real (real64) c
  real (real64) eps
  real (real64), parameter :: g = 0.577215664901533D+00
  real (real64) j
  real (real64) psi
  real (real64) s
  real (real64) t
  real (real64) tol
  real (real64) x

  bx = b * x
  t = x - bx

  if ( b * eps <= 0.02D+00 ) then
    c = log ( x ) + psi ( b ) + g + t
  else
    c = log ( bx ) + g + t
  end if

  tol = 5.0D+00 * eps * abs ( c )
  j = 1.0D+00
  s = 0.0D+00

  do

    j = j + 1.0D+00
    t = t * ( x - bx / j )
    aj = t / j
    s = s + aj

    if ( abs ( aj ) <= tol ) then
      exit
    end if

  end do

  apser = -a * ( c + s )

  return
end
function bcorr ( a0, b0 )

!*****************************************************************************80
!
!! BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0).
!
!  Discussion:
!
!    The function DEL(A) is a remainder term that is used in the expression:
!
!      ln ( Gamma ( A ) ) = ( A - 0.5 ) * ln ( A ) 
!        - A + 0.5 * ln ( 2 * PI ) + DEL ( A ),
!
!    or, in other words, DEL ( A ) is defined as:
!
!      DEL ( A ) = ln ( Gamma ( A ) ) - ( A - 0.5 ) * ln ( A ) 
!        + A + 0.5 * ln ( 2 * PI ).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A0, B0, the arguments.
!    It is assumed that 8 <= A0 and 8 <= B0.
!
!    Output, real (real64) BCORR, the value of the function.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a0
  real (real64) b
  real (real64) b0
  real (real64) bcorr
  real (real64) c
  real (real64), parameter :: c0 =  0.833333333333333D-01
  real (real64), parameter :: c1 = -0.277777777760991D-02
  real (real64), parameter :: c2 =  0.793650666825390D-03
  real (real64), parameter :: c3 = -0.595202931351870D-03
  real (real64), parameter :: c4 =  0.837308034031215D-03
  real (real64), parameter :: c5 = -0.165322962780713D-02
  real (real64) h
  real (real64) s11
  real (real64) s3
  real (real64) s5
  real (real64) s7
  real (real64) s9
  real (real64) t
  real (real64) w
  real (real64) x
  real (real64) x2

  a = min ( a0, b0 )
  b = max ( a0, b0 )

  h = a / b
  c = h / ( 1.0D+00 + h )
  x = 1.0D+00 / ( 1.0D+00 + h )
  x2 = x * x
!
!  Set SN = (1 - X**N)/(1 - X)
!
  s3 = 1.0D+00 + ( x + x2 )
  s5 = 1.0D+00 + ( x + x2 * s3 )
  s7 = 1.0D+00 + ( x + x2 * s5 )
  s9 = 1.0D+00 + ( x + x2 * s7 )
  s11 = 1.0D+00 + ( x + x2 * s9 )
!
!  Set W = DEL(B) - DEL(A + B)
!
  t = ( 1.0D+00 / b )**2

  w = (((( &
       c5 * s11  * t &
     + c4 * s9 ) * t &
     + c3 * s7 ) * t &
     + c2 * s5 ) * t &
     + c1 * s3 ) * t &
     + c0

  w = w * ( c / b )
!
!  Compute  DEL(A) + W.
!
  t = ( 1.0D+00 / a )**2

  bcorr = ((((( &
         c5   * t &
       + c4 ) * t &
       + c3 ) * t &
       + c2 ) * t &
       + c1 ) * t &
       + c0 ) / a + w

  return
end
function beta ( a, b )

!*****************************************************************************80
!
!! BETA evaluates the beta function.
!
!  Modified:
!
!    03 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real (real64) A, B, the arguments of the beta function.
!
!    Output, real (real64) BETA, the value of the beta function.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) b
  real (real64) beta
  real (real64) beta_log

  beta = exp ( beta_log ( a, b ) )

  return
end
function beta_asym ( a, b, lambda, eps )

!*****************************************************************************80
!
!! BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters of the function.
!    A and B should be nonnegative.  It is assumed that both A and B
!    are greater than or equal to 15.
!
!    Input, real (real64) LAMBDA, the value of ( A + B ) * Y - B.
!    It is assumed that 0 <= LAMBDA.
!
!    Input, real (real64) EPS, the tolerance.
!
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: num = 20

  real (real64) a
  real (real64) a0(num+1)
  real (real64) b
  real (real64) b0(num+1)
  real (real64) bcorr
  real (real64) beta_asym
  real (real64) bsum
  real (real64) c(num+1)
  real (real64) d(num+1)
  real (real64) dsum
  real (real64), parameter :: e0 = 1.12837916709551D+00
  real (real64), parameter :: e1 = 0.353553390593274D+00
  real (real64) eps
  real (real64) error_fc
  real (real64) f
  real (real64) h
  real (real64) h2
  real (real64) hn
  integer i
  integer j
  real (real64) j0
  real (real64) j1
  real (real64) lambda
  integer m
  integer mm1
  integer mmj
  integer n
  integer np1
  real (real64) r
  real (real64) r0
  real (real64) r1
  real (real64) rlog1
  real (real64) s
  real (real64) sum1
  real (real64) t
  real (real64) t0
  real (real64) t1
  real (real64) u
  real (real64) w
  real (real64) w0
  real (real64) z
  real (real64) z0
  real (real64) z2
  real (real64) zn
  real (real64) znm1

  beta_asym = 0.0D+00

  if ( a < b ) then
    h = a / b
    r0 = 1.0D+00 / ( 1.0D+00 + h )
    r1 = ( b - a ) / b
    w0 = 1.0D+00 / sqrt ( a * ( 1.0D+00 + h ))
  else
    h = b / a
    r0 = 1.0D+00 / ( 1.0D+00 + h )
    r1 = ( b - a ) / a
    w0 = 1.0D+00 / sqrt ( b * ( 1.0D+00 + h ))
  end if

  f = a * rlog1 ( - lambda / a ) + b * rlog1 ( lambda / b )
  t = exp ( - f )
  if ( t == 0.0D+00 ) then
    return
  end if

  z0 = sqrt ( f )
  z = 0.5D+00 * ( z0 / e1 )
  z2 = f + f

  a0(1) = ( 2.0D+00 / 3.0D+00 ) * r1
  c(1) = -0.5D+00 * a0(1)
  d(1) = -c(1)
  j0 = ( 0.5D+00 / e0 ) * error_fc ( 1, z0 )
  j1 = e1
  sum1 = j0 + d(1) * w0 * j1

  s = 1.0D+00
  h2 = h * h
  hn = 1.0D+00
  w = w0
  znm1 = z
  zn = z2

  do n = 2, num, 2

    hn = h2 * hn
    a0(n) = 2.0D+00 * r0 * ( 1.0D+00 + h * hn ) &
      / ( n + 2.0D+00 )
    np1 = n + 1
    s = s + hn
    a0(np1) = 2.0D+00 * r1 * s / ( n + 3.0D+00 )

    do i = n, np1

      r = -0.5D+00 * ( i + 1.0D+00 )
      b0(1) = r * a0(1)
      do m = 2, i
        bsum = 0.0D+00
        mm1 = m - 1
        do j = 1, mm1
          mmj = m - j
          bsum = bsum + ( j * r - mmj ) * a0(j) * b0(mmj)
        end do
        b0(m) = r * a0(m) + bsum / m
      end do

      c(i) = b0(i) / ( i + 1.0D+00 )

      dsum = 0.0
      do j = 1, i-1
        dsum = dsum + d(i-j) * c(j)
      end do
      d(i) = - ( dsum + c(i) )

    end do

    j0 = e1 * znm1 + ( n - 1.0D+00 ) * j0
    j1 = e1 * zn + n * j1
    znm1 = z2 * znm1
    zn = z2 * zn
    w = w0 * w
    t0 = d(n) * w * j0
    w = w0 * w
    t1 = d(np1) * w * j1
    sum1 = sum1 + ( t0 + t1 )

    if ( ( abs ( t0 ) + abs ( t1 )) <= eps * sum1 ) then
      u = exp ( - bcorr ( a, b ) )
      beta_asym = e0 * t * u * sum1
      return
    end if

  end do

  u = exp ( - bcorr ( a, b ) )
  beta_asym = e0 * t * u * sum1

  return
end
function beta_frac ( a, b, x, y, lambda, eps )

!*****************************************************************************80
!
!! BETA_FRAC evaluates a continued fraction expansion for IX(A,B).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters of the function.
!    A and B should be nonnegative.  It is assumed that both A and
!    B are greater than 1.
!
!    Input, real (real64) X, Y.  X is the argument of the
!    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
!
!    Input, real (real64) LAMBDA, the value of ( A + B ) * Y - B.
!
!    Input, real (real64) EPS, a tolerance.
!
!    Output, real (real64) BETA_FRAC, the value of the continued
!    fraction approximation for IX(A,B).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) alpha
  real (real64) an
  real (real64) anp1
  real (real64) b
  real (real64) beta
  real (real64) beta_frac
  real (real64) beta_rcomp
  real (real64) bn
  real (real64) bnp1
  real (real64) c
  real (real64) c0
  real (real64) c1
  real (real64) e
  real (real64) eps
  real (real64) lambda
  real (real64) n
  real (real64) p
  real (real64) r
  real (real64) r0
  real (real64) s
  real (real64) t
  real (real64) w
  real (real64) x
  real (real64) y
  real (real64) yp1

  beta_frac = beta_rcomp ( a, b, x, y )

  if ( beta_frac == 0.0D+00 ) then
    return
  end if

  c = 1.0D+00 + lambda
  c0 = b / a
  c1 = 1.0D+00 + 1.0D+00 / a
  yp1 = y + 1.0D+00

  n = 0.0D+00
  p = 1.0D+00
  s = a + 1.0D+00
  an = 0.0D+00
  bn = 1.0D+00
  anp1 = 1.0D+00
  bnp1 = c / c1
  r = c1 / c
!
!  Continued fraction calculation.
!
  do

    n = n + 1.0D+00
    t = n / a
    w = n * ( b - n ) * x
    e = a / s
    alpha = ( p * ( p + c0 ) * e * e ) * ( w * x )
    e = ( 1.0D+00 + t ) / ( c1 + t + t )
    beta = n + w / s + e * ( c + n * yp1 )
    p = 1.0D+00 + t
    s = s +  2.0D+00 
!
!  Update AN, BN, ANP1, and BNP1.
!
    t = alpha * an + beta * anp1
    an = anp1
    anp1 = t
    t = alpha * bn + beta * bnp1
    bn = bnp1
    bnp1 = t

    r0 = r
    r = anp1 / bnp1

    if ( abs ( r - r0 ) <= eps * r ) then
      beta_frac = beta_frac * r
      exit
    end if
!
!  Rescale AN, BN, ANP1, and BNP1.
!
    an = an / bnp1
    bn = bn / bnp1
    anp1 = r
    bnp1 = 1.0D+00

  end do

  return
end
subroutine beta_grat ( a, b, x, y, w, eps, ierr )

!*****************************************************************************80
!
!! BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters of the function.
!    A and B should be nonnegative.  It is assumed that 15 <= A 
!    and B <= 1, and that B is less than A.
!
!    Input, real (real64) X, Y.  X is the argument of the
!    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
!
!    Input/output, real (real64) W, a quantity to which the
!    result of the computation is to be added on output.
!
!    Input, real (real64) EPS, a tolerance.
!
!    Output, integer IERR, an error flag, which is 0 if no error
!    was detected.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) algdiv
  real (real64) alnrel
  real (real64) b
  real (real64) bm1
  real (real64) bp2n
  real (real64) c(30)
  real (real64) cn
  real (real64) coef
  real (real64) d(30)
  real (real64) dj
  real (real64) eps
  real (real64) gam1
  integer i
  integer ierr
  real (real64) j
  real (real64) l
  real (real64) lnx
  integer n
  real (real64) n2
  real (real64) nu
  real (real64) p
  real (real64) q
  real (real64) r
  real (real64) s
  real (real64) sum1
  real (real64) t
  real (real64) t2
  real (real64) u
  real (real64) v
  real (real64) w
  real (real64) x
  real (real64) y
  real (real64) z

  bm1 = ( b - 0.5D+00 ) - 0.5D+00
  nu = a + 0.5D+00 * bm1

  if ( y <= 0.375D+00 ) then
    lnx = alnrel ( - y )
  else
    lnx = log ( x )
  end if

  z = -nu * lnx

  if ( b * z == 0.0D+00 ) then
    ierr = 1
    return
  end if
!
!  Computation of the expansion.
!
!  Set R = EXP(-Z)*Z**B/GAMMA(B)
!
  r = b * ( 1.0D+00 + gam1 ( b ) ) * exp ( b * log ( z ))
  r = r * exp ( a * lnx ) * exp ( 0.5D+00 * bm1 * lnx )
  u = algdiv ( b, a ) + b * log ( nu )
  u = r * exp ( - u )

  if ( u == 0.0D+00 ) then
    ierr = 1
    return
  end if

  call gamma_rat1 ( b, z, r, p, q, eps )

  v = 0.25D+00 * ( 1.0D+00 / nu )**2
  t2 = 0.25D+00 * lnx * lnx
  l = w / u
  j = q / r
  sum1 = j
  t = 1.0D+00
  cn = 1.0D+00
  n2 = 0.0D+00

  do n = 1, 30

    bp2n = b + n2
    j = ( bp2n * ( bp2n + 1.0D+00 ) * j &
      + ( z + bp2n + 1.0D+00 ) * t ) * v
    n2 = n2 +  2.0D+00 
    t = t * t2
    cn = cn / ( n2 * ( n2 + 1.0D+00 ))
    c(n) = cn
    s = 0.0D+00

    coef = b - n
    do i = 1, n-1
      s = s + coef * c(i) * d(n-i)
      coef = coef + b
    end do

    d(n) = bm1 * cn + s / n
    dj = d(n) * j
    sum1 = sum1 + dj

    if ( sum1 <= 0.0D+00 ) then
      ierr = 1
      return
    end if

    if ( abs ( dj ) <= eps * ( sum1 + l ) ) then
      ierr = 0
      w = w + u * sum1
      return
    end if

  end do

  ierr = 0
  w = w + u * sum1

  return
end
subroutine beta_inc ( a, b, x, y, w, w1, ierr )

!*****************************************************************************80
!
!! BETA_INC evaluates the incomplete beta function IX(A,B).
!
!  Author:
!
!    Alfred Morris,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!  Modified February 2024 (error handling):
!    Seth Siefken
!    U.S. Geological Survey
!    Helena, Montana
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters of the function.
!    A and B should be nonnegative.
!
!    Input, real (real64) X, Y.  X is the argument of the
!    function, and should satisy 0 <= X <= 1.  Y should equal 1 - X.
!
!    Output, real (real64) W, W1, the values of IX(A,B) and
!    1-IX(A,B).
!
!    Output, integer IERR, the error flag.
!    0, no error was detected.
!    1, A or B is negative;
!    2, A = B = 0;
!    3, X < 0 or 1 < X;
!    4, Y < 0 or 1 < Y;
!    5, X + Y /= 1;
!    6, X = A = 0;
!    7, Y = B = 0.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a0
  real (real64) apser
  real (real64) b
  real (real64) b0
  real (real64) beta_asym
  real (real64) beta_frac
  real (real64) beta_pser
  real (real64) beta_up
  real (real64) eps
  real (real64) fpser
  integer ierr
  integer ierr1
  integer ind
  real (real64) lambda
  integer n
  real (real64) t
  real (real64) w
  real (real64) w1
  real (real64) x
  real (real64) x0
  real (real64) y
  real (real64) y0
  real (real64) z

  eps = epsilon ( eps )
  w = 0.0D+00
  w1 = 0.0D+00

  if ( a < 0.0D+00 .or. b < 0.0D+00 ) then
    error stop 'BETA_INC - Fatal error: A or B is negative'
  end if

  if ( a == 0.0D+00 .and. b == 0.0D+00 ) then
    error stop 'BETA_INC - Fatal error: A = B = 0'
  end if

  if ( x < 0.0D+00 .or. 1.0D+00 < x ) then
    error stop  'BETA_INC - Fatal error: X < 0 or 1 < X'
  end if

  if ( y < 0.0D+00 .or. 1.0D+00 < y ) then
    error stop  'BETA_INC - Fatal error: Y < 0 or 1 < Y'
  end if

  z = ( ( x + y ) - 0.5D+00 ) - 0.5D+00

  if ( 3.0D+00 * eps < abs ( z ) ) then
    error stop  'BETA_INC - Fatal error: X + Y /= 1'
  end if

  ierr = 0

  if ( x == 0.0D+00 ) then
    w = 0.0D+00
    w1 = 1.0D+00
    if ( a == 0.0D+00 ) then
    error stop  'BETA_INC - Fatal error: X = A = 0'
    end if
    return
  end if

  if ( y == 0.0D+00 ) then
    if ( b == 0.0D+00 ) then
    error stop  'BETA_INC - Fatal error: Y = B = 0'
    end if
    w = 1.0D+00
    w1 = 0.0D+00
    return
  end if

  if ( a == 0.0D+00 ) then
    w = 1.0D+00
    w1 = 0.0D+00
    return
  end if

  if ( b == 0.0D+00 ) then
    w = 0.0D+00
    w1 = 1.0D+00
    return
  end if

  eps = max ( eps, 1.0D-15 )

  if ( max ( a, b ) < 0.001D+00 * eps ) then
    go to 260
  end if

  ind = 0
  a0 = a
  b0 = b
  x0 = x
  y0 = y

  if ( 1.0D+00 < min ( a0, b0 ) ) then
    go to 40
  end if
!
!  Procedure for A0 <= 1 or B0 <= 1
!
  if ( 0.5D+00 < x ) then
    ind = 1
    a0 = b
    b0 = a
    x0 = y
    y0 = x
  end if

  if ( b0 < min ( eps, eps * a0 ) ) then
    go to 90
  end if

  if ( a0 < min ( eps, eps * b0 ) .and. b0 * x0 <= 1.0D+00 ) then
    go to 100
  end if

  if ( 1.0D+00 < max ( a0, b0 ) ) then
    go to 20
  end if

  if ( min ( 0.2D+00, b0 ) <= a0 ) then
    go to 110
  end if

  if ( x0**a0 <= 0.9D+00 ) then
    go to 110
  end if

  if ( 0.3D+00 <= x0 ) then
    go to 120
  end if

  n = 20
  go to 140

20 continue

  if ( b0 <= 1.0D+00 ) then
    go to 110
  end if

  if ( 0.3D+00 <= x0 ) then
    go to 120
  end if

  if ( 0.1D+00 <= x0 ) then
    go to 30
  end if

  if ( ( x0 * b0 )**a0 <= 0.7D+00 ) then
    go to 110
  end if

30 continue

  if ( 15.0D+00 < b0 ) then
    go to 150
  end if

  n = 20
  go to 140
!
!  PROCEDURE for 1 < A0 and 1 < B0.
!
40 continue

  if ( a <= b ) then
    lambda = a - ( a + b ) * x
  else
    lambda = ( a + b ) * y - b
  end if

  if ( lambda < 0.0D+00 ) then
    ind = 1
    a0 = b
    b0 = a
    x0 = y
    y0 = x
    lambda = abs ( lambda )
  end if

70 continue

  if ( b0 < 40.0D+00 .and. b0 * x0 <= 0.7D+00 ) then
    go to 110
  end if

  if ( b0 < 40.0D+00 ) then
    go to 160
  end if

  if ( b0 < a0 ) then
    go to 80
  end if

  if ( a0 <= 100.0D+00 ) then
    go to 130
  end if

  if ( 0.03D+00 * a0 < lambda ) then
    go to 130
  end if

  go to 200

80 continue

  if ( b0 <= 100.0D+00 ) then
    go to 130
  end if

  if ( 0.03D+00 * b0 < lambda ) then
    go to 130
  end if

  go to 200
!
!  Evaluation of the appropriate algorithm.
!
90 continue

  w = fpser ( a0, b0, x0, eps )
  w1 = 0.5D+00 + ( 0.5D+00 - w )
  go to 250

100 continue

  w1 = apser ( a0, b0, x0, eps )
  w = 0.5D+00 + ( 0.5D+00 - w1 )
  go to 250

110 continue

  w = beta_pser ( a0, b0, x0, eps )
  w1 = 0.5D+00 + ( 0.5D+00 - w )
  go to 250

120 continue

  w1 = beta_pser ( b0, a0, y0, eps )
  w = 0.5D+00 + ( 0.5D+00 - w1 )
  go to 250

130 continue

  w = beta_frac ( a0, b0, x0, y0, lambda, 15.0D+00 * eps )
  w1 = 0.5D+00 + ( 0.5D+00 - w )
  go to 250

140 continue

  w1 = beta_up ( b0, a0, y0, x0, n, eps )
  b0 = b0 + n

150 continue

  call beta_grat ( b0, a0, y0, x0, w1, 15.0D+00 * eps, ierr1 )
  w = 0.5D+00 + ( 0.5D+00 - w1 )
  go to 250

160 continue

  n = b0
  b0 = b0 - n

  if ( b0 == 0.0D+00 ) then
    n = n - 1
    b0 = 1.0D+00
  end if

170 continue

  w = beta_up ( b0, a0, y0, x0, n, eps )

  if ( x0 <= 0.7D+00 ) then
    w = w + beta_pser ( a0, b0, x0, eps )
    w1 = 0.5D+00 + ( 0.5D+00 - w )
    go to 250
  end if

  if ( a0 <= 15.0D+00 ) then
    n = 20
    w = w + beta_up ( a0, b0, x0, y0, n, eps )
    a0 = a0 + n
  end if

190 continue

  call beta_grat ( a0, b0, x0, y0, w, 15.0D+00 * eps, ierr1 )
  w1 = 0.5D+00 + ( 0.5D+00 - w )
  go to 250

200 continue

  w = beta_asym ( a0, b0, lambda, 100.0D+00 * eps )
  w1 = 0.5D+00 + ( 0.5D+00 - w )
  go to 250
!
!  Termination of the procedure.
!
250 continue

  if ( ind /= 0 ) then
    t = w
    w = w1
    w1 = t
  end if

  return
!
!  Procedure for A and B < 0.001 * EPS
!
260 continue

  w = b / ( a + b )
  w1 = a / ( a + b )

  return
end
subroutine beta_inc_values ( n_data, a, b, x, fx )

!*****************************************************************************80
!
!! BETA_INC_VALUES returns some values of the incomplete Beta function.
!
!  Discussion:
!
!    The incomplete Beta function may be written
!
!      BETA_INC(A,B,X) = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
!                      / Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
!
!    Thus,
!
!      BETA_INC(A,B,0.0) = 0.0
!      BETA_INC(A,B,1.0) = 1.0
!
!    Note that in Mathematica, the expressions:
!
!      BETA[A,B]   = Integral (0 to 1) T**(A-1) * (1-T)**(B-1) dT
!      BETA[X,A,B] = Integral (0 to X) T**(A-1) * (1-T)**(B-1) dT
!
!    and thus, to evaluate the incomplete Beta function requires:
!
!      BETA_INC(A,B,X) = BETA[X,A,B] / BETA[A,B]
!
!  Modified:
!
!    17 February 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    US Department of Commerce, 1964.
!
!    Karl Pearson,
!    Tables of the Incomplete Beta Function,
!    Cambridge University Press, 1968.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real (real64) A, B, X, the arguments of the function.
!
!    Output, real (real64) FX, the value of the function.
!
  use iso_fortran_env, only: real64
  implicit none

  integer, parameter :: n_max = 30

  real (real64) a
  real (real64), save, dimension ( n_max ) :: a_vec = (/ &
     0.5D+00,  0.5D+00,  0.5D+00,  1.0D+00, &
     1.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, &
     2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
     2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
     2.0D+00,  5.5D+00, 10.0D+00, 10.0D+00, &
    10.0D+00, 10.0D+00, 20.0D+00, 20.0D+00, &
    20.0D+00, 20.0D+00, 20.0D+00, 30.0D+00, &
    30.0D+00, 40.0D+00 /)
  real (real64) b
  real (real64), save, dimension ( n_max ) :: b_vec = (/ &
     0.5D+00,  0.5D+00,  0.5D+00,  0.5D+00, &
     0.5D+00,  0.5D+00,  0.5D+00,  1.0D+00, &
     2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
     2.0D+00,  2.0D+00,  2.0D+00,  2.0D+00, &
     2.0D+00,  5.0D+00,  0.5D+00,  5.0D+00, &
     5.0D+00, 10.0D+00,  5.0D+00, 10.0D+00, &
    10.0D+00, 20.0D+00, 20.0D+00, 10.0D+00, &
    10.0D+00, 20.0D+00 /)
  real (real64) fx
  real (real64), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0637686D+00, 0.2048328D+00, 1.0000000D+00, 0.0D+00,       &
    0.0050126D+00, 0.0513167D+00, 0.2928932D+00, 0.5000000D+00, &
    0.028D+00,     0.104D+00,     0.216D+00,     0.352D+00,     &
    0.500D+00,     0.648D+00,     0.784D+00,     0.896D+00,     &
    0.972D+00,     0.4361909D+00, 0.1516409D+00, 0.0897827D+00, &
    1.0000000D+00, 0.5000000D+00, 0.4598773D+00, 0.2146816D+00, &
    0.9507365D+00, 0.5000000D+00, 0.8979414D+00, 0.2241297D+00, &
    0.7586405D+00, 0.7001783D+00 /)
  integer n_data
  real (real64) x
  real (real64), save, dimension ( n_max ) :: x_vec = (/ &
    0.01D+00, 0.10D+00, 1.00D+00, 0.0D+00,  &
    0.01D+00, 0.10D+00, 0.50D+00, 0.50D+00, &
    0.1D+00,  0.2D+00,  0.3D+00,  0.4D+00,  &
    0.5D+00,  0.6D+00,  0.7D+00,  0.8D+00,  &
    0.9D+00,  0.50D+00, 0.90D+00, 0.50D+00, &
    1.00D+00, 0.50D+00, 0.80D+00, 0.60D+00, &
    0.80D+00, 0.50D+00, 0.60D+00, 0.70D+00, &
    0.80D+00, 0.70D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    a = 0.0D+00
    b = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    a = a_vec(n_data)
    b = b_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
function beta_log ( a0, b0 )

!*****************************************************************************80
!
!! BETA_LOG evaluates the logarithm of the beta function.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software, 
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A0, B0, the parameters of the function.
!    A0 and B0 should be nonnegative.
!
!    Output, real (real64) BETA_LOG, the value of the logarithm
!    of the Beta function.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a0
  real (real64) algdiv
  real (real64) alnrel
  real (real64) b
  real (real64) b0
  real (real64) bcorr
  real (real64) beta_log
  real (real64) c
  real (real64), parameter :: e = 0.918938533204673D+00
  real (real64) gamma_log
  real (real64) gsumln
  real (real64) h
  integer i
  integer n
  real (real64) u
  real (real64) v
  real (real64) w
  real (real64) z

  a = min ( a0, b0 )
  b = max ( a0, b0 )
!
!  8 < A.
!
  if ( 8.0D+00 <= a ) then

    w = bcorr ( a, b )
    h = a / b
    c = h / ( 1.0D+00 + h )
    u = - ( a - 0.5D+00 ) * log ( c )
    v = b * alnrel ( h )

    if ( v < u ) then
      beta_log = ((( -0.5D+00 * log ( b ) + e ) + w ) - v ) - u
    else
      beta_log = ((( -0.5D+00 * log ( b ) + e ) + w ) - u ) - v
    end if

    return
  end if
!
!  Procedure when A < 1
!
  if ( a < 1.0D+00 ) then

    if ( b < 8.0D+00 ) then
      beta_log = gamma_log ( a ) + ( gamma_log ( b ) - gamma_log ( a + b ) )
    else
      beta_log = gamma_log ( a ) + algdiv ( a, b )
    end if

    return

  end if
!
!  Procedure when 1 <= A < 8
!
  if ( 2.0D+00 < a ) then
    go to 40
  end if

  if ( b <= 2.0D+00 ) then
    beta_log = gamma_log ( a ) + gamma_log ( b ) - gsumln ( a, b )
    return
  end if

  w = 0.0D+00

  if ( b < 8.0D+00 ) then
    go to 60
  end if

  beta_log = gamma_log ( a ) + algdiv ( a, b )
  return

40 continue
!
!  Reduction of A when 1000 < B.
!
  if ( 1000.0D+00 < b ) then

    n = a - 1.0D+00
    w = 1.0D+00
    do i = 1, n
      a = a - 1.0D+00
      w = w * ( a / ( 1.0D+00 + a / b ))
    end do

    beta_log = ( log ( w ) - n * log ( b ) ) &
      + ( gamma_log ( a ) + algdiv ( a, b ) )

    return
  end if

  n = a - 1.0D+00
  w = 1.0D+00
  do i = 1, n
    a = a - 1.0D+00
    h = a / b
    w = w * ( h / ( 1.0D+00 + h ) )
  end do
  w = log ( w )

  if ( 8.0D+00 <= b ) then
    beta_log = w + gamma_log ( a ) + algdiv ( a, b )
    return
  end if
!
!  Reduction of B when B < 8.
!
60 continue

  n = b - 1.0D+00
  z = 1.0D+00
  do i = 1, n
    b = b - 1.0D+00
    z = z * ( b / ( a + b ))
  end do

  beta_log = w + log ( z ) + ( gamma_log ( a ) + ( gamma_log ( b ) &
    - gsumln ( a, b ) ) )

  return
end
function beta_pser ( a, b, x, eps )

!*****************************************************************************80
!
!! BETA_PSER uses a power series expansion to evaluate IX(A,B)(X).
!
!  Discussion:
!
!    BETA_PSER is used when B <= 1 or B*X <= 0.7.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters.
!
!    Input, real (real64) X, the point where the function
!    is to be evaluated.
!
!    Input, real (real64) EPS, the tolerance.
!
!    Output, real (real64) BETA_PSER, the approximate value of IX(A,B)(X).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a0
  real (real64) algdiv
  real (real64) apb
  real (real64) b
  real (real64) b0
  real (real64) beta_log
  real (real64) beta_pser
  real (real64) c
  real (real64) eps
  real (real64) gam1
  real (real64) gamma_ln1
  integer i
  integer m
  real (real64) n
  real (real64) sum1
  real (real64) t
  real (real64) tol
  real (real64) u
  real (real64) w
  real (real64) x
  real (real64) z

  beta_pser = 0.0D+00

  if ( x == 0.0D+00 ) then
    return
  end if
!
!  Compute the factor X**A/(A*BETA(A,B))
!
  a0 = min ( a, b )

  if ( 1.0D+00 <= a0 ) then

    z = a * log ( x ) - beta_log ( a, b )
    beta_pser = exp ( z ) / a
  
  else

    b0 = max ( a, b )

    if ( b0 <= 1.0D+00 ) then

      beta_pser = x**a
      if ( beta_pser == 0.0D+00 ) then
        return
      end if

      apb = a + b

      if ( apb <= 1.0D+00 ) then
        z = 1.0D+00 + gam1 ( apb )
      else
        u = a + b - 1.0D+00
        z = ( 1.0D+00 + gam1 ( u ) ) / apb
      end if

      c = ( 1.0D+00 + gam1 ( a ) ) &
        * ( 1.0D+00 + gam1 ( b ) ) / z
      beta_pser = beta_pser * c * ( b / apb )

    else if ( b0 < 8.0D+00 ) then

      u = gamma_ln1 ( a0 )
      m = b0 - 1.0D+00

      c = 1.0D+00
      do i = 1, m
        b0 = b0 - 1.0D+00
        c = c * ( b0 / ( a0 + b0 ))
      end do

      u = log ( c ) + u
      z = a * log ( x ) - u
      b0 = b0 - 1.0D+00
      apb = a0 + b0

      if ( apb <= 1.0D+00 ) then
        t = 1.0D+00 + gam1 ( apb )
      else
        u = a0 + b0 - 1.0D+00
        t = ( 1.0D+00 + gam1 ( u ) ) / apb
      end if

      beta_pser = exp ( z ) * ( a0 / a ) &
        * ( 1.0D+00 + gam1 ( b0 )) / t

    else if ( 8.0D+00 <= b0 ) then

      u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
      z = a * log ( x ) - u
      beta_pser = ( a0 / a ) * exp ( z )

    end if

  end if

  if ( beta_pser == 0.0D+00 .or. a <= 0.1D+00 * eps ) then
    return
  end if
!
!  Compute the series.
!
  sum1 = 0.0D+00
  n = 0.0D+00
  c = 1.0D+00
  tol = eps / a

  do

    n = n + 1.0D+00
    c = c * ( 0.5D+00 + ( 0.5D+00 - b / n ) ) * x
    w = c / ( a + n )
    sum1 = sum1 + w

    if ( abs ( w ) <= tol ) then
      exit
    end if

  end do

  beta_pser = beta_pser * ( 1.0D+00 + a * sum1 )

  return
end
function beta_rcomp ( a, b, x, y )

!*****************************************************************************80
!
!! BETA_RCOMP evaluates X**A * Y**B / Beta(A,B).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters of the Beta function.
!    A and B should be nonnegative.
!
!    Input, real (real64) X, Y, define the numerator of the fraction.
!
!    Output, real (real64) BETA_RCOMP, the value of X**A * Y**B / Beta(A,B).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a0
  real (real64) algdiv
  real (real64) alnrel
  real (real64) apb
  real (real64) b
  real (real64) b0
  real (real64) bcorr
  real (real64) beta_log
  real (real64) beta_rcomp
  real (real64) c
  real (real64), parameter :: const = 0.398942280401433D+00
  real (real64) e
  real (real64) gam1
  real (real64) gamma_ln1
  real (real64) h
  integer i
  real (real64) lambda
  real (real64) lnx
  real (real64) lny
  integer n
  real (real64) rlog1
  real (real64) t
  real (real64) u
  real (real64) v
  real (real64) x
  real (real64) x0
  real (real64) y
  real (real64) y0
  real (real64) z

  beta_rcomp = 0.0D+00
  if ( x == 0.0D+00 .or. y == 0.0D+00 ) then
    return
  end if

  a0 = min ( a, b )

  if ( a0 < 8.0D+00 ) then

    if ( x <= 0.375D+00 ) then
      lnx = log ( x )
      lny = alnrel ( - x )      
    else if ( y <= 0.375D+00 ) then
      lnx = alnrel ( - y )
      lny = log ( y )
    else
      lnx = log ( x )
      lny = log ( y )
    end if

    z = a * lnx + b * lny

    if ( 1.0D+00 <= a0 ) then
      z = z - beta_log ( a, b )
      beta_rcomp = exp ( z )
      return
    end if
!
!  Procedure for A < 1 or B < 1
!
    b0 = max ( a, b )

    if ( b0 <= 1.0D+00 ) then

      beta_rcomp = exp ( z )
      if ( beta_rcomp == 0.0D+00 ) then
        return
      end if

      apb = a + b

      if ( apb <= 1.0D+00 ) then
        z = 1.0D+00 + gam1 ( apb )
      else
        u = a + b - 1.0D+00
        z = ( 1.0D+00 + gam1 ( u ) ) / apb
      end if

      c = ( 1.0D+00 + gam1 ( a ) ) &
        * ( 1.0D+00 + gam1 ( b ) ) / z
      beta_rcomp = beta_rcomp * ( a0 * c ) &
        / ( 1.0D+00 + a0 / b0 )

    else if ( b0 < 8.0D+00 ) then

      u = gamma_ln1 ( a0 )
      n = b0 - 1.0D+00

      c = 1.0D+00
      do i = 1, n
        b0 = b0 - 1.0D+00
        c = c * ( b0 / ( a0 + b0 ))
      end do
      u = log ( c ) + u

      z = z - u
      b0 = b0 - 1.0D+00
      apb = a0 + b0

      if ( apb <= 1.0D+00 ) then
        t = 1.0D+00 + gam1 ( apb )
      else
        u = a0 + b0 - 1.0D+00
        t = ( 1.0D+00 + gam1 ( u ) ) / apb
      end if

      beta_rcomp = a0 * exp ( z ) * ( 1.0D+00 + gam1 ( b0 ) ) / t

    else if ( 8.0D+00 <= b0 ) then

      u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
      beta_rcomp = a0 * exp ( z - u ) 

    end if

  else

    if ( a <= b ) then
      h = a / b
      x0 = h / ( 1.0D+00 + h )
      y0 = 1.0D+00 / (  1.0D+00 + h )
      lambda = a - ( a + b ) * x
    else
      h = b / a
      x0 = 1.0D+00 / ( 1.0D+00 + h )
      y0 = h / ( 1.0D+00 + h )
      lambda = ( a + b ) * y - b
    end if

    e = -lambda / a

    if ( abs ( e ) <= 0.6D+00 ) then
      u = rlog1 ( e )
    else
      u = e - log ( x / x0 )
    end if

    e = lambda / b
 
    if ( abs ( e ) <= 0.6D+00 ) then
      v = rlog1 ( e )
    else
      v = e - log ( y / y0 )
    end if

    z = exp ( - ( a * u + b * v ) )
    beta_rcomp = const * sqrt ( b * x0 ) * z * exp ( - bcorr ( a, b ))

  end if

  return
end
function beta_rcomp1 ( mu, a, b, x, y )

!*****************************************************************************80
!
!! BETA_RCOMP1 evaluates exp(MU) * X**A * Y**B / Beta(A,B).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, integer MU, ?
!
!    Input, real (real64) A, B, the parameters of the Beta function.
!    A and B should be nonnegative.
!
!    Input, real (real64) X, Y, quantities whose powers form part of
!    the expression.
!
!    Output, real (real64) BETA_RCOMP1, the value of
!    exp(MU) * X**A * Y**B / Beta(A,B).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a0
  real (real64) algdiv
  real (real64) alnrel
  real (real64) apb
  real (real64) b
  real (real64) b0
  real (real64) bcorr
  real (real64) beta_log
  real (real64) beta_rcomp1
  real (real64) c
  real (real64), parameter :: const = 0.398942280401433D+00
  real (real64) e
  real (real64) esum
  real (real64) gam1
  real (real64) gamma_ln1
  real (real64) h
  integer i
  real (real64) lambda
  real (real64) lnx
  real (real64) lny
  integer mu
  integer n
  real (real64) rlog1
  real (real64) t
  real (real64) u
  real (real64) v
  real (real64) x
  real (real64) x0
  real (real64) y
  real (real64) y0
  real (real64) z

  a0 = min ( a, b )
!
!  Procedure for 8 <= A and 8 <= B.
!
  if ( 8.0D+00 <= a0 ) then

    if ( a <= b ) then
      h = a / b
      x0 = h / ( 1.0D+00 + h )
      y0 = 1.0D+00 / ( 1.0D+00 + h )
      lambda = a - ( a + b ) * x
    else
      h = b / a
      x0 = 1.0D+00 / ( 1.0D+00 + h )
      y0 = h / ( 1.0D+00 + h )
      lambda = ( a + b ) * y - b
    end if

    e = -lambda / a

    if ( abs ( e ) <= 0.6D+00 ) then
      u = rlog1 ( e )
    else
      u = e - log ( x / x0 )
    end if

    e = lambda / b

    if ( abs ( e ) <= 0.6D+00 ) then
      v = rlog1 ( e )
    else
      v = e - log ( y / y0 )
    end if

    z = esum ( mu, - ( a * u + b * v ))
    beta_rcomp1 = const * sqrt ( b * x0 ) * z * exp ( - bcorr ( a, b ) )
!
!  Procedure for A < 8 or B < 8.
!
  else

    if ( x <= 0.375D+00 ) then
      lnx = log ( x )
      lny = alnrel ( - x )
    else if ( y <= 0.375D+00 ) then
      lnx = alnrel ( - y )
      lny = log ( y )
    else
      lnx = log ( x )
      lny = log ( y )
    end if
  
    z = a * lnx + b * lny

    if ( 1.0D+00 <= a0 ) then
      z = z - beta_log ( a, b )
      beta_rcomp1 = esum ( mu, z )
      return
    end if
!
!  Procedure for A < 1 or B < 1.
!
    b0 = max ( a, b )

    if ( 8.0D+00 <= b0 ) then
      u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
      beta_rcomp1 = a0 * esum ( mu, z-u )
      return
    end if

    if ( 1.0D+00 < b0 ) then
!
!  Algorithm for 1 < B0 < 8
!
      u = gamma_ln1 ( a0 )
      n = b0 - 1.0D+00

      c = 1.0D+00
      do i = 1, n
        b0 = b0 - 1.0D+00
        c = c * ( b0 / ( a0 + b0 ) )
      end do
      u = log ( c ) + u

      z = z - u
      b0 = b0 - 1.0D+00
      apb = a0 + b0

      if ( apb <= 1.0D+00 ) then
        t = 1.0D+00 + gam1 ( apb )
      else
        u = a0 + b0 - 1.0D+00
        t = ( 1.0D+00 + gam1 ( u ) ) / apb
      end if

      beta_rcomp1 = a0 * esum ( mu, z ) &
        * ( 1.0D+00 + gam1 ( b0 ) ) / t
!
!  Algorithm for B0 <= 1
!
    else

      beta_rcomp1 = esum ( mu, z )
      if ( beta_rcomp1 == 0.0D+00 ) then
        return
      end if

      apb = a + b

      if ( apb <= 1.0D+00 ) then
        z = 1.0D+00 + gam1 ( apb )
      else
        u = real ( a, kind = real64 ) + real ( b, kind = real64 ) - 1.0D+00
        z = ( 1.0D+00 + gam1 ( u )) / apb
      end if

      c = ( 1.0D+00 + gam1 ( a ) ) &
        * ( 1.0D+00 + gam1 ( b ) ) / z
      beta_rcomp1 = beta_rcomp1 * ( a0 * c ) / ( 1.0D+00 + a0 / b0 )

    end if

  end if

  return
end
function beta_up ( a, b, x, y, n, eps )

!*****************************************************************************80
!
!! BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, the parameters of the function.
!    A and B should be nonnegative.
!
!    Input, real (real64) X, Y, ?
!
!    Input, integer N, the increment to the first argument of IX.
!
!    Input, real (real64) EPS, the tolerance.
!
!    Output, real (real64) BETA_UP, the value of IX(A,B) - IX(A+N,B).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) ap1
  real (real64) apb
  real (real64) b
  real (real64) beta_rcomp1
  real (real64) beta_up
  real (real64) d
  real (real64) eps
  real (real64) exparg
  integer i
  integer k
  real (real64) l
  integer mu
  integer n
  real (real64) r
  real (real64) t
  real (real64) w
  real (real64) x
  real (real64) y
!
!  Obtain the scaling factor EXP(-MU) AND
!  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
!
  apb = a + b
  ap1 = a + 1.0D+00
  mu = 0
  d = 1.0D+00

  if ( n /= 1 ) then

    if ( 1.0D+00 <= a ) then

      if ( 1.1D+00 * ap1 <= apb ) then
        mu = abs ( exparg ( 1 ) )
        k = exparg ( 0 )
        if ( k < mu ) then
          mu = k
        end if
        t = mu
        d = exp ( - t )
      end if

    end if

  end if

  beta_up = beta_rcomp1 ( mu, a, b, x, y ) / a

  if ( n == 1 .or. beta_up == 0.0D+00 ) then
    return
  end if

  w = d
!
!  Let K be the index of the maximum term.
!
  k = 0

  if ( 1.0D+00 < b ) then

    if ( y <= 0.0001D+00 ) then

      k = n - 1

    else

      r = ( b - 1.0D+00 ) * x / y - a

      if ( 1.0D+00 <= r ) then
        k = n - 1
        t = n - 1
        if ( r < t ) then
          k = r
        end if
      end if

    end if
!
!  Add the increasing terms of the series.
!
    do i = 1, k
      l = i - 1
      d = ( ( apb + l ) / ( ap1 + l ) ) * x * d
      w = w + d
    end do

  end if
!
!  Add the remaining terms of the series.
!
  do i = k+1, n-1
    l = i - 1
    d = ( ( apb + l ) / ( ap1 + l ) ) * x * d
    w = w + d
    if ( d <= eps * w ) then
      beta_up = beta_up * w
      return
    end if
  end do

  beta_up = beta_up * w

  return
end
subroutine cdfbet ( which, p, q, x, y, a, b, status, bound )

!*****************************************************************************80
!
!! CDFBET evaluates the CDF of the Beta Distribution.
!
!  Discussion:
!
!    This routine calculates any one parameter of the beta distribution 
!    given the others.
!
!    The value P of the cumulative distribution function is calculated 
!    directly by code associated with the reference.
!
!    Computation of the other parameters involves a seach for a value that
!    produces the desired value of P.  The search relies on the
!    monotonicity of P with respect to the other parameters.
!
!    The beta density is proportional to t^(A-1) * (1-t)^(B-1).
!
!  Modified:
!
!    14 April 2007
!
!    Modified February 2024 (error handling):
!    Seth Siefken
!    U.S. Geological Survey
!    Helena, Montana
!
!  Author:
!
!    Barry Brown, James Lovato, Kathy Russell
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software, 
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, integer WHICH, indicates which of the next four argument
!    values is to be calculated from the others.
!    1: Calculate P and Q from X, Y, A and B;
!    2: Calculate X and Y from P, Q, A and B;
!    3: Calculate A from P, Q, X, Y and B;
!    4: Calculate B from P, Q, X, Y and A.
!
!    Input/output, real (real64) P, the integral from 0 to X of the
!    chi-square distribution.  Input range: [0, 1].
!
!    Input/output, real (real64) Q, equals 1-P.  Input range: [0, 1].
!
!    Input/output, real (real64) X, the upper limit of integration 
!    of the beta density.  If it is an input value, it should lie in
!    the range [0,1].  If it is an output value, it will be searched for
!    in the range [0,1].
!
!    Input/output, real (real64) Y, equal to 1-X.  If it is an input
!    value, it should lie in the range [0,1].  If it is an output value,
!    it will be searched for in the range [0,1].
!
!    Input/output, real (real64) A, the first parameter of the beta
!    density.  If it is an input value, it should lie in the range
!    (0, +infinity).  If it is an output value, it will be searched
!    for in the range [1D-300,1D300].
!
!    Input/output, real (real64) B, the second parameter of the beta
!    density.  If it is an input value, it should lie in the range
!    (0, +infinity).  If it is an output value, it will be searched
!    for in the range [1D-300,1D300].
!
!    Output, integer STATUS, reports the status of the computation.
!     0, if the calculation completed correctly;
!    -I, if the input parameter number I is out of range;
!    +1, if the answer appears to be lower than lowest search bound;
!    +2, if the answer appears to be higher than greatest search bound;
!    +3, if P + Q /= 1;
!    +4, if X + Y /= 1.
!
!    Output, real (real64) BOUND, is only defined if STATUS is nonzero.
!    If STATUS is negative, then this is the value exceeded by parameter I.
!    if STATUS is 1 or 2, this is the search bound that was exceeded.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64), parameter :: atol = 1.0D-10
  real (real64) b
  real (real64) bound
  real (real64) ccum
  real (real64) cum
  real (real64) fx
  real (real64), parameter :: inf = 1.0D+300
  real (real64) p
  real (real64) q
  logical qhi
  logical qleft
  integer status
  real (real64), parameter :: tol = 1.0D-08
  integer which
  real (real64) x
  real (real64) xhi
  real (real64) xlo
  real (real64) y

  status = 0
  bound = 0.0D+00

  if ( which < 1 ) then
    bound = 1.0D+00
    status = -1
    error stop 'CDFBET - Input WHICH outside legal range between 1 and 4.'
  end if

  if ( 4 < which ) then
    bound = 4.0D+00
    status = -1
    error stop 'CDFBET - Input WHICH outside legal range between 1 and 4.'
  end if
!
!  Unless P is to be computed, make sure it is legal.
!
  if ( which /= 1 ) then
    if ( p < 0.0D+00 ) then
      bound = 0.0D+00
      status = -2
      error stop 'CDFBET - Input parameter P is out of range.'
    else if ( 1.0D+00 < p ) then
      bound = 1.0D+00
      status = -2
      error stop 'CDFBET - Input parameter P is out of range.'
    end if
  end if
!
!  Unless Q is to be computed, make sure it is legal.
!
  if ( which /= 1 ) then
    if ( q < 0.0D+00 ) then
      bound = 0.0D+00
      status = -3
      error stop 'CDFBET - Input parameter Q is out of range.'
    else if ( 1.0D+00 < q ) then
      bound = 1.0D+00
      status = -3
      error stop 'CDFBET - Input parameter Q is out of range.'
    end if
  end if
!
!  Unless X is to be computed, make sure it is legal.
!
  if ( which /= 2 ) then
    if ( x < 0.0D+00 ) then
      bound = 0.0D+00
      status = -4
      error stop 'CDFBET - Input parameter X is out of range.'
    else if ( 1.0D+00 < x ) then
      bound = 1.0D+00
      status = -4
      error stop 'CDFBET - Input parameter X is out of range.'
    end if
  end if
!
!  Unless Y is to be computed, make sure it is legal.
!
  if ( which /= 2 ) then
    if ( y < 0.0D+00 ) then
      bound = 0.0D+00
      status = -5
      error stop 'CDFBET - Input parameter Y is out of range.'
    else if ( 1.0D+00 < y ) then
      bound = 1.0D+00
      status = -5
      error stop 'CDFBET - Input parameter Y is out of range.'
    end if
  end if
!
!  Unless A is to be computed, make sure it is legal.
!
  if ( which /= 3 ) then
    if ( a <= 0.0D+00 ) then
      bound = 0.0D+00
      status = -6
      error stop 'CDFBET - Input parameter A is out of range.'
    end if
  end if
!
!  Unless B is to be computed, make sure it is legal.
!
  if ( which /= 4 ) then
    if ( b <= 0.0D+00 ) then
      bound = 0.0D+00
      status = -7
      error stop 'CDFBET - Input parameter A is out of range.'
    end if
  end if
!
!  Check that P + Q = 1.
!
  if ( which /= 1 ) then
    if ( 3.0D+00 * epsilon ( p ) < abs ( ( p + q ) - 1.0D+00 ) ) then
      status = 3
      error stop 'CDFBET - Fatal error: P + Q /= 1'
    end if
  end if
!
!  Check that X + Y = 1.
!
  if ( which /= 2 ) then
    if ( 3.0D+00 * epsilon ( x ) < abs ( ( x + y ) - 1.0D+00 ) ) then
      status = 4
      error stop 'CDFBET - Fatal error: X + Y /= 1'
    end if
  end if
!
!  Compute P and Q.
!
  if ( which == 1 ) then

    call cumbet ( x, y, a, b, p, q )
    status = 0
!
!  Compute X and Y.
!
  else if ( which == 2 ) then

    call dstzr ( 0.0D+00, 1.0D+00, atol, tol )

    if ( p <= q ) then

      status = 0
      fx = 0.0D+00
      call dzror ( status, x, fx, xlo, xhi, qleft, qhi )
      y = 1.0D+00 - x

      do while ( status == 1 )
 
        call cumbet ( x, y, a, b, cum, ccum )
        fx = cum - p
        call dzror ( status, x, fx, xlo, xhi, qleft, qhi )
        y = 1.0D+00 - x
 
      end do

    else

      status = 0
      fx = 0.0D+00
      call dzror ( status, y, fx, xlo, xhi, qleft, qhi )
      x = 1.0D+00 - y

      do while ( status == 1 )

        call cumbet ( x, y, a, b, cum, ccum )
        fx = ccum - q
        call dzror ( status, y, fx, xlo, xhi, qleft, qhi )
        x = 1.0D+00 - y

      end do

    end if

    if ( status == -1 ) then
      if ( qleft ) then
        status = 1
        bound = 0.0D+00
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'CDFBET - Warning!'
!       write ( *, '(a)' ) '  The desired answer appears to be lower than'
!       write ( *, '(a,g14.6)' ) '  the search bound of ', bound
      else
        status = 2
        bound = 1.0D+00
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'CDFBET - Warning!'
!       write ( *, '(a)' ) '  The desired answer appears to be higher than'
!       write ( *, '(a,g14.6)' ) '  the search bound of ', bound
      end if
    end if
!
!  Compute A.
!
  else if ( which == 3 ) then

    call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )

    status = 0
    a = 5.0D+00
    fx = 0.0D+00

    call dinvr ( status, a, fx, qleft, qhi )

    do while ( status == 1 )

      call cumbet ( x, y, a, b, cum, ccum )

      if ( p <= q ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, a, fx, qleft, qhi )

    end do

    if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0D+00
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'CDFBET - Warning!'
!       write ( *, '(a)' ) '  The desired answer appears to be lower than'
!       write ( *, '(a,g14.6)' ) '  the search bound of ', bound
      else
        status = 2
        bound = inf
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'CDFBET - Warning!'
!       write ( *, '(a)' ) '  The desired answer appears to be higher than'
!       write ( *, '(a,g14.6)' ) '  the search bound of ', bound
      end if

    end if
!
!  Compute B.
!
  else if ( which == 4 ) then

    call dstinv ( 0.0D+00, inf, 0.5D+00, 0.5D+00, 5.0D+00, atol, tol )

    status = 0
    b = 5.0D+00
    fx = 0.0D+00

    call dinvr ( status, b, fx, qleft, qhi )

    do while ( status == 1 )

      call cumbet ( x, y, a, b, cum, ccum )

      if ( p <= q ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, b, fx, qleft, qhi )

    end do

    if ( status == -1 ) then
      if ( qleft ) then
        status = 1
        bound = 0.0D+00
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'CDFBET - Warning!'
!       write ( *, '(a)' ) '  The desired answer appears to be lower than'
!       write ( *, '(a,g14.6)' ) '  the search bound of ', bound
      else
        status = 2
        bound = inf
!       write ( *, '(a)' ) ' '
!       write ( *, '(a)' ) 'CDFBET - Warning!'
!       write ( *, '(a)' ) '  The desired answer appears to be higher than'
!       write ( *, '(a,g14.6)' ) '  the search bound of ', bound
      end if
    end if

  end if

  return
end

subroutine cumbet ( x, y, a, b, cum, ccum )

!*****************************************************************************80
!
!! CUMBET evaluates the cumulative incomplete beta distribution.
!
!  Discussion:
!
!    This routine calculates the CDF to X of the incomplete beta distribution
!    with parameters A and B.  This is the integral from 0 to x
!    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
!
!  Author:
!
!    Barry Brown, James Lovato, Kathy Russell
!
!  Parameters:
!
!    Input, real (real64) X, the upper limit of integration.
!
!    Input, real (real64) Y, the value of 1-X.
!
!    Input, real (real64) A, B, the parameters of the distribution.
!
!    Output, real (real64) CUM, CCUM, the values of the cumulative
!    density function and complementary cumulative density function.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) b
  real (real64) ccum
  real (real64) cum
  integer ierr
  real (real64) x
  real (real64) y

  if ( x <= 0.0D+00 ) then

    cum = 0.0
    ccum = 1.0D+00

  else if ( y <= 0.0D+00 ) then

    cum = 1.0D+00
    ccum = 0.0

  else

    call beta_inc ( a, b, x, y, cum, ccum, ierr )

  end if

  return
end
subroutine dinvr ( status, x, fx, qleft, qhi )

!*****************************************************************************80
!
!! DINVR bounds the zero of the function and invokes DZROR.
!
!  Discussion:
!
!    This routine seeks to find bounds on a root of the function and 
!    invokes DZROR to perform the zero finding.  DSTINV must have been 
!    called before this routine in order to set its parameters.
!
!  Reference:
!
!    JCP Bus, TJ Dekker,
!    Two Efficient Algorithms with Guaranteed Convergence for 
!    Finding a Zero of a Function,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 4, pages 330-345, 1975.
!
!  Parameters:
!
!    Input/output, integer STATUS.  At the beginning of a zero finding 
!    problem, STATUS should be set to 0 and this routine invoked.  The value
!    of parameters other than X will be ignored on this call.
!    If this routine needs the function to be evaluated, it will set STATUS 
!    to 1 and return.  The value of the function should be set in FX and 
!    this routine again called without changing any of its other parameters.
!    If this routine finishes without error, it returns with STATUS 0, 
!    and X an approximate root of F(X).
!    If this routine cannot bound the function, it returns a negative STATUS and 
!    sets QLEFT and QHI.
!
!    Output, real (real64) X, the value at which F(X) is to be evaluated.
!
!    Input, real (real64) FX, the value of F(X) calculated by the user
!    on the previous call, when this routine returned with STATUS = 1.
!
!    Output, logical QLEFT, is defined only if QMFINV returns FALSE.  In that
!    case, QLEFT is TRUE if the stepping search terminated unsucessfully 
!    at SMALL, and FALSE if the search terminated unsucessfully at BIG.
!
!    Output, logical QHI, is defined only if QMFINV returns FALSE.  In that
!    case, it is TRUE if Y < F(X) at the termination of the search and FALSE
!    if F(X) < Y.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) :: absstp
  real (real64) :: abstol
  real (real64) :: big
  real (real64) fbig
  real (real64) fsmall
  real (real64) fx
  integer i99999
  logical qbdd
  logical qcond
  logical qdum1
  logical qdum2
  logical qhi
  logical qincr
  logical qleft
  logical qlim
  logical qup
  real (real64) :: relstp
  real (real64) :: reltol
  real (real64) :: small
  integer status
  real (real64) step
  real (real64) :: stpmul
  real (real64) x
  real (real64) xhi
  real (real64) xlb
  real (real64) xlo
  real (real64) xsave
  real (real64) xub
  real (real64) yy
  real (real64) zabsst
  real (real64) zabsto
  real (real64) zbig
  real (real64) zrelst
  real (real64) zrelto
  real (real64) zsmall
  real (real64) zstpmu

  save

  if ( 0 < status ) then
!    go to i99999
     if ( i99999 == 10 ) go to 10
     if ( i99999 == 20 ) go to 20
     if ( i99999 == 90 ) go to 90
     if ( i99999 == 130 ) go to 130
     if ( i99999 == 200 ) go to 200
     if ( i99999 == 270 ) go to 270
  end if

  qcond = .not. ( small <= x .and. x <= big )

  if ( .not. ( small <= x .and. x <= big ) ) then
    error stop 'DINVR - The values SMALL, X, BIG are not monotone.'
  end if

  xsave = x
!
!  See that SMALL and BIG bound the zero and set QINCR.
!
  x = small
!
!  GET-function-VALUE
!
!  assign 10 to i99999
  i99999 = 10
  status = 1
  return

   10 continue

  fsmall = fx
  x = big
!
!  GET-function-VALUE
!
!  assign 20 to i99999
  i99999 = 20
  status = 1
  return

   20 continue

  fbig = fx

  qincr = ( fsmall < fbig )

  if ( fsmall <= fbig ) then

    if ( 0.0D+00 < fsmall ) then
      status = -1
      qleft = .true.
      qhi = .true.
      return
    end if

    if ( fbig < 0.0D+00 ) then
      status = -1
      qleft = .false.
      qhi = .false.
      return
    end if

  else if ( fbig < fsmall ) then

    if ( fsmall < 0.0D+00 ) then
      status = -1
      qleft = .true.
      qhi = .false.
      return
    end if

    if ( 0.0D+00 < fbig ) then
      status = -1
      qleft = .false.
      qhi = .true.
      return
    end if

  end if

  x = xsave
  step = max ( absstp, relstp * abs ( x ) )
!
!  YY = F(X) - Y
!  GET-function-VALUE
!
!  assign 90 to i99999
  i99999 = 90
  status = 1
  return

   90 continue

  yy = fx

  if ( yy == 0.0D+00 ) then
    status = 0
    return
  end if

  100 continue

  qup = ( qincr .and. ( yy < 0.0D+00 ) ) .or. &
        ( .not. qincr .and. ( 0.0D+00 < yy ) )
!
!  Handle case in which we must step higher.
!
  if (.not. qup ) then
    go to 170
  end if

  xlb = xsave
  xub = min ( xlb + step, big )
  go to 120

  110 continue

  if ( qcond ) then
    go to 150
  end if
!
!  YY = F(XUB) - Y
!
  120 continue

  x = xub
!
!  GET-function-VALUE
!
!  assign 130 to i99999
  i99999 = 130
  status = 1
  return

  130 continue

  yy = fx
  qbdd = ( qincr .and. ( 0.0D+00 <= yy ) ) .or. &
    ( .not. qincr .and. ( yy <= 0.0D+00 ) )
  qlim = ( big <= xub )
  qcond = qbdd .or. qlim

  if ( .not. qcond ) then
    step = stpmul * step
    xlb = xub
    xub = min ( xlb + step, big )
  end if

  go to 110

  150 continue

  if ( qlim .and. .not. qbdd ) then
    status = -1
    qleft = .false.
    qhi = .not. qincr
    x = big
    return
  end if

  160 continue

  go to 240
!
!  Handle the case in which we must step lower.
!
  170 continue

  xub = xsave
  xlb = max ( xub - step, small )
  go to 190

  180 continue

  if ( qcond ) then
    go to 220
  end if
!
!  YY = F(XLB) - Y
!
  190 continue

  x = xlb
!
!  GET-function-VALUE
!
!  assign 200 to i99999
  i99999 = 200
  status = 1
  return

  200 continue

  yy = fx
  qbdd = ( qincr .and. ( yy <= 0.0D+00 ) ) .or. &
    ( .not. qincr .and. ( 0.0D+00 <= yy ) )
  qlim = xlb <= small
  qcond = qbdd .or. qlim

  if ( .not. qcond ) then
    step = stpmul * step
    xub = xlb
    xlb = max ( xub - step, small )
  end if

  go to 180

  220 continue

  if ( qlim .and. ( .not. qbdd ) ) then
    status = -1
    qleft = .true.
    qhi = qincr
    x = small
    return
  end if

  230 continue
  240 continue

  call dstzr ( xlb, xub, abstol, reltol )
!
!  If we reach here, XLB and XUB bound the zero of F.
!
  status = 0
  go to 260

  250 continue

    if ( status /= 1 ) then
      x = xlo
      status = 0
      return
    end if

  260 continue

  call dzror ( status, x, fx, xlo, xhi, qdum1, qdum2 )

  if ( status /= 1 ) then
    go to 250
  end if
!
!  GET-function-VALUE
!
!  assign 270 to i99999
  i99999 = 270
  status = 1
  return

  270 continue
  go to 250

entry dstinv ( zsmall, zbig, zabsst, zrelst, zstpmu, zabsto, zrelto )

!*****************************************************************************80
!
!! DSTINV SeT INverse finder - Reverse Communication
!
!  Discussion:
!
!    This routine is given a monotone function F, and a value Y, 
!    and seeks an argument value X such that F(X) = Y.  
!
!    This routine uses reverse communication -- see DINVR.
!    This routine sets quantities needed by DINVR.
!
!    F must be a monotone function, the results of QMFINV are
!    otherwise undefined.  QINCR must be TRUE if F is nondecreasing 
!    and FALSE if F is nonincreasing.
!
!    QMFINV will return TRUE if and only if F(SMALL) and
!    F(BIG) bracket Y, i. e.,
!      QINCR is TRUE and F(SMALL) <= Y <= F(BIG) or
!      QINCR is FALSE and F(BIG) <= Y <= F(SMALL)
!
!    If QMFINV returns TRUE, then the X returned satisfies
!    the following condition.  Let
!      TOL(X) = MAX ( ABSTOL, RELTOL * ABS ( X ) )
!    then if QINCR is TRUE,
!      F(X-TOL(X)) <= Y <= F(X+TOL(X))
!    and if QINCR is FALSE
!      F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
!
!    Compares F(X) with Y for the input value of X then uses QINCR
!    to determine whether to step left or right to bound the
!    desired X.  The initial step size is
!
!      max ( ABSSTP, RELSTP * ABS ( S ) ) 
!
!    for the input value of X.
!
!    Iteratively steps right or left until it bounds X.
!    At each step which doesn't bound X, the step size is doubled.
!    The routine is careful never to step beyond SMALL or BIG.  If
!    it hasn't bounded X at SMALL or BIG, QMFINV returns FALSE
!    after setting QLEFT and QHI.
!
!    If X is successfully bounded then Algorithm R of the paper
!    Bus and Dekker is employed to find the zero of the function F(X)-Y. 
!    This is routine QRZERO.
!
!  Reference:
!
!    JCP Bus, TJ Dekker,
!    Two Efficient Algorithms with Guaranteed Convergence for 
!    Finding a Zero of a Function,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 4, pages 330-345, 1975.
!
!  Parameters:
!
!    Input, real (real64) ZSMALL, ZBIG, the left and right endpoints 
!    of the interval to be searched for a solution.
!
!    Input, real (real64) ZABSST, ZRELSTP, the initial step size in 
!    the search is max ( ZABSST, ZRELST * abs ( X ) ).
!
!    Input, real (real64) STPMUL.  When a step doesn't bound the zero, 
!    the stepsize is multiplied by STPMUL and another step taken.  A 
!    popular value is 2.0.
!
!    Input, real (real64) ABSTOL, RELTOL, two numbers that determine 
!    the accuracy of the solution
!
  small = zsmall
  big = zbig
  absstp = zabsst
  relstp = zrelst
  stpmul = zstpmu
  abstol = zabsto
  reltol = zrelto

  return
end
subroutine dzror ( status, x, fx, xlo, xhi, qleft, qhi )

!*****************************************************************************80
!
!! DZROR seeks a zero of a function, using reverse communication.
!
!  Discussion:
!
!    This routine performs the zero finding.  STZROR must have been called
!    before this routine in order to set its parameters.
!
!  Modified:
!
!    09 June 2004
!
!  Reference:
!
!    JCP Bus, TJ Dekker,
!    Two Efficient Algorithms with Guaranteed Convergence for 
!    Finding a Zero of a Function,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 4, pages 330-345, 1975.
!
!  Parameters:
!
!    Input/output, integer STATUS.  At the beginning of a zero finding problem,
!    STATUS should be set to 0 and ZROR invoked.  The value of other 
!    parameters will be ignored on this call.
!    When ZROR needs the function evaluated, it will set
!    STATUS to 1 and return.  The value of the function
!    should be set in FX and ZROR again called without
!    changing any of its other parameters.
!    When ZROR has finished without error, it will return
!    with STATUS 0.  In that case (XLO,XHI) bound the answe
!    If ZROR finds an error (which implies that F(XLO)-Y an
!    F(XHI)-Y have the same sign, it returns STATUS -1.  In
!    this case, XLO and XHI are undefined.
!
!    Output, real (real64) X, the value of X at which F(X) is to 
!    be evaluated.
!
!    Input, real (real64) FX, the value of F(X), which must be calculated 
!    by the user when ZROR has returned on the previous call with STATUS = 1.
!
!    Output, real (real64) XLO, XHI, are lower and upper bounds for the
!    solution when ZROR returns with STATUS = 0.
!
!    Output, logical QLEFT,is TRUE if the stepping search terminated
!    unsucessfully at XLO.  If it is FALSE, the search terminated
!    unsucessfully at XHI.
!
!    Output, logical QHI, is TRUE if Y < F(X) at the termination of the
!    search and FALSE if F(X) < Y at the termination of the search.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) abstol
  real (real64) b
  real (real64) c
  real (real64) d
  integer ext
  real (real64) fa
  real (real64) fb
  real (real64) fc
  real (real64) fd
  real (real64) fda
  real (real64) fdb
  logical first
  real (real64) ftol
  real (real64) fx
  integer i99999
  real (real64) m
  real (real64) mb
  real (real64) p
  real (real64) q
  logical qhi
  logical qleft
  logical qrzero
  real (real64) reltol
  integer status
  real (real64) tol
  real (real64) w
  real (real64) x
  real (real64) xhi
  real (real64) xlo
  real (real64) :: xxhi = 0.0D+00
  real (real64) :: xxlo = 0.0D+00
  real (real64) zabstl
  real (real64) zreltl
  real (real64) zx
  real (real64) zxhi
  real (real64) zxlo

  save

  ftol(zx) = 0.5D+00 * max ( abstol, reltol * abs ( zx ) )

  if ( 0 < status ) then
    go to 280
  end if

  xlo = xxlo
  xhi = xxhi
  b = xlo
  x = xlo
!
!     GET-function-VALUE
!
!  assign 10 to i99999
  i99999 = 10
  go to 270

10 continue

  fb = fx
  xlo = xhi
  a = xlo
  x = xlo
!
!     GET-function-VALUE
!
!  assign 20 to i99999
  i99999 = 20
  go to 270
!
!  Check that F(ZXLO) < 0 < F(ZXHI)  or F(ZXLO) > 0 > F(ZXHI)
!
20 continue

  if ( fb < 0.0D+00 ) then
    if ( fx < 0.0D+00 ) then
      status = -1
      qleft = ( fx < fb )
      qhi = .false.
      return
    end if
  end if

  if ( 0.0D+00 < fb ) then
    if ( 0.0D+00 < fx ) then
      status = -1
      qleft = ( fb < fx )
      qhi = .true.
      return
    end if
  end if

  fa = fx
  first = .true.

70 continue

  c = a
  fc = fa
  ext = 0

80 continue

  if ( abs ( fc ) < abs ( fb ) ) then

    if ( c == a ) then
      d = a
      fd = fa
    end if

    a = b
    fa = fb
    xlo = c
    b = xlo
    fb = fc
    c = a
    fc = fa

  end if

  tol = ftol ( xlo )
  m = ( c + b ) * 0.5D+00
  mb = m - b

  if (.not. ( tol < abs ( mb ) ) ) then
    go to 240
  end if

  if ( 3 < ext ) then
    w = mb
    go to 190
  end if

  110 continue

  tol = sign ( tol, mb )
  p = ( b - a ) * fb
!
!  I had to insert a rudimentary check on the divisions here
!  to avoid ninny errors, JVB, 09 June 2004.
!
  if ( first ) then

    q = fa - fb
    first = .false.

  else

    if ( d == b ) then
      fdb = 1.0D+00
    else
      fdb = ( fd - fb ) / ( d - b )
    end if

    if ( d == a ) then
      fda = 1.0D+00
    else
      fda = ( fd - fa ) / ( d - a )
    end if

    p = fda * p
    q = fdb * fa - fda * fb

  end if

  130 continue

  if ( p < 0.0D+00 ) then
    p = -p
    q = -q
  end if

  140 continue

  if ( ext == 3 ) then
    p = p *  2.0D+00 
  end if

  if (.not. ( ( p * 1.0D+00 ) == 0.0D+00 .or. p <= ( q * tol ) ) ) then
    go to 150
  end if

  w = tol
  go to 180

  150 continue

  if ( p < mb * q ) then
    w = p / q
  else
    w = mb
  end if

  180 continue
  190 continue

  d = a
  fd = fa
  a = b
  fa = fb
  b = b + w
  xlo = b
  x = xlo
!
!  GET-function-VALUE
!
!  assign 200 to i99999
  i99999 = 200
  go to 270

  200 continue

  fb = fx

  if ( 0.0D+00 <= fc * fb ) then

    go to 70

  else

    if ( w == mb ) then
      ext = 0
    else
      ext = ext + 1
    end if

    go to 80

  end if

  240 continue

  xhi = c
  qrzero = ( 0.0D+00 <= fc .and. fb <= 0.0D+00 ) .or. &
    ( fc < 0.0D+00 .and. fb >= 0.0D+00 )

  if ( qrzero ) then
    status = 0
  else
    status = -1
  end if

  return

entry dstzr ( zxlo, zxhi, zabstl, zreltl )

!*****************************************************************************80
!
!! DSTZR - SeT ZeRo finder - Reverse communication version
!
!  Discussion:
!
!    This routine sets quantities needed by ZROR.  The function of ZROR
!    and the quantities set is given here.
!
!    Given a function F, find XLO such that F(XLO) = 0.
!
!     Input condition. F is a real (real64) function of a single
!     real (real64) argument and XLO and XHI are such that
!          F(XLO)*F(XHI)  <=  0.0
!
!     If the input condition is met, QRZERO returns .TRUE.
!     and output values of XLO and XHI satisfy the following
!          F(XLO)*F(XHI)  <= 0.
!          ABS ( F(XLO) ) <= ABS ( F(XHI) )
!          ABS ( XLO - XHI ) <= TOL(X)
!     where
!          TOL(X) = MAX ( ABSTOL, RELTOL * ABS ( X ) )
!
!     If this algorithm does not find XLO and XHI satisfying
!     these conditions then QRZERO returns .FALSE.  This
!     implies that the input condition was not met.
!
!  Reference:
!
!    JCP Bus, TJ Dekker,
!    Two Efficient Algorithms with Guaranteed Convergence for 
!    Finding a Zero of a Function,
!    ACM Transactions on Mathematical Software,
!    Volume 1, Number 4, pages 330-345, 1975.
!
!  Parameters:
!
!    Input, real (real64) XLO, XHI, the left and right endpoints of the 
!    interval to be searched for a solution.
!
!    Input, real (real64) ABSTOL, RELTOL, two numbers that determine 
!    the accuracy of the solution.
!
  xxlo = zxlo
  xxhi = zxhi
  abstol = zabstl
  reltol = zreltl
  return
!
!     TO GET-function-VALUE
!
  270 status = 1
  return

  280 continue
!  go to i99999
  if ( i99999 == 10 ) go to 10
  if ( i99999 == 20 ) go to 20
  if ( i99999 == 200 ) go to 200

end

function error_f ( x )

!*****************************************************************************80
!
!! ERROR_F evaluates the error function.
!
!  Discussion:
!
!    Since some compilers already supply a routine named ERF which evaluates
!    the error function, this routine has been given a distinct, if
!    somewhat unnatural, name.
!    
!    The function is defined by:
!
!      ERF(X) = ( 2 / sqrt ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( - T**2 ) dT.
!
!    Properties of the function include:
!
!      Limit ( X -> -Infinity ) ERF(X) =          -1.0;
!                               ERF(0) =           0.0;
!                               ERF(0.476936...) = 0.5;
!      Limit ( X -> +Infinity ) ERF(X) =          +1.0.
!
!      0.5D+00 * ( ERF(X/sqrt(2)) + 1 ) = Normal_01_CDF(X)
!
!  Modified:
!
!    17 November 2006
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) X, the argument.
!
!    Output, real (real64) ERF, the value of the error function at X.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64), parameter, dimension ( 5 ) :: a = (/ &
    0.771058495001320D-04, &
    -0.133733772997339D-02, &
    0.323076579225834D-01, &
    0.479137145607681D-01, &
    0.128379167095513D+00 /)
  real (real64) ax
  real (real64), parameter, dimension ( 3 ) :: b = (/ &
    0.301048631703895D-02, &
    0.538971687740286D-01, &
    0.375795757275549D+00 /)
  real (real64) bot
  real (real64), parameter :: c = 0.564189583547756D+00
  real (real64) error_f
  real (real64), dimension ( 8 ) :: p = (/   &
   -1.36864857382717D-07, 5.64195517478974D-01, &
    7.21175825088309D+00, 4.31622272220567D+01, &
    1.52989285046940D+02, 3.39320816734344D+02, &
    4.51918953711873D+02, 3.00459261020162D+02 /)
  real (real64), dimension ( 8 ) :: q = (/ &
    1.00000000000000D+00, 1.27827273196294D+01, &
    7.70001529352295D+01, 2.77585444743988D+02, &
    6.38980264465631D+02, 9.31354094850610D+02, & 
    7.90950925327898D+02, 3.00459260956983D+02 /)
  real (real64), dimension ( 5 ) :: r = (/ &
    2.10144126479064D+00, 2.62370141675169D+01, &
    2.13688200555087D+01, 4.65807828718470D+00, &
    2.82094791773523D-01 /)
  real (real64), parameter, dimension ( 4 ) :: s = (/ &
    9.41537750555460D+01, 1.87114811799590D+02, &
    9.90191814623914D+01, 1.80124575948747D+02 /)
  real (real64) t
  real (real64) top
  real (real64) x
  real (real64) x2

  ax = abs ( x )

  if ( ax <= 0.5D+00 ) then

    t = x * x

    top = (((( a(1)   * t &
             + a(2) ) * t &
             + a(3) ) * t &
             + a(4) ) * t &
             + a(5) ) + 1.0D+00

    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00
    error_f = ax * ( top / bot )

  else if ( ax <= 4.0D+00 ) then

    top = (((((( p(1)   * ax &
               + p(2) ) * ax &
               + p(3) ) * ax &
               + p(4) ) * ax &
               + p(5) ) * ax &
               + p(6) ) * ax &
               + p(7) ) * ax &
               + p(8)

    bot = (((((( q(1) * ax + q(2) ) * ax + q(3) ) * ax + q(4) ) * ax &
      + q(5) ) * ax + q(6) ) * ax + q(7) ) * ax + q(8)

    error_f = 0.5D+00 &
      + ( 0.5D+00 - exp ( - x * x ) * top / bot )

  else if ( ax < 5.8D+00 ) then

    x2 = x * x
    t = 1.0D+00 / x2

    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
      + 1.0D+00

    error_f = ( c - top / ( x2 * bot )) / ax
    error_f = 0.5D+00 &
      + ( 0.5D+00 - exp ( - x2 ) * error_f )

  else

    error_f = 1.0D+00

  end if

  if ( x < 0.0D+00 ) then
    error_f = -error_f
  end if

  return
end
function error_fc ( ind, x )

!*****************************************************************************80
!
!! ERROR_FC evaluates the complementary error function.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, integer IND, chooses the scaling.
!    If IND is nonzero, then the value returned has been multiplied by
!    EXP(X*X).
!
!    Input, real (real64) X, the argument of the function.
!
!    Output, real (real64) ERROR_FC, the value of the complementary 
!    error function.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64), dimension ( 5 ) :: a = (/ &
     0.771058495001320D-04,  -0.133733772997339D-02, &
     0.323076579225834D-01,   0.479137145607681D-01, &
     0.128379167095513D+00 /)
  real (real64) ax
  real (real64), dimension(3) :: b = (/ &
    0.301048631703895D-02, &
    0.538971687740286D-01, &
    0.375795757275549D+00 /)
  real (real64) bot
  real (real64), parameter :: c = 0.564189583547756D+00
  real (real64) e
  real (real64) error_fc
  real (real64) exparg
  integer ind
  real (real64), dimension ( 8 ) :: p = (/ &
    -1.36864857382717D-07, 5.64195517478974D-01, &
     7.21175825088309D+00, 4.31622272220567D+01, &
     1.52989285046940D+02, 3.39320816734344D+02, &
     4.51918953711873D+02, 3.00459261020162D+02 /)
  real (real64), dimension ( 8 ) :: q = (/  &
    1.00000000000000D+00, 1.27827273196294D+01, &
    7.70001529352295D+01, 2.77585444743988D+02, &
    6.38980264465631D+02, 9.31354094850610D+02, &
    7.90950925327898D+02, 3.00459260956983D+02 /)
  real (real64), dimension ( 5 ) :: r = (/ &
    2.10144126479064D+00, 2.62370141675169D+01, &
    2.13688200555087D+01, 4.65807828718470D+00, &
    2.82094791773523D-01 /)
  real (real64), dimension ( 4 ) :: s = (/ &
    9.41537750555460D+01, 1.87114811799590D+02, &
    9.90191814623914D+01, 1.80124575948747D+02 /)
  real (real64) t
  real (real64) top
  real (real64) w
  real (real64) x
!
!  ABS ( X ) <= 0.5
!
  ax = abs ( x )

  if ( ax <= 0.5D+00 ) then

    t = x * x

    top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t + a(5) ) &
      + 1.0D+00

    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0D+00

    error_fc = 0.5D+00 + ( 0.5D+00 &
      - x * ( top / bot ) )

    if ( ind /= 0 ) then
      error_fc = exp ( t ) * error_fc
    end if

    return

  end if
!
!  0.5 < abs ( X ) <= 4
!
  if ( ax <= 4.0D+00 ) then

    top = (((((( p(1) * ax + p(2)) * ax + p(3)) * ax + p(4)) * ax &
      + p(5)) * ax + p(6)) * ax + p(7)) * ax + p(8)

    bot = (((((( q(1) * ax + q(2)) * ax + q(3)) * ax + q(4)) * ax &
      + q(5)) * ax + q(6)) * ax + q(7)) * ax + q(8)

    error_fc = top / bot
!
!  4 < ABS ( X )
!
  else

    if ( x <= -5.6D+00 ) then

      if ( ind == 0 ) then
        error_fc =  2.0D+00 
      else
        error_fc =  2.0D+00  * exp ( x * x )
      end if

      return

    end if

    if ( ind == 0 ) then

      if ( 100.0D+00 < x ) then
        error_fc = 0.0D+00
        return
      end if

      if ( -exparg ( 1 ) < x * x ) then
        error_fc = 0.0D+00
        return
      end if

    end if

    t = ( 1.0D+00 / x )**2

    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t &
      + 1.0D+00

    error_fc = ( c - t * top / bot ) / ax

  end if
!
!  Final assembly.
!
  if ( ind /= 0 ) then

    if ( x < 0.0D+00 ) then
      error_fc =  2.0D+00  * exp ( x * x ) - error_fc
    end if

  else

    w = x * x
    t = w
    e = w - t
    error_fc = (( 0.5D+00 &
      + ( 0.5D+00 - e ) ) * exp ( - t ) ) * error_fc

    if ( x < 0.0D+00 ) then
      error_fc =  2.0D+00  - error_fc
    end if

  end if

  return
end
function esum ( mu, x )

!*****************************************************************************80
!
!! ESUM evaluates exp ( MU + X ).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, integer MU, part of the argument.
!
!    Input, real (real64) X, part of the argument.
!
!    Output, real (real64) ESUM, the value of exp ( MU + X ).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) esum
  integer mu
  real (real64) w
  real (real64) x

  if ( x <= 0.0D+00 ) then
    if ( 0 <= mu ) then
      w = mu + x
      if ( w <= 0.0D+00 ) then
        esum = exp ( w )
        return
      end if
    end if
  else if ( 0.0D+00 < x ) then
    if ( mu <= 0 ) then
      w = mu + x
      if ( 0.0D+00 <= w ) then
        esum = exp ( w )
        return
      end if
    end if
  end if

  w = mu
  esum = exp ( w ) * exp ( x )

  return
end
function exparg ( l )

!*****************************************************************************80
!
!! EXPARG returns the largest or smallest legal argument for EXP.
!
!  Discussion:
!
!    Only an approximate limit for the argument of EXP is desired.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, integer L, indicates which limit is desired.
!    If L = 0, then the largest positive argument for EXP is desired.
!    Otherwise, the largest negative argument for EXP for which the
!    result is nonzero is desired.
!
!    Output, real (real64) EXPARG, the desired value.
!
  use iso_fortran_env, only: real64
  implicit none

  integer b
  real (real64) exparg
  integer ipmpar
  integer l
  real (real64) lnb
  integer m
!
!  Get the arithmetic base.
!
  b = ipmpar(4)
!
!  Compute the logarithm of the arithmetic base.
!
  if ( b == 2 ) then
    lnb = 0.69314718055995D+00
  else if ( b == 8 ) then
    lnb = 2.0794415416798D+00
  else if ( b == 16 ) then
    lnb = 2.7725887222398D+00
  else
    lnb = log ( real ( b, kind = real64 ) )
  end if

  if ( l /= 0 ) then
    m = ipmpar(9) - 1
    exparg = 0.99999D+00 * ( m * lnb )
  else
    m = ipmpar(10)
    exparg = 0.99999D+00 * ( m * lnb )
  end if

  return
end
function fpser ( a, b, x, eps )

!*****************************************************************************80
!
!! FPSER evaluates IX(A,B)(X) for very small B.
!
!  Discussion:
!
!    This routine is appropriate for use when 
!
!      B < min ( EPS, EPS * A ) 
!
!    and 
!
!      X <= 0.5.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, parameters of the function.
!
!    Input, real (real64) X, the point at which the function is to
!    be evaluated.
!
!    Input, real (real64) EPS, a tolerance.
!
!    Output, real (real64) FPSER, the value of IX(A,B)(X).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) an
  real (real64) b
  real (real64) c
  real (real64) eps
  real (real64) exparg
  real (real64) fpser
  real (real64) s
  real (real64) t
  real (real64) tol
  real (real64) x

  fpser = 1.0D+00

  if ( 1.0D-03 * eps < a ) then
    fpser = 0.0D+00
    t = a * log ( x )
    if ( t < exparg ( 1 ) ) then
      return
    end if
    fpser = exp ( t )
  end if
!
!  1/B(A,B) = B
!
  fpser = ( b / a ) * fpser
  tol = eps / a
  an = a + 1.0D+00
  t = x
  s = t / an

  do

    an = an + 1.0D+00
    t = x * t
    c = t / an
    s = s + c

    if ( abs ( c ) <= tol ) then
      exit
    end if

  end do

  fpser = fpser * ( 1.0D+00 + a * s )

  return
end
function gam1 ( a )

!*****************************************************************************80
!
!! GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5 <= A <= 1.5
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, forms the argument of the Gamma function.
!
!    Output, real (real64) GAM1, the value of 1 / GAMMA ( A + 1 ) - 1.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) bot
  real (real64) d
  real (real64) gam1
  real (real64), parameter, dimension ( 7 ) :: p = (/ &
     0.577215664901533D+00, -0.409078193005776D+00, &
    -0.230975380857675D+00,  0.597275330452234D-01, &
     0.766968181649490D-02, -0.514889771323592D-02, &
     0.589597428611429D-03 /)
  real (real64), dimension ( 5 ) :: q = (/ &
    0.100000000000000D+01, 0.427569613095214D+00, &
    0.158451672430138D+00, 0.261132021441447D-01, &
    0.423244297896961D-02 /)
  real (real64), dimension ( 9 ) :: r = (/ &
    -0.422784335098468D+00, -0.771330383816272D+00, &
    -0.244757765222226D+00,  0.118378989872749D+00, &
     0.930357293360349D-03, -0.118290993445146D-01, &
     0.223047661158249D-02,  0.266505979058923D-03, &
    -0.132674909766242D-03 /)
  real (real64), parameter :: s1 = 0.273076135303957D+00
  real (real64), parameter :: s2 = 0.559398236957378D-01
  real (real64) t
  real (real64) top
  real (real64) w

  d = a - 0.5D+00

  if ( 0.0D+00 < d ) then
    t = d - 0.5D+00
  else
    t = a
  end if

  if ( t == 0.0D+00 ) then

    gam1 = 0.0D+00

  else if ( 0.0D+00 < t ) then

    top = (((((    &
            p(7)   &
      * t + p(6) ) &
      * t + p(5) ) &
      * t + p(4) ) &
      * t + p(3) ) &
      * t + p(2) ) &
      * t + p(1)

    bot = ((( q(5) * t + q(4) ) * t + q(3) ) * t + q(2) ) * t &
      + 1.0D+00

    w = top / bot

    if ( d <= 0.0D+00 ) then
      gam1 = a * w
    else
      gam1 = ( t / a ) * ( ( w - 0.5D+00 ) &
        - 0.5D+00 )
    end if

  else if ( t < 0.0D+00 ) then

    top = (((((((  &
            r(9)   &
      * t + r(8) ) & 
      * t + r(7) ) &
      * t + r(6) ) &
      * t + r(5) ) &
      * t + r(4) ) &
      * t + r(3) ) &
      * t + r(2) ) &
      * t + r(1)

    bot = ( s2 * t + s1 ) * t + 1.0D+00
    w = top / bot

    if ( d <= 0.0D+00 ) then
      gam1 = a * ( ( w + 0.5D+00 ) + 0.5D+00 )
    else
      gam1 = t * w / a
    end if

  end if

  return
end
function gamma_ln1 ( a )

!*****************************************************************************80
!
!! GAMMA_LN1 evaluates ln ( Gamma ( 1 + A ) ), for -0.2 <= A <= 1.25.
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, defines the argument of the function.
!
!    Output, real (real64) GAMMA_LN1, the value of ln ( Gamma ( 1 + A ) ).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) bot
  real (real64) gamma_ln1
  real (real64), parameter :: p0 =  0.577215664901533D+00
  real (real64), parameter :: p1 =  0.844203922187225D+00
  real (real64), parameter :: p2 = -0.168860593646662D+00
  real (real64), parameter :: p3 = -0.780427615533591D+00
  real (real64), parameter :: p4 = -0.402055799310489D+00
  real (real64), parameter :: p5 = -0.673562214325671D-01
  real (real64), parameter :: p6 = -0.271935708322958D-02
  real (real64), parameter :: q1 =  0.288743195473681D+01
  real (real64), parameter :: q2 =  0.312755088914843D+01
  real (real64), parameter :: q3 =  0.156875193295039D+01
  real (real64), parameter :: q4 =  0.361951990101499D+00
  real (real64), parameter :: q5 =  0.325038868253937D-01
  real (real64), parameter :: q6 =  0.667465618796164D-03
  real (real64), parameter :: r0 = 0.422784335098467D+00
  real (real64), parameter :: r1 = 0.848044614534529D+00
  real (real64), parameter :: r2 = 0.565221050691933D+00
  real (real64), parameter :: r3 = 0.156513060486551D+00
  real (real64), parameter :: r4 = 0.170502484022650D-01
  real (real64), parameter :: r5 = 0.497958207639485D-03
  real (real64), parameter :: s1 = 0.124313399877507D+01
  real (real64), parameter :: s2 = 0.548042109832463D+00
  real (real64), parameter :: s3 = 0.101552187439830D+00
  real (real64), parameter :: s4 = 0.713309612391000D-02
  real (real64), parameter :: s5 = 0.116165475989616D-03
  real (real64) top
  real (real64) x

  if ( a < 0.6D+00 ) then

    top = (((((  &
            p6   &
      * a + p5 ) &
      * a + p4 ) &
      * a + p3 ) &
      * a + p2 ) &
      * a + p1 ) &
      * a + p0  

    bot = (((((  &
            q6   &
      * a + q5 ) &
      * a + q4 ) &
      * a + q3 ) &
      * a + q2 ) &
      * a + q1 ) &
      * a + 1.0D+00

    gamma_ln1 = -a * ( top / bot )

  else

    x = ( a - 0.5D+00 ) - 0.5D+00

    top = ((((( r5 * x + r4 ) * x + r3 ) * x + r2 ) * x + r1 ) * x + r0 ) 

    bot = ((((( s5 * x + s4 ) * x + s3 ) * x + s2 ) * x + s1 ) * x + 1.0D+00 )

    gamma_ln1 = x * ( top / bot )

  end if

  return
end
function gamma_log ( a )

!*****************************************************************************80
!
!! GAMMA_LOG evaluates ln ( Gamma ( A ) ) for positive A.
!
!  Author:
!
!    Alfred Morris,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, the argument of the function.
!    A should be positive.
!
!    Output, real (real64), GAMMA_LOG, the value of ln ( Gamma ( A ) ).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64), parameter :: c0 =  0.833333333333333D-01
  real (real64), parameter :: c1 = -0.277777777760991D-02
  real (real64), parameter :: c2 =  0.793650666825390D-03
  real (real64), parameter :: c3 = -0.595202931351870D-03
  real (real64), parameter :: c4 =  0.837308034031215D-03
  real (real64), parameter :: c5 = -0.165322962780713D-02
  real (real64), parameter :: d  =  0.418938533204673D+00
  real (real64) gamma_log
  real (real64) gamma_ln1
  integer i
  integer n
  real (real64) t
  real (real64) w

  if ( a <= 0.8D+00 ) then

    gamma_log = gamma_ln1 ( a ) - log ( a )

  else if ( a <= 2.25D+00 ) then

    t = ( a - 0.5D+00 ) - 0.5D+00
    gamma_log = gamma_ln1 ( t )

  else if ( a < 10.0D+00 ) then

    n = a - 1.25D+00
    t = a
    w = 1.0D+00
    do i = 1, n
      t = t - 1.0D+00
      w = t * w
    end do

    gamma_log = gamma_ln1 ( t - 1.0D+00 ) + log ( w )

  else

    t = ( 1.0D+00 / a )**2

    w = ((((( c5 * t + c4 ) * t + c3 ) * t + c2 ) * t + c1 ) * t + c0 ) / a

    gamma_log = ( d + w ) + ( a - 0.5D+00 ) &
      * ( log ( a ) - 1.0D+00 )

  end if

  return
end
subroutine gamma_rat1 ( a, x, r, p, q, eps )

!*****************************************************************************80
!
!! GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, X, the parameters of the functions.
!    It is assumed that A <= 1.
!
!    Input, real (real64) R, the value exp(-X) * X**A / Gamma(A).
!
!    Output, real (real64) P, Q, the values of P(A,X) and Q(A,X).
!
!    Input, real (real64) EPS, the tolerance.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) a2n
  real (real64) a2nm1
  real (real64) am0
  real (real64) an
  real (real64) an0
  real (real64) b2n
  real (real64) b2nm1
  real (real64) c
  real (real64) cma
  real (real64) eps
  real (real64) error_f
  real (real64) error_fc
  real (real64) g
  real (real64) gam1
  real (real64) h
  real (real64) j
  real (real64) l
  real (real64) p
  real (real64) q
  real (real64) r
  real (real64) rexp
  real (real64) sum1
  real (real64) t
  real (real64) tol
  real (real64) w
  real (real64) x
  real (real64) z

  if ( a * x == 0.0D+00 ) then

    if ( x <= a ) then
      p = 0.0D+00
      q = 1.0D+00
    else
      p = 1.0D+00
      q = 0.0D+00
    end if

    return
  end if

  if ( a == 0.5D+00 ) then

    if ( x < 0.25D+00 ) then
      p = error_f ( sqrt ( x ) )
      q = 0.5D+00 + ( 0.5D+00 - p )
    else
      q = error_fc ( 0, sqrt ( x ) )
      p = 0.5D+00 + ( 0.5D+00 - q )
    end if

    return

  end if
!
!  Taylor series for P(A,X)/X**A
!
  if ( x < 1.1D+00 ) then

    an = 3.0
    c = x
    sum1 = x / ( a + 3.0D+00 )
    tol = 0.1D+00 * eps / ( a + 1.0D+00 )

    do

      an = an + 1.0D+00
      c = -c * ( x / an )
      t = c / ( a + an )
      sum1 = sum1 + t

      if ( abs ( t ) <= tol ) then
        exit
      end if

    end do

    j = a * x * ( ( sum1 / 6.0D+00 - 0.5D+00 &
      / ( a +  2.0D+00  ) ) &
      * x + 1.0D+00 / ( a + 1.0D+00 ) )

    z = a * log ( x )
    h = gam1 ( a )
    g = 1.0D+00 + h

    if ( x < 0.25D+00 ) then
      go to 30
    end if

    if ( a < x / 2.59D+00 ) then
      go to 50
    else
      go to 40
    end if

30 continue

    if ( -0.13394D+00 < z ) then
      go to 50
    end if

40 continue

    w = exp ( z )
    p = w * g * ( 0.5D+00 + ( 0.5D+00 - j ))
    q = 0.5D+00 + ( 0.5D+00 - p )
    return

50 continue

    l = rexp ( z )
    w = 0.5D+00 + ( 0.5D+00 + l )
    q = ( w * j - l ) * g - h

    if  ( q < 0.0D+00 ) then
      p = 1.0D+00
      q = 0.0D+00
    else
      p = 0.5D+00 + ( 0.5D+00 - q )
    end if
!
!  Continued fraction expansion.
!
  else

    a2nm1 = 1.0D+00
    a2n = 1.0D+00
    b2nm1 = x
    b2n = x + ( 1.0D+00 - a )
    c = 1.0D+00

    do

      a2nm1 = x * a2n + c * a2nm1
      b2nm1 = x * b2n + c * b2nm1
      am0 = a2nm1 / b2nm1
      c = c + 1.0D+00
      cma = c - a
      a2n = a2nm1 + cma * a2n
      b2n = b2nm1 + cma * b2n
      an0 = a2n / b2n

      if ( abs ( an0 - am0 ) < eps * an0 ) then
        exit
      end if

    end do

    q = r * an0
    p = 0.5D+00 + ( 0.5D+00 - q )

  end if

  return
end
function gsumln ( a, b )

!*****************************************************************************80
!
!! GSUMLN evaluates the function ln(Gamma(A + B)).
!
!  Discussion:
!
!    GSUMLN is used for 1 <= A <= 2 and 1 <= B <= 2
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) A, B, values whose sum is the argument of
!    the Gamma function.
!
!    Output, real (real64) GSUMLN, the value of ln(Gamma(A+B)).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) a
  real (real64) alnrel
  real (real64) b
  real (real64) gamma_ln1
  real (real64) gsumln
  real (real64) x

  x = a + b - 2.0D+00

  if ( x <= 0.25D+00 ) then
    gsumln = gamma_ln1 ( 1.0D+00 + x )
  else if ( x <= 1.25D+00 ) then
    gsumln = gamma_ln1 ( x ) + alnrel ( x )
  else
    gsumln = gamma_ln1 ( x - 1.0D+00 ) + log ( x * ( 1.0D+00 + x ) )
  end if

  return
end
function ipmpar ( i )

!*****************************************************************************80
!
!! IPMPAR returns integer machine constants. 
!
!  Discussion:
!
!    Input arguments 1 through 3 are queries about integer arithmetic.
!    We assume integers are represented in the N-digit, base A form
!
!      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
!
!    where 0 <= X(0:N-1) < A.
!
!    Then:
!
!      IPMPAR(1) = A, the base of integer arithmetic;
!      IPMPAR(2) = N, the number of base A digits;
!      IPMPAR(3) = A**N - 1, the largest magnitude.
!
!    It is assumed that the single and real (real64) floating
!    point arithmetics have the same base, say B, and that the
!    nonzero numbers are represented in the form
!
!      sign * (B**E) * (X(1)/B + ... + X(M)/B**M)
!
!    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
!    EMIN <= E <= EMAX.
!
!    Input argument 4 is a query about the base of real arithmetic:
!
!      IPMPAR(4) = B, the base of single and real (real64) arithmetic.
!
!    Input arguments 5 through 7 are queries about single precision
!    floating point arithmetic:
!
!     IPMPAR(5) = M, the number of base B digits for single precision.
!     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
!     IPMPAR(7) = EMAX, the largest exponent E for single precision.
!
!    Input arguments 8 through 10 are queries about real (real64)
!    floating point arithmetic:
!
!     IPMPAR(8) = M, the number of base B digits for real (real64).
!     IPMPAR(9) = EMIN, the smallest exponent E for real (real64).
!     IPMPAR(10) = EMAX, the largest exponent E for real (real64).
!
!  Reference:
!
!    Phyllis Fox, Andrew Hall, Norman Schryer,
!    Algorithm 528:
!    Framework for a Portable FORTRAN Subroutine Library,
!    ACM Transactions on Mathematical Software,
!    Volume 4, 1978, pages 176-188.
!
!  Parameters:
!
!    Input, integer I, the index of the desired constant.
!
!    Output, integer IPMPAR, the value of the desired constant.
!
  use iso_fortran_env, only: real64
  implicit none

  integer i
  integer imach(10)
  integer ipmpar
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!     data imach( 1) /   2 /
!     data imach( 2) /  31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /  16 /
!     data imach( 5) /   6 /
!     data imach( 6) / -64 /
!     data imach( 7) /  63 /
!     data imach( 8) /  14 /
!     data imach( 9) / -64 /
!     data imach(10) /  63 /
!
!     Machine constants for the AT&T 3B SERIES, AT&T
!     PC 7300, AND AT&T 6300.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the BURROUGHS 1700 SYSTEM.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   33 /
!     data imach( 3) / 8589934591 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -256 /
!     data imach( 7) /  255 /
!     data imach( 8) /   60 /
!     data imach( 9) / -256 /
!     data imach(10) /  255 /
!
!     Machine constants for the BURROUGHS 5700 SYSTEM.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /    8 /
!     data imach( 5) /   13 /
!     data imach( 6) /  -50 /
!     data imach( 7) /   76 /
!     data imach( 8) /   26 /
!     data imach( 9) /  -50 /
!     data imach(10) /   76 /
!
!     Machine constants for the BURROUGHS 6700/7700 SYSTEMS.
!
!     data imach( 1) /      2 /
!     data imach( 2) /     39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /      8 /
!     data imach( 5) /     13 /
!     data imach( 6) /    -50 /
!     data imach( 7) /     76 /
!     data imach( 8) /     26 /
!     data imach( 9) / -32754 /
!     data imach(10) /  32780 /
!
!     Machine constants for the CDC 6000/7000 SERIES
!     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS OPERATING SYSTEM).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   48 /
!     data imach( 3) / 281474976710655 /
!     data imach( 4) /    2 /
!     data imach( 5) /   48 /
!     data imach( 6) / -974 /
!     data imach( 7) / 1070 /
!     data imach( 8) /   95 /
!     data imach( 9) / -926 /
!     data imach(10) / 1070 /
!
!     Machine constants for the CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS/VE OPERATING SYSTEM).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    48 /
!     data imach( 6) / -4096 /
!     data imach( 7) /  4095 /
!     data imach( 8) /    96 /
!     data imach( 9) / -4096 /
!     data imach(10) /  4095 /
!
!     Machine constants for the CRAY 1, XMP, 2, AND 3.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    47 /
!     data imach( 6) / -8189 /
!     data imach( 7) /  8190 /
!     data imach( 8) /    94 /
!     data imach( 9) / -8099 /
!     data imach(10) /  8190 /
!
!     Machine constants for the data GENERAL ECLIPSE S/200.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     Machine constants for the HARRIS 220.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   23 /
!     data imach( 3) / 8388607 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   38 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the HONEYWELL 600/6000
!     AND DPS 8/70 SERIES.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   63 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 2100
!     3 WORD real (real64) OPTION WITH FTN4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   39 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 2100
!     4 WORD real (real64) OPTION WITH FTN4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   55 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 9000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -126 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the IBM 360/370 SERIES,
!     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
!     5/7/9 AND THE SEL SYSTEMS 85/86.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     Machine constants for the IBM PC.
!
!      data imach(1)/2/
!      data imach(2)/31/
!      data imach(3)/2147483647/
!      data imach(4)/2/
!      data imach(5)/24/
!      data imach(6)/-125/
!      data imach(7)/128/
!      data imach(8)/53/
!      data imach(9)/-1021/
!      data imach(10)/1024/
!
!     Machine constants for the MACINTOSH II - ABSOFT
!     MACFORTRAN II.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the MICROVAX - VMS FORTRAN.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the PDP-11 FORTRAN SUPPORTING
!     32-BIT integer ARITHMETIC.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the SEQUENT BALANCE 8000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the SILICON GRAPHICS IRIS-4D
!     SERIES (MIPS R3000 PROCESSOR).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
  data imach( 1) /     2 /
  data imach( 2) /    31 /
  data imach( 3) / 2147483647 /
  data imach( 4) /     2 /
  data imach( 5) /    24 /
  data imach( 6) /  -125 /
  data imach( 7) /   128 /
  data imach( 8) /    53 /
  data imach( 9) / -1021 /
  data imach(10) /  1024 /
!
!     Machine constants for the UNIVAC 1100 SERIES.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   60 /
!     data imach( 9) /-1024 /
!     data imach(10) / 1023 /
!
!     Machine constants for the VAX 11/780.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
  ipmpar = imach(i)

  return
end
function psi ( xx )

!*****************************************************************************80
!
!! PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
!
!  Discussion:
!
!    The main computation involves evaluation of rational Chebyshev
!    approximations.  PSI was written at Argonne National Laboratory 
!    for FUNPACK, and subsequently modified by A. H. Morris of NSWC.
!
!  Reference:
!
!    William Cody, Anthony Strecok, Henry Thacher,
!    Chebyshev Approximations for the Psi Function,
!    Mathematics of Computation,
!    Volume 27, 1973, pages 123-127.
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) XX, the argument of the psi function.
!
!    Output, real (real64) PSI, the value of the psi function.  PSI 
!    is assigned the value 0 when the psi function is undefined.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64) aug
  real (real64) den
  real (real64), parameter :: dx0 = &
    1.461632144968362341262659542325721325D+00
  integer i
  integer ipmpar
  integer m
  integer n
  integer nq
  real (real64), parameter, dimension ( 7 ) :: p1 = (/ &
   0.895385022981970D-02, &
   0.477762828042627D+01, &
   0.142441585084029D+03, &
   0.118645200713425D+04, &
   0.363351846806499D+04, &
   0.413810161269013D+04, &
   0.130560269827897D+04/)
  real (real64), dimension ( 4 ) :: p2 = (/ &
    -0.212940445131011D+01, &
    -0.701677227766759D+01, &
    -0.448616543918019D+01, &
    -0.648157123766197D+00 /)
  real (real64), parameter :: piov4 = 0.785398163397448D+00
  real (real64) psi
!
!  Coefficients for rational approximation of
!  PSI(X) / (X - X0),  0.5D+00 <= X <= 3.0D+00
!
  real (real64), dimension ( 6 ) :: q1 = (/ &
    0.448452573429826D+02, &
    0.520752771467162D+03, &
    0.221000799247830D+04, &
    0.364127349079381D+04, &
    0.190831076596300D+04, &
    0.691091682714533D-05 /)
  real (real64), dimension ( 4 ) :: q2 = (/ &
    0.322703493791143D+02, &
    0.892920700481861D+02, &
    0.546117738103215D+02, &
    0.777788548522962D+01 /)
  real (real64) sgn
  real (real64) upper
  real (real64) w
  real (real64) x
  real (real64) xmax1
  real (real64) xmx0
  real (real64) xsmall
  real (real64) xx
  real (real64) z
!
!  XMAX1 is the largest positive floating point constant with entirely 
!  integer representation.  It is also used as negative of lower bound 
!  on acceptable negative arguments and as the positive argument beyond which
!  psi may be represented as LOG(X).
!
  xmax1 = real ( ipmpar(3), kind = real64 )
  xmax1 = min ( xmax1, 1.0D+00 / epsilon ( xmax1 ) )
!
!  XSMALL is the absolute argument below which PI*COTAN(PI*X)
!  may be represented by 1/X.
!
  xsmall = 1.0D-09

  x = xx
  aug = 0.0D+00

  if ( x == 0.0D+00 ) then
    psi = 0.0D+00
    return
  end if
!
!  X < 0.5,  Use reflection formula PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!
  if ( x < 0.5D+00 ) then
!
!  0 < ABS ( X ) <= XSMALL.  Use 1/X as a substitute for PI*COTAN(PI*X)
!
    if ( abs ( x ) <= xsmall ) then
      aug = -1.0D+00 / x
      go to 40
    end if
!
!  Reduction of argument for cotangent.
!
    w = -x
    sgn = piov4

    if ( w <= 0.0D+00 ) then
      w = -w
      sgn = -sgn
    end if
!
!  Make an error exit if X <= -XMAX1
!
    if ( xmax1 <= w ) then
      psi = 0.0D+00
      return
    end if

    nq = int ( w )
    w = w - real ( nq, kind = real64 )
    nq = int ( w * 4.0D+00 )
    w = 4.0D+00 * ( w - real ( nq, kind = real64 ) * 0.25D+00 )
!
!  W is now related to the fractional part of 4.0D+00 * X.
!  Adjust argument to correspond to values in first
!  quadrant and determine sign.
!
    n = nq / 2
    if ( n + n /= nq ) then
      w = 1.0D+00 - w
    end if

    z = piov4 * w
    m = n / 2

    if ( m + m /= n ) then
      sgn = -sgn
    end if
!
!  Determine final value for -PI * COTAN(PI*X).
!
    n = ( nq + 1 ) / 2
    m = n / 2
    m = m + m

    if ( m == n ) then

      if ( z == 0.0D+00 ) then
        psi = 0.0D+00
        return
      end if

      aug = 4.0D+00 * sgn * ( cos(z) / sin(z) )

    else

      aug = 4.0D+00 * sgn * ( sin(z) / cos(z) )

    end if

   40   continue

    x = 1.0D+00 - x

  end if
!
!  0.5 <= X <= 3 
!
  if ( x <= 3.0D+00 ) then

    den = x
    upper = p1(1) * x

    do i = 1, 5
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do

    den = ( upper + p1(7) ) / ( den + q1(6) )
    xmx0 = real ( x, kind = real64 ) - dx0
    psi = den * xmx0 + aug
!
!  3 < X < XMAX1
!
  else if ( x < xmax1 ) then

    w = 1.0D+00 / x**2
    den = w
    upper = p2(1) * w

    do i = 1, 3
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do

    aug = upper / ( den + q2(4) ) - 0.5D+00 / x + aug
    psi = aug + log ( x )
!
!  XMAX1 <= X
!
  else

    psi = aug + log ( x )

  end if

  return
end
function rexp ( x )

!*****************************************************************************80
!
!! REXP evaluates the function EXP(X) - 1.
!
!  Modified:
!
!    09 December 1999
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) X, the argument of the function.
!
!    Output, real (real64) REXP, the value of EXP(X)-1.
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64), parameter :: p1 =  0.914041914819518D-09
  real (real64), parameter :: p2 =  0.238082361044469D-01
  real (real64), parameter :: q1 = -0.499999999085958D+00
  real (real64), parameter :: q2 =  0.107141568980644D+00
  real (real64), parameter :: q3 = -0.119041179760821D-01
  real (real64), parameter :: q4 =  0.595130811860248D-03
  real (real64) rexp
  real (real64) w
  real (real64) x

  if ( abs ( x ) <= 0.15D+00 ) then

    rexp = x * ( ( ( p2 * x + p1 ) * x + 1.0D+00 ) &
      / ( ( ( ( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0D+00 ) )

  else

    w = exp ( x )

    if ( x <= 0.0D+00 ) then
      rexp = ( w - 0.5D+00 ) - 0.5D+00
    else
      rexp = w * ( 0.5D+00 + ( 0.5D+00 - 1.0D+00 / w ) )
    end if

  end if

  return
end
function rlog1 ( x )

!*****************************************************************************80
!
!! RLOG1 evaluates the function X - ln ( 1 + X ).
!
!  Author:
!
!    Armido DiDinato, Alfred Morris
!
!  Reference:
!
!    Armido DiDinato, Alfred Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!    Input, real (real64) X, the argument.
!
!    Output, real (real64) RLOG1, the value of X - ln ( 1 + X ).
!
  use iso_fortran_env, only: real64
  implicit none

  real (real64), parameter :: a = 0.566749439387324D-01
  real (real64), parameter :: b = 0.456512608815524D-01
  real (real64) h
  real (real64), parameter :: half = 0.5D+00
  real (real64), parameter :: p0 = 0.333333333333333D+00
  real (real64), parameter :: p1 = -0.224696413112536D+00
  real (real64), parameter :: p2 = 0.620886815375787D-02
  real (real64), parameter :: q1 = -0.127408923933623D+01
  real (real64), parameter :: q2 = 0.354508718369557D+00
  real (real64) r
  real (real64) rlog1
  real (real64) t
  real (real64), parameter :: two =  2.0D+00 
  real (real64) w
  real (real64) w1
  real (real64) x

  if ( x < -0.39D+00 ) then

    w = ( x + half ) + half
    rlog1 = x - log ( w )

  else if ( x < -0.18D+00 ) then

    h = x + 0.3D+00
    h = h / 0.7D+00
    w1 = a - h * 0.3D+00

    r = h / ( h + 2.0D+00 )
    t = r * r
    w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
    rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

  else if ( x <= 0.18D+00 ) then

    h = x
    w1 = 0.0D+00

    r = h / ( h + two )
    t = r * r
    w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
    rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

  else if ( x <= 0.57D+00 ) then

    h = 0.75D+00 * x - 0.25D+00
    w1 = b + h / 3.0D+00

    r = h / ( h + 2.0D+00 )
    t = r * r
    w = ( ( p2 * t + p1 ) * t + p0 ) / ( ( q2 * t + q1 ) * t + 1.0D+00 )
    rlog1 = two * t * ( 1.0D+00 / ( 1.0D+00 - r ) - r * w ) + w1

  else 

    w = ( x + half ) + half
    rlog1 = x - log ( w )

  end if

  return
end
