module compute
  use fbauhube
  implicit none

  real(8), parameter :: pi = acos(-1d0)
  private :: pi

contains
  subroutine frill(N, aC, rC, x_full, y_full, n_xy)
      real(8) :: aC, rC
      integer :: N
      integer, intent(out) :: n_xy
      real(8), dimension(4*N), intent(out) :: x_full, y_full

      real(8) :: rho, sa, sb, sc, sx, sy
      real(8) :: dif, pha, sa2, phi
      real(8) :: r, q, p
      real(8), dimension(0:NMAX) :: ai,ar
      real(8), dimension(0:NMAX) :: rootsi, rootsr, val 
      real(8) :: t1, t2
      integer :: k,j,lol

      ! write(*,*) rC, aC

      sx = rC*cos(aC)
      sy = rC*sin(aC)

      x_full = 0d0; y_full = 0d0
      n_xy = 0

      do j=1,N
        sy = rC*(j*2d0/N-1d0+1.1d0/N)
        sx = sqrt(rC*rC-sy*sy)
        phi = asin(j*2d0/N-1d0+1d0/N)
        r = 1d0 - rC*rC
        q = 2d0*sx*sx- 2d0*sy*sy
        p = r-3d0
        ! ar = [1d0, 0d0, p, q, r]
        ar(:4) = [r, q, p, 0d0, 1d0]
        ai = 0d0
        lol = bauhub (0, 0, 4, ar, ai, rootsr, rootsi, val)
        ! write(*,*) sum(val(:3))/4
        ! P = poly1d([1,0,p,q,r])
        do k=0,3
          if(abs(rootsi(k))<1d-30 .and. (1d-14<= rootsr(k) .and. rootsr(k)<=1d0)) then
            rho = sqrt(rootsr(k))
            n_xy = n_xy + 1
            x_full(n_xy) = rho*cos(aC+phi)
            y_full(n_xy) = rho*sin(aC+phi)
            n_xy = n_xy + 1
            x_full(n_xy) = -rho*cos(aC+phi)
            y_full(n_xy) = -rho*sin(aC+phi)
          endif
        enddo
      enddo

  end subroutine frill


  subroutine fundamental(S, C, M, N, a, mabs)
    complex(8) :: S, C, M
    integer :: N

    integer :: i 
    complex(8) :: an
    real(8), intent(out) :: mabs
    complex(8), dimension(N), intent(out) :: a

    a = 0d0
    an = 1d0
    mabs = 1d0
    do i=1,N
      a(i) = an
      if(abs(an)> mabs) then
        mabs = abs(an)
      endif
      an = an *(S + C*sin(i*M))
    enddo
  end subroutine fundamental


end module compute