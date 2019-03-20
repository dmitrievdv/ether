module compute
  implicit none

  real(8), parameter :: pi = acos(-1d0)
  private :: pi

contains
  subroutine frill(N, aC, rC, x_full, y_full, n_xy)
      real(8) :: aC, rC
      integer :: N
      integer, intent(out) :: n_xy
      real(8), dimension(N*N), intent(out) :: x_full, y_full

      real(8) :: rho, sa, sb, sc, sx, sy
      real(8) :: dif, pha, sa2, phi
      integer :: k,j

      ! write(*,*) rC, aC

      sx = rC*cos(aC)
      sy = rC*sin(aC)

      x_full = 0d0; y_full = 0d0
      n_xy = 0

      do k=1,N
        rho = k*1d0/N
        sa = 1d0+rho*rho
        sb = 1d0-rho*rho
        sc = sqrt(sa*sa-sb*sb)
        dif = 10d0
        pha = 0d0
        do j=1,N
          phi = j*2d0*pi/N
          sa2 = sqrt((sc*cos(phi) - sx)**2 + (sc*sin(phi) - sy)**2)  + sqrt((sc*cos(phi) + sx)**2 + (sc*sin(phi) + sy)**2)
          if(abs(sa2-2d0*sa)<1d0/4d0/N) then
            n_xy = n_xy + 1
            x_full(n_xy) = rho*cos(phi)
            y_full(n_xy) = rho*sin(phi)
          endif
        enddo
      enddo

      ! write(*,*) n_xy, x_full(:n_xy)

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