module linear_aperture_mod

implicit none

contains

subroutine linear_aperture(ring,da_config)
  !works OK with dispersion at injection point.  TESTED!
  !assumes all apertures are elliptical
  use bmad
  use custom_dynamic_aperture_mod, only: custom_aperture_scan_struct, local_cosphi, local_sinphi
  use transfer_map_mod

  implicit none

  type(lat_struct) ring
  type(custom_aperture_scan_struct) da_config

  integer i, j, p
  real(rp) tmat(6,6), t1(6,6)
  real(rp) tvec(6), v1(6)
  real(rp) delta_angle
  real(rp) x1,y1,x2,y2
  real(rp) dE, theta
  real(rp) a,b
  real(rp) x1ds,y1ds,x2ds,y2ds
  real(rp) m,k,pl2,det,Ex1,Ey1,Ex2,Ey2
  real(rp) da2
  real(rp) xa, ya, xb, yb, xpick, ypick
  real(rp) vec_length
  real(rp) offx, offy
  logical mask_x(1:ring%n_ele_track)
  logical mask_y(1:ring%n_ele_track)
  integer check

  delta_angle = (da_config%max_angle - da_config%min_angle)/(da_config%n_angle -1)
  if(delta_angle .ne. delta_angle) delta_angle = 0.0d0  !case of n_angle==1.  In fortran NaN does not equal itself.
  a = ring%ele(1)%value(x1_limit$)
  b = ring%ele(1)%value(y1_limit$)

  mask_x = abs(ring%ele(1:ring%n_ele_track)%value(x1_limit$)) > 0.0001
  mask_y = abs(ring%ele(1:ring%n_ele_track)%value(y1_limit$)) > 0.0001
  da_config%Sx = minval( ring%ele(1:ring%n_ele_track)%value(x1_limit$) / sqrt(ring%ele(1:ring%n_ele_track)%a%beta), mask_x ) * sqrt(ring%ele(1)%a%beta)
  da_config%Sy = minval( ring%ele(1:ring%n_ele_track)%value(y1_limit$) / sqrt(ring%ele(1:ring%n_ele_track)%b%beta), mask_y ) * sqrt(ring%ele(1)%b%beta)

  !initialize to physical aperture at track point
  do i=1,da_config%n_angle
    theta = (i-1)*delta_angle + da_config%min_angle
    vec_length = sqrt(a**2 * cos(theta)**2 + a**2 * sin(theta)**2)
    da_config%aperture(i)%x = vec_length * local_cosphi(theta,da_config%Sx,da_config%Sy)
    da_config%aperture(i)%y = vec_length * local_sinphi(theta,da_config%Sx,da_config%Sy)
    da_config%aperture(i)%i_turn = 1  !this routine does not calculate the turn number at which the particle was lost.
  enddo

  do i=1,ring%n_ele_track
    a = ring%ele(i)%value(x1_limit$)  !horizontal axis of ellipse
    b = ring%ele(i)%value(y1_limit$)  !vertical axis of ellipse

    call transfer_matrix_calc (ring, t1, v1, ix1=i, one_turn=.true.)
    call transfer_matrix_calc (ring, tmat, tvec, ix1=0, ix2=i)
    do j=1,100  !simulate 100 turns
      if (j>1) then
        call concat_transfer_mat (t1, v1, tmat, tvec, tmat, tvec)
      endif
      do p=1,da_config%n_angle
        theta = (p-1)*delta_angle + da_config%min_angle  !unscaled angle at injection point

        x1 = local_cosphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(1)  ! cos (scaled theta at injection point)
        y1 = local_sinphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(3)  ! sin (scaled theta at injection point)
        x2 = 0.5*local_cosphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(1)  ! cos (scaled theta at injection point)
        y2 = 0.5*local_sinphi(theta,da_config%Sx,da_config%Sy) + da_config%param%closed_orbit%vec(3)  ! sin (scaled theta at injection point)

        offx = tmat(1,6)*da_config%param%closed_orbit%vec(6) + tvec(1)
        offy = tmat(3,6)*da_config%param%closed_orbit%vec(6) + tvec(3)

        x1ds = tmat(1,1)*x1 + tmat(1,3)*y1 + offx  ! horiz component of (xth,yth) transported downstream
        y1ds = tmat(3,1)*x1 + tmat(3,3)*y1 + offy  ! vertical component of (xth,yth) transported downstream
        x2ds = tmat(1,1)*x2 + tmat(1,3)*y2 + offx  ! horiz component of (xth,yth) transported downstream
        y2ds = tmat(3,1)*x2 + tmat(3,3)*y2 + offy  ! vertical component of (xth,yth) transported downstream

        ! y=mx+k form of downstream vector
        m = (y2ds-y1ds)/(x2ds-x1ds)
        k = y1ds - m*x1ds

        !Calculate where downstream vector intersects elliptical beam pipe.
        !These equations are the cos and sin of the intersection of an ellipse
        !and a line described by y=mx+k.  The axes of the ellipse are a and b.
        Ex1 = a * (-a*k*m-b*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)
        Ex2 = a * (-a*k*m+b*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)
        Ey1 = b * (b*k-a*m*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)
        Ey2 = b * (b*k+a*m*sqrt(b*b-k*k+a*a*m*m))/(b*b+a*a*m*m)

        !Transport vector back to the injection point
        det = tmat(1,1)*tmat(3,3) - tmat(1,3)*tmat(3,1)
        xa = (tmat(3,3)*(Ex1-offx) - tmat(1,3)*(Ey1-offy))/det  !aperture projected back to injection point
        ya = (-tmat(3,1)*(Ex1-offx) + tmat(1,1)*(Ey1-offy))/det  !aperture projected back to injection point
        xb = (tmat(3,3)*(Ex2-offx) - tmat(1,3)*(Ey2-offy))/det  !aperture projected back to injection point
        yb = (-tmat(3,1)*(Ex2-offx) + tmat(1,1)*(Ey2-offy))/det  !aperture projected back to injection point

        if( ya > 0 .and. yb < 0) then
          xpick = xa
          ypick = ya
          check=1
        elseif( yb > 0 .and. ya < 0) then
          xpick = xb
          ypick = yb
          check=2
        else
          !ERROR.  Happens at last element in lattice.  Not sure why.
          xpick = 1.0d0
          ypick = 1.0d0
          check=3
        endif
        !write(*,'(a,2es14.5,i6)') "CHECK: ", xpick, ypick, check

        pl2 = (xpick-da_config%param%closed_orbit%vec(1))**2 + (ypick-da_config%param%closed_orbit%vec(3))**2
        da2 = (da_config%aperture(p)%x-da_config%param%closed_orbit%vec(1))**2 + (da_config%aperture(p)%y-da_config%param%closed_orbit%vec(3))**2

        if( pl2 .lt. da2 ) then
          da_config%aperture(p)%x = xpick
          da_config%aperture(p)%y = ypick
        endif
      enddo
    enddo
  enddo
end subroutine

end module





