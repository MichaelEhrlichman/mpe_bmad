module analytic_rm_mod

contains
!----------------------------------------------------
!+
! Subroutine analytic_orm
!-
subroutine analytic_orm(lat,I1,plane,A,bpms,correctors)
  use bmad
  implicit none

  type(lat_struct) lat
  real(rp) I1
  character(1) plane
  real(rp) A(:,:)
  type(ele_pointer_struct) :: bpms(:)
  type(ele_pointer_struct) :: correctors(:)
  type(ele_struct) ele_at_s

  integer :: n_bpms
  integer :: n_correctors
  integer status

  integer i,j
  real(rp) beta_1, beta_2
  real(rp) eta_1, eta_2
  real(rp) phi_1, phi_2
  real(rp) tune_factor
  real(rp) alpha_c, alpha_factor, gamma
  real(rp) midpoint

  n_bpms = size(bpms)
  n_correctors = size(correctors)

  if( plane == 'x' ) then
    tune_factor = lat%ele(lat%n_ele_track)%a%phi / 2
  elseif( plane == 'y' ) then
    tune_factor = lat%ele(lat%n_ele_track)%b%phi / 2
  endif

  alpha_c = I1/lat%param%total_length
  call convert_total_energy_to(lat%ele(1)%value(e_tot$),electron$,gamma=gamma)
  alpha_factor = (alpha_c - 1.0d0/gamma/gamma)*lat%param%total_length

  do i=1,n_bpms
    if( plane == 'x' ) then
      beta_1 = bpms(i)%ele%a%beta
      eta_1  = bpms(i)%ele%a%eta
      phi_1  = bpms(i)%ele%a%phi
    elseif( plane == 'y' ) then
      beta_1 = bpms(i)%ele%b%beta
      eta_1  = bpms(i)%ele%b%eta
      phi_1  = bpms(i)%ele%b%phi
    endif
    do j=1,n_correctors
      midpoint = correctors(j)%ele%s - correctors(j)%ele%value(l$)/2.0d0
      call twiss_and_track_at_s(lat,midpoint,ele_at_s)
      if( plane == 'x' ) then
        beta_2 = ele_at_s%a%beta
        eta_2  = ele_at_s%a%eta
        phi_2  = ele_at_s%a%phi
        !beta_2 = correctors(j)%ele%a%beta
        !eta_2  = correctors(j)%ele%a%eta
        !phi_2  = correctors(j)%ele%a%phi
      elseif( plane == 'y' ) then
        beta_2 = ele_at_s%b%beta
        eta_2  = ele_at_s%b%eta
        phi_2  = ele_at_s%b%phi
        !beta_2 = correctors(j)%ele%b%beta
        !eta_2  = correctors(j)%ele%b%eta
        !phi_2  = correctors(j)%ele%b%phi
      endif
      A(i,j) = sqrt(beta_1*beta_2)/2/sin(tune_factor) * cos( abs(phi_1-phi_2) - tune_factor ) &
               - eta_1*eta_2/alpha_factor
    enddo
  enddo
end subroutine analytic_orm

!----------------------------------------------------
!+
! Subroutine analytic_prm
!-
subroutine analytic_prm(lat,A,bpms,corr_quads)
  use bmad
  implicit none

  type(lat_struct) lat
  real(rp) A(:,:)
  type(ele_pointer_struct) :: bpms(:)
  type(ele_pointer_struct) :: corr_quads(:)

  integer :: n_bpms
  integer :: n_correctors
  integer status

  integer i,j
  real(rp) beta_i, beta_j
  real(rp) phi_i, phi_j
  real(rp) nu, l_k

  n_bpms = size(bpms)
  n_correctors = size(corr_quads)

  !phase response
  nu = lat%ele(lat%n_ele_track)%a%phi
  do i=1,n_bpms
    beta_i = bpms(i)%ele%a%beta
    phi_i  = bpms(i)%ele%a%phi
    do j=1,n_correctors
      beta_j = corr_quads(j)%ele%a%beta
      phi_j  = corr_quads(j)%ele%a%phi
      A(i,j) = beta_j/4.0d0/sin(nu) * ( sin(nu)+sin(2*phi_j-nu) + sign(1.0d0,phi_i-phi_j)*(sin(nu)+sin(2*abs(phi_i-phi_j)-nu ) ) )
    enddo
  enddo
  nu = lat%ele(lat%n_ele_track)%b%phi
  do i=1,n_bpms
    beta_i = bpms(i)%ele%b%beta
    phi_i  = bpms(i)%ele%b%phi
    do j=1,n_correctors
      beta_j = corr_quads(j)%ele%b%beta
      phi_j  = corr_quads(j)%ele%b%phi
      A(i+n_bpms,j) = beta_j/4.0d0/sin(nu) * ( sin(nu)+sin(2*phi_j-nu) + sign(1.0d0,phi_i-phi_j)*(sin(nu)+sin(2*abs(phi_i-phi_j)-nu ) ) )
    enddo
  enddo
end subroutine analytic_prm

!+
! Subroutine numerical_drm
!-
subroutine numerical_drm(lat,A,bpms,corr_quads)
  use bmad
  implicit none

  type(lat_struct) lat
  real(rp) A(:,:)
  type(ele_pointer_struct) :: bpms(:)
  type(ele_pointer_struct) :: corr_quads(:)

  type(coord_struct), allocatable :: co(:)

  logical err

  integer n_bpms
  integer n_correctors
  integer i,j
  integer status

  real(rp) meas0a(size(bpms))
  real(rp) meas0b(size(bpms))
  real(rp), parameter :: delta = 0.0001

  n_bpms = size(bpms)
  n_correctors = size(corr_quads)

  allocate(co(0:lat%n_ele_track))

  !dispersion response
  do i=1,n_bpms
    meas0a(i) = bpms(i)%ele%a%eta
    meas0b(i) = bpms(i)%ele%b%eta
  enddo
  do j=1,n_correctors
    corr_quads(j)%ele%value(k1$) = corr_quads(j)%ele%value(k1$) + delta
    call set_flags_for_changed_attribute(corr_quads(j)%ele, corr_quads(j)%ele%value(k1$))
    call lattice_bookkeeper(lat)
    call twiss_and_track(lat,co,status)
    do i=1,n_bpms
      A(i,j) = -(meas0a(i)-bpms(i)%ele%a%eta)/delta
      !A(i+n_bpms,j) = -(meas0b(i)-bpms(i)%ele%b%eta)/delta
    enddo
    corr_quads(j)%ele%value(k1$) = corr_quads(j)%ele%value(k1$) - delta
    call set_flags_for_changed_attribute(corr_quads(j)%ele, corr_quads(j)%ele%value(k1$))
  enddo
  call lattice_bookkeeper(lat)
  call twiss_and_track(lat,co,status)
end subroutine

end module






