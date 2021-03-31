MODULE mc_ibs_mod

USE bmad

IMPLICIT none

TYPE grid_struct
  REAL(rp), ALLOCATABLE :: x(:)
  REAL(rp), ALLOCATABLE :: y(:)
  REAL(rp), ALLOCATABLE :: z(:)
END TYPE grid_struct

TYPE cell_to_particles_struct
  INTEGER, ALLOCATABLE :: offsets(:)
  INTEGER, ALLOCATABLE :: N(:)
  INTEGER, ALLOCATABLE :: particles(:)
END TYPE cell_to_particles_struct

TYPE xyz_struct
  REAL(rp) vec(1:3)
END TYPE xyz_struct

REAL(rp), PARAMETER :: me_kgs = 9.109d-31

CONTAINS

!+
! subroutine make_grid
!-
SUBROUTINE make_grid(dims,nx,ny,nz,grid)

REAL(rp) dims(1:6)
INTEGER nx, ny, nz
TYPE(grid_struct) :: grid

INTEGER i

IF( .not. ALLOCATED(grid%x)) ALLOCATE(grid%x(1:nx+1))
IF( .not. ALLOCATED(grid%y)) ALLOCATE(grid%y(1:ny+1))
IF( .not. ALLOCATED(grid%z)) ALLOCATE(grid%z(1:nz+1))

DO i=1,nx+1
  grid%x(i) = dims(1) + (dims(2)-dims(1))/nx*(i-1)
ENDDO
DO i=1,ny+1
  grid%y(i) = dims(3) + (dims(4)-dims(3))/ny*(i-1)
ENDDO
DO i=1,nz+1
  grid%z(i) = dims(5) + (dims(6)-dims(5))/nz*(i-1)
ENDDO

END SUBROUTINE make_grid

!+
! subroutine bin_dist
!-
SUBROUTINE bin_dist(dist,nx,ny,nz,cell_volume,particle_to_cell,cell_to_particles)

USE eigen_mod

TYPE(grid_struct) grid
INTEGER nx, ny, nz
REAL(rp) dx,dy,dz
REAL(rp) cell_volume
TYPE(coord_struct) dist(:)
TYPE(xyz_struct), ALLOCATABLE, SAVE :: edist(:)
TYPE(xyz_strucT) temp_dist
INTEGER Nmp
INTEGER particle_to_cell(:)
TYPE(cell_to_particles_struct) cell_to_particles
REAL(rp), PARAMETER :: border = 1.0d-9

INTEGER, ALLOCATABLE :: num_so_far(:)
INTEGER cell
INTEGER Ncell
INTEGER x_loc, y_loc, z_loc
INTEGER i, j

REAL(rp) xyz_sigma_mat(1:3,1:3)
REAL(rp) evals(1:3)
REAL(rp) evecs(1:3,1:3)
REAL(rp) evecsT(1:3,1:3)
REAL(rp) dims(1:6)
INTEGER  etypes(1:3)

logical error

Nmp = SIZE(dist)
IF( .not. ALLOCATED(edist) ) ALLOCATE(edist(Nmp))

CALL make_xyz_sigma_mat(dist,xyz_sigma_mat)
CALL eigensys(xyz_sigma_mat, evals, evecs, etypes, 3, error)
evecsT = TRANSPOSE(evecs)

DO i=1,Nmp
  temp_dist%vec(1) = dist(i)%vec(1)
  temp_dist%vec(2) = dist(i)%vec(3)
  temp_dist%vec(3) = dist(i)%vec(5)
  edist(i)%vec = MATMUL(evecsT,temp_dist%vec)
ENDDO

! WRITE(1000,'(A)') "# Coordinates are BMAD canonical."
! WRITE(1000,'(7A14)') "# particle", "x", "px", "y", "py", "z", "pz"
! DO i=1,Nmp
!   WRITE(1000,'(I14,6ES14.4)') i, edist(i)%vec(1:6)
! ENDDO
! CLOSE(1000)
! STOP

!Make the grid based on beam dimensions.
!Update this routine if outlayers are a problem. (bin based on sigma matrix?)
dims(1) = MINVAL(edist(:)%vec(1)) - border
dims(2) = MAXVAL(edist(:)%vec(1)) + border
dims(3) = MINVAL(edist(:)%vec(2)) - border
dims(4) = MAXVAL(edist(:)%vec(2)) + border
dims(5) = MINVAL(edist(:)%vec(3)) - border
dims(6) = MAXVAL(edist(:)%vec(3)) + border

CALL make_grid(dims,nx,ny,nz,grid)

dx = (dims(2)-dims(1)) / nx
dy = (dims(4)-dims(3)) / ny
dz = (dims(6)-dims(5)) / nz
cell_volume = dx*dy*dz

!Create map from particle number to cell number (many to 1)
DO i=1,Nmp
  DO j=2,nx+1
    IF( edist(i)%vec(1) .le. grid%x(j) ) THEN
      x_loc = j-1
      EXIT
    ENDIF
  ENDDO
  DO j=2,ny+1
    IF( edist(i)%vec(2) .le. grid%y(j) ) THEN
      y_loc = j-1
      EXIT
    ENDIF
  ENDDO
  DO j=2,nz+1
    IF( edist(i)%vec(3) .le. grid%z(j) ) THEN
      z_loc = j-1
      EXIT
    ENDIF
  ENDDO
  particle_to_cell(i) = x_loc + nx*(y_loc-1) + nx*ny*(z_loc-1)
ENDDO

!Count number of particles in each cell
cell_to_particles%N(:) = 0
DO i=1,Nmp
  cell_to_particles%N(particle_to_cell(i)) = cell_to_particles%N(particle_to_cell(i)) + 1
ENDDO

Ncell = nx*ny*nz
!Count up offsets
cell_to_particles%offsets(1) = 1
DO i=2,Ncell
  cell_to_particles%offsets(i) = cell_to_particles%offsets(i-1) + cell_to_particles%N(i-1)
ENDDO

!Create map from cell number to particle number (1 to many)
ALLOCATE(num_so_far(1:Ncell))
num_so_far(:) = 0
DO i=1,Nmp
  cell = particle_to_cell(i)
  cell_to_particles%particles( cell_to_particles%offsets(cell) + num_so_far(cell) ) = i
  num_so_far(cell) = num_so_far(cell) + 1
ENDDO
DEALLOCATE(num_so_far)

END SUBROUTINE bin_dist

!+
!-
SUBROUTINE calc_variance(density,tau_rest_frame,sigma_y,delta_u,delta_t,variance)
  REAL(rp) density
  REAL(rp) tau_rest_frame
  REAL(rp) sigma_y
  REAL(rp) delta_u
  REAL(rp) delta_t
  REAL(rp) variance

  REAL(rp) bmin
  REAL(rp) Clog

  REAL(rp), PARAMETER :: echarge4 = 6.589333674d-76
  REAL(rp), PARAMETER :: eps_0_vac2 = 7.839664195d-23
  REAL(rp), PARAMETER :: four_pi = 12.56637061d0
  REAL(rp), PARAMETER :: me_kgs2 = 8.2973881d-61
  REAL(rp), PARAMETER :: consts = echarge4 / me_kgs2 / eps_0_vac2 / four_pi

  IF( density .lt. 1 ) THEN  !check for zero density, i.e. no IBS
    variance = 0.0d0
  ELSE
    !transparent code
    !bmin = SQRT( 1 / tau_rest_frame / density / pi / delta_u  )
    !Clog = LOG(sigma_y/bmin)

    !faster code
    bmin = tau_rest_frame * density * pi * delta_u
    Clog = 0.5_rp*LOG(sigma_y*sigma_y*bmin)

    variance = consts * density * Clog / delta_u**3 * delta_t
    variance = SQRT(variance)
  ENDIF
END SUBROUTINE calc_variance

!------------------------------------------------------------------------

!+
!-
SUBROUTINE make_sigma_mat(canonical_dist,sigma_mat,sigma_beta_mat)
  USE bmad

  IMPLICIT none

  TYPE(coord_struct) canonical_dist(:)
  REAL(rp) sigma_mat(:,:)
  REAL(rp), OPTIONAL :: sigma_beta_mat(:,:)

  INTEGER j, k
  INTEGER Nparts

  Nparts = SIZE(canonical_dist)

  DO j=1,6
    DO k=j,6
      !Fortran * on two arrays is element-wise multiplication
      sigma_mat(j,k) = SUM(canonical_dist(:)%vec(j)*canonical_dist(:)%vec(k)) / Nparts
      IF( j .ne. k ) THEN
        sigma_mat(k,j) = sigma_mat(j,k)
      ENDIF
    ENDDO
  ENDDO

  IF( PRESENT(sigma_beta_mat) ) THEN
    DO j=1,6
      DO k=j,6
        sigma_beta_mat(j,k) = sigma_mat(j,k) - sigma_mat(j,6)*sigma_mat(k,6)/sigma_mat(6,6)
        IF( j .ne. k ) THEN
          sigma_beta_mat(k,j) = sigma_beta_mat(j,k)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

END SUBROUTINE make_sigma_mat
!------------------------------------------------------------------------

!+
!-
SUBROUTINE normal_sigma_mat(sigma_mat,normal)
  USE bmad
  USE eigen_mod
  
  IMPLICIT none

  REAL(rp) sigma_mat(1:6,1:6)
  REAL(rp) normal(1:3)
  complex(rp) eval(1:6)
  complex(rp) evec(1:6,1:6)
  REAL(rp) sigmaS(1:6,1:6)

  REAL(rp) S(1:6,1:6)
  real(rp), parameter :: S2(2,2)= reshape( [0.0d0, -1.0d0, 1.0d0, 0.0d0], [2,2] )

  INTEGER i

  LOGICAL error

  S = 0.0d0
  S(1:2,1:2) = S2
  S(3:4,3:4) = S2
  S(5:6,5:6) = S2

  sigmaS = MATMUL(sigma_mat,S)

  CALL mat_eigen(sigmaS,eval,evec,error)

  normal(1) = ABS(aimag(eval(1)))
  normal(2) = ABS(aimag(eval(3)))
  normal(3) = ABS(aimag(eval(5)))
END SUBROUTINE normal_sigma_mat

!+
!-
SUBROUTINE make_xyz_sigma_mat(spatial_dist,xyz_sigma_mat)
  USE bmad

  IMPLICIT none

  TYPE(coord_struct) spatial_dist(:)
  REAL(rp) xyz_sigma_mat(:,:)

  INTEGER j, k
  INTEGER ix1, ix2
  INTEGER Nparts

  Nparts = SIZE(spatial_dist)

  DO j=1,3
    DO k=j,3
      ix1 = (j*2)-1
      ix2 = (k*2)-1
      !Fortran * on two arrays is element-wise multiplication
      xyz_sigma_mat(j,k) = SUM(spatial_dist(:)%vec(ix1)*spatial_dist(:)%vec(ix2)) / Nparts
      IF( j .ne. k ) THEN
        xyz_sigma_mat(k,j) = xyz_sigma_mat(j,k)
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE make_xyz_sigma_mat
!------------------------------------------------------------------------

SUBROUTINE gaussian_dist_by_J (P0, ele, mode, coord)
  USE bmad

  IMPLICIT none

  REAL(rp) P0
  TYPE (normal_modes_struct) mode
  TYPE (coord_struct) coord(:)
  TYPE (ele_struct) ele

  REAL(rp) P
  REAL(rp) rn

  REAL(rp) J_a, J_b, J_z
  REAL(rp) theta_a, theta_b, theta_z
  REAL(rp) Delta_s, delta_p

  INTEGER i, N

  N=SIZE(coord)

  DO i=1,N
    !- Obtain J from inverse transform sampling.
    CALL RANDOM_NUMBER(rn)
    J_a = -mode%a%emittance * LOG(rn)
    CALL RANDOM_NUMBER(rn)
    J_b = -mode%b%emittance * LOG(rn)
    CALL RANDOM_NUMBER(rn)
    J_z = -mode%sig_z * mode%sigE_E * LOG(rn)

    !Uniformly distributed phases
    CALL RANDOM_NUMBER(rn)
    theta_a = rn*twopi
    CALL RANDOM_NUMBER(rn)
    theta_b = rn*twopi
    CALL RANDOM_NUMBER(rn)
    theta_z = rn*twopi

    Delta_s = SQRT(2.0d0 * J_z * ele%z%beta)* SIN(theta_z) 
    delta_p = SQRT(2.0d0 * J_z / ele%z%beta)* COS(theta_z)
    P = P0 * (1.0d0 + delta_p)

    coord(i)%vec(1) = SQRT(2.0d0 * J_a * ele%a%beta) * SIN(theta_a) + &
                         delta_p * ele%a%eta

    coord(i)%vec(2) = -SQRT(2.0d0 * J_a / ele%a%beta)*ele%a%alpha*SIN(theta_a) + SQRT(2.0d0 * J_a / ele%a%beta)*COS(theta_a) + &
                         delta_p * ele%a%etap
    coord(i)%vec(2) = coord(i)%vec(2) * P

    coord(i)%vec(3) = SQRT(2.0d0 * J_b * ele%b%beta) * SIN(theta_b) + &
                         delta_p * ele%b%eta

    coord(i)%vec(4) = -SQRT(2.0d0 * J_b / ele%b%beta)*ele%b%alpha*SIN(theta_b) + SQRT(2.0d0 * J_b / ele%b%beta)*COS(theta_b) + &
                         delta_p * ele%b%etap
    coord(i)%vec(4) = coord(i)%vec(4) * P

    coord(i)%vec(5) = Delta_s
    coord(i)%vec(6) = SQRT(P**2 - coord(i)%vec(2)**2 - coord(i)%vec(4)**2)
  ENDDO

END SUBROUTINE gaussian_dist_by_J

!---------------------------------------------------------------------

SUBROUTINE shuffle(N,shuffled)
  ! Durstenfeld's O(n) algorithm for randomizing an array of numbers
  INTEGER N
  INTEGER shuffled(:)

  INTEGER i
  INTEGER j, temp
  REAL(rp) rn

  DO i=1,N
    shuffled(i) = i
  ENDDO

  DO i=N,2,-1
    CALL RANDOM_NUMBER(rn)
    j = CEILING(rn*i)
    temp = shuffled(i)
    shuffled(i) = shuffled(j)
    shuffled(j) = temp
  ENDDO
END SUBROUTINE shuffle

!---------------------------------------------------------------------

!+
! Subroutine xyz_to_canonical(gamma,xyz_dist,canonical_dist)
!
! Subroutine uses the spatial z coordinate to drift particles into the canonical coordinates
! where z refers to beta*c*dt and pz refers to (P-P0)/P0
!
! Input:
!   xyz_vec(1:6)    -- real(rp): Coordinates where z refers to spatial coordinate
!
! Output:
!   canonical_vec(1:6) -- real(rp): Coordinates where z refers to temporal coordinate
!-

SUBROUTINE xyz_to_canonical(P0, xyz_dist, canonical_dist)
!input xyz_vec is (x(z), Px, y(z), Py, z, Pz)
!output canonical_vec is (x(0), Px/P0, y(0), Py/P0, beta*c*dt, (P-P0)/P0)

IMPLICIT none

REAL(rp) P0
TYPE(coord_struct) xyz_dist(:)
TYPE(coord_struct) canonical_dist(:)

INTEGER N, i
REAL(rp) distance
REAL(rp) rel_pc

N = SIZE(xyz_dist)

DO i=1,N
  distance = -xyz_dist(i)%vec(5)

  !Convert pz from Pz into (P-P0)/P0
  canonical_dist(i)%vec(6) = SQRT( (xyz_dist(i)%vec(6)/P0)**2 + (xyz_dist(i)%vec(2)/P0)**2 + (xyz_dist(i)%vec(4)/P0)**2 ) - 1.0d0
  rel_pc = 1.0d0 + canonical_dist(i)%vec(6)

  canonical_dist(i)%vec(1) = xyz_dist(i)%vec(1) + distance * xyz_dist(i)%vec(2)/P0/rel_pc
  canonical_dist(i)%vec(2) = xyz_dist(i)%vec(2)/P0
  canonical_dist(i)%vec(3) = xyz_dist(i)%vec(3) + distance * xyz_dist(i)%vec(4)/P0/rel_pc
  canonical_dist(i)%vec(4) = xyz_dist(i)%vec(4)/P0
  canonical_dist(i)%vec(5) = -distance * ( 1.0d0 + ((xyz_dist(i)%vec(2)/P0)**2 + (xyz_dist(i)%vec(4)/P0)**2) / (2.0d0 * rel_pc**2) )
ENDDO

END SUBROUTINE xyz_to_canonical

!+
! Subroutine canonical_to_xyz(P0, canonical_dist,xyz_dist)
!
! Subroutine uses the BMAD z coordinate to transform particles from the BMAD coordinate
! system to a real longitudinal spatial distribution.
! The particles are propagated through a drift for a time -z/beta/c_light
!
! Modules needed:
!   use precision_def
!
! Input:
!   P0                 -- real(rp): Momentum of reference particle
!   canonical_vec(1:6) -- real(rp): Coordinates where z refers to temporal coordinate
!
! Output:
!   xyz_vec(1:6)    -- real(rp): Coordinates where z refers to spatial coordinate
!-

SUBROUTINE canonical_to_xyz(P0, canonical_dist, xyz_dist)
!input canonical_vec is (x(0), Px/P0, y(0), Py/P0, beta*c*dt, (P-P0)/P0)
!output xyz_vec is (x(z), Px, y(z), Py, z, Pz)

IMPLICIT none

TYPE(coord_struct) canonical_dist(:)
TYPE(coord_struct) xyz_dist(:)

REAL(rp) P0
REAL(rp) beta
REAL(rp) g2b2  !gamma*gamma*beta*beta
REAL(rp) denom
REAL(rp) ps
REAL(rp) dx_dt, dy_dt, ds_dt
REAL(rp) delta_t, rel_pc
INTEGER N, i

N = SIZE(canonical_dist)

DO i=1,N
  !Convert pz from (P-P0)/P0 into Ps
  rel_pc = 1 + canonical_dist(i)%vec(6)
  ps = SQRT(rel_pc**2 - canonical_dist(i)%vec(2)**2 - canonical_dist(i)%vec(4)**2)
  xyz_dist(i)%vec(6) = ps * P0

  ! Calculate the relativistic beta of the particle
  beta = rel_pc*P0 / SQRT( m_electron**2 + (rel_pc*P0)**2 )
  g2b2 = beta*beta/(1.0d0-beta*beta)

  denom = SQRT( 1.0d0/g2b2 + canonical_dist(i)%vec(2)**2 + &
                             canonical_dist(i)%vec(4)**2 + &
                             ps**2 )

  ! True in paraxial approximation and ultra-relativistic approximation
  dx_dt = canonical_dist(i)%vec(2)/denom * c_light
  dy_dt = canonical_dist(i)%vec(4)/denom * c_light
  ds_dt =                       ps/denom * c_light

  delta_t = canonical_dist(i)%vec(5) / beta / c_light

  xyz_dist(i)%vec(1) = canonical_dist(i)%vec(1) + delta_t * dx_dt
  xyz_dist(i)%vec(2) = canonical_dist(i)%vec(2) * P0  
  xyz_dist(i)%vec(3) = canonical_dist(i)%vec(3) + delta_t * dy_dt
  xyz_dist(i)%vec(4) = canonical_dist(i)%vec(4) * P0
  xyz_dist(i)%vec(5) = delta_t * ds_dt 

ENDDO

END SUBROUTINE canonical_to_xyz

!-  !+
!-  !-
!-  
!-  SUBROUTINE boost_sigma_mat(Bx, By, Bz, dir, sigma_mat, boosted_sigma_mat)
!-    !To boost forward, dir = +1.  To boost backward, dir = -1.
!-    IMPLICIT NONE
!-  
!-    REAL(rp) Bx, By, Bz
!-    REAL(rp) Gfactor
!-    INTEGER dir
!-    REAL(rp) :: sigma_mat(6,6)
!-    REAL(rp) :: boosted_sigma_mat(6,6)
!-    TYPE(coord_struct) :: distp(:)
!-  
!-    REAL(rp) B2, Gdir
!-    REAL(rp) gamma_mat(1:4,1:4)
!-  
!-    INTEGER i
!-    REAL(rp) four_momentum(1:4)
!-    REAL(rp) temp(1:4)
!-    REAL(rp), PARAMETER :: m2 = m_electron*m_electron
!-    REAL(rp) BxByGfactor, BxBzGfactor, ByBzGfactor
!-    REAL(rp) ByByBzBzBxBxGdir, BxBxBzBzByByGdir, BxBxByByBzBzGdir
!-  
!-    B2 = Bx*Bx + By*By + Bz*Bz
!-    Gdir = SQRT(1-B2)**(-dir)
!-  
!-    CALL make_gamma_mat(dir*Bx, dir*By, dir*Bz, gamma_mat)
!-  
!-    Gfactor = (Gdir-1)
!-  
!-    BxByGfactor = Bx*By*Gfactor/B2
!-    BxBzGfactor = Bx*Bz*Gfactor/B2
!-    ByBzGfactor = By*Bz*Gfactor/B2
!-    ByByBzBzBxBxGdir = (By*By + Bz*Bz + Bx*Bx*Gdir)/B2
!-    BxBxBzBzByByGdir = (Bx*Bx + Bz*Bz + By*By*Gdir)/B2
!-    BxBxByByBzBzGdir = (Bx*Bx + By*By + Bz*Bz*Gdir)/B2
!-  
!-    DO i=1,6
!-      !To boost spatial coordinates, simply apply length contraction along {Bx, By, Bz}
!-      distp(i)%vec(1) = (dist(i)%vec(3)*BxByGfactor + dist(i)%vec(5)*BxBzGfactor + dist(i)%vec(1)*(ByByBzBzBxBxGdir))
!-      distp(i)%vec(3) = (dist(i)%vec(1)*BxByGfactor + dist(i)%vec(5)*ByBzGfactor + dist(i)%vec(3)*(BxBxBzBzByByGdir))
!-      distp(i)%vec(5) = (dist(i)%vec(1)*BxBzGfactor + dist(i)%vec(3)*ByBzGfactor + dist(i)%vec(5)*(BxBxByByBzBzGdir))
!-  
!-      !To boost momentum, use a lorentz boost
!-      four_momentum(1) = SQRT( dist(i)%vec(2)**2 + dist(i)%vec(4)**2 + dist(i)%vec(6)**2 + m2)
!-      four_momentum(2) = dist(i)%vec(2)
!-      four_momentum(3) = dist(i)%vec(4)
!-      four_momentum(4) = dist(i)%vec(6)
!-  
!-      temp = MATMUL(gamma_mat,four_momentum)
!-      distp(i)%vec(2) = temp(2)
!-      distp(i)%vec(4) = temp(3)
!-      distp(i)%vec(6) = temp(4)
!-    ENDDO
!-  
!-  END SUBROUTINE boost_sigma_mat

!+
!-

SUBROUTINE boost_dist(Bx, By, Bz, dir, dist, distp)
  !To boost forward, dir = +1.  To boost backward, dir = -1.
  IMPLICIT NONE

  REAL(rp) Bx, By, Bz
  REAL(rp) Gfactor
  INTEGER dir
  TYPE(coord_struct) :: dist(:)
  TYPE(coord_struct) :: distp(:)

  REAL(rp) B2, Gdir
  REAL(rp) gamma_mat(1:4,1:4)

  INTEGER i
  INTEGER Nmp
  REAL(rp) four_momentum(1:4)
  REAL(rp) temp(1:4)
  REAL(rp), PARAMETER :: m2 = m_electron*m_electron
  REAL(rp) BxByGfactor, BxBzGfactor, ByBzGfactor
  REAL(rp) ByByBzBzBxBxGdir, BxBxBzBzByByGdir, BxBxByByBzBzGdir

  B2 = Bx*Bx + By*By + Bz*Bz
  Gdir = SQRT(1-B2)**(-dir)

  CALL make_gamma_mat(dir*Bx, dir*By, dir*Bz, gamma_mat)

  Nmp = SIZE(dist)

  Gfactor = (Gdir-1)

  BxByGfactor = Bx*By*Gfactor/B2
  BxBzGfactor = Bx*Bz*Gfactor/B2
  ByBzGfactor = By*Bz*Gfactor/B2
  ByByBzBzBxBxGdir = (By*By + Bz*Bz + Bx*Bx*Gdir)/B2
  BxBxBzBzByByGdir = (Bx*Bx + Bz*Bz + By*By*Gdir)/B2
  BxBxByByBzBzGdir = (Bx*Bx + By*By + Bz*Bz*Gdir)/B2

  DO i=1,Nmp
    !To boost spatial coordinates, simply apply length contraction along {Bx, By, Bz}
    distp(i)%vec(1) = (dist(i)%vec(3)*BxByGfactor + dist(i)%vec(5)*BxBzGfactor + dist(i)%vec(1)*(ByByBzBzBxBxGdir))
    distp(i)%vec(3) = (dist(i)%vec(1)*BxByGfactor + dist(i)%vec(5)*ByBzGfactor + dist(i)%vec(3)*(BxBxBzBzByByGdir))
    distp(i)%vec(5) = (dist(i)%vec(1)*BxBzGfactor + dist(i)%vec(3)*ByBzGfactor + dist(i)%vec(5)*(BxBxByByBzBzGdir))

    !To boost momentum, use a lorentz boost
    four_momentum(1) = SQRT( dist(i)%vec(2)**2 + dist(i)%vec(4)**2 + dist(i)%vec(6)**2 + m2)
    four_momentum(2) = dist(i)%vec(2)
    four_momentum(3) = dist(i)%vec(4)
    four_momentum(4) = dist(i)%vec(6)

    temp = MATMUL(gamma_mat,four_momentum)
    distp(i)%vec(2) = temp(2)
    distp(i)%vec(4) = temp(3)
    distp(i)%vec(6) = temp(4)
  ENDDO

END SUBROUTINE boost_dist

!+
!-

SUBROUTINE make_gamma_mat(Bx, By, Bz, gamma_mat)
  IMPLICIT NONE

  REAL(rp) Bx, By, Bz
  REAL(rp) gamma_mat(1:4,1:4)

  REAL(rp) G, B2

  B2 = Bx*Bx + By*By + Bz*Bz
  G = 1.0d0/SQRT(1.0d0-B2)

  gamma_mat(1,1) = G
  gamma_mat(1,2) = -G * Bx
  gamma_mat(1,3) = -G * By
  gamma_mat(1,4) = -G * Bz
  
  gamma_mat(2,1) = gamma_mat(1,2)
  gamma_mat(2,2) = 1.0d0 + (G-1.0d0)*Bx*Bx/B2
  gamma_mat(2,3) = (G-1.0d0)*Bx*By/B2
  gamma_mat(2,4) = (G-1.0d0)*Bx*Bz/B2

  gamma_mat(3,1) = gamma_mat(1,3)
  gamma_mat(3,2) = gamma_mat(2,3)
  gamma_mat(3,3) = 1.0d0 + (G-1.0d0)*By*By/B2
  gamma_mat(3,4) = (G-1.0d0)*By*Bz/B2

  gamma_mat(4,1) = gamma_mat(1,4)
  gamma_mat(4,2) = gamma_mat(2,4)
  gamma_mat(4,3) = gamma_mat(3,4)
  gamma_mat(4,4) = 1.0d0 + (G-1.0d0)*Bz*Bz/B2
END SUBROUTINE make_gamma_mat

!-----------------------------------------------------------------------------

SUBROUTINE collide(coords_a, coords_b, delta_t, sigma_y, tau_rest_frame, local_density, change_px, change_py, change_pz)
  IMPLICIT NONE

  TYPE(coord_struct) coords_a
  TYPE(coord_struct) coords_b
  REAL(rp) delta_t
  REAL(rp) sigma_y
  REAL(rp) tau_rest_frame     !damping rate
  REAL(rp) local_density
  REAL(rp) change_px
  REAL(rp) change_py
  REAL(rp) change_pz
  REAL(rp), PARAMETER :: m_electron2 = m_electron*m_electron

  REAL(rp) va_x, va_y, va_z
  REAL(rp) vb_x, vb_y, vb_z
  REAL(rp) u_x,  u_y,  u_z
  REAL(rp) u, uperp

  REAL(rp) phi, theta
  REAL(rp) variance
  REAL(rp) azimuth, zenith
  REAL(rp) rnum
  REAL(rp) delta

  REAL(rp) change_ux
  REAL(rp) change_uy
  REAL(rp) change_uz

  REAL(rp) Ea_c, Eb_c

  REAL(rp) :: conversion_factor = e_charge/ c_light / me_kgs

  !Calculate particle velocity in units of meters per second
  Ea_c = SQRT(coords_a%vec(2)**2 + coords_a%vec(4)**2 + coords_a%vec(6)**2 + m_electron2) / c_light
  Eb_c = SQRT(coords_b%vec(2)**2 + coords_b%vec(4)**2 + coords_b%vec(6)**2 + m_electron2) / c_light
  va_x = coords_a%vec(2) / Ea_c
  va_y = coords_a%vec(4) / Ea_c
  va_z = coords_a%vec(6) / Ea_c
  vb_x = coords_b%vec(2) / Eb_c
  vb_y = coords_b%vec(4) / Eb_c
  vb_z = coords_b%vec(6) / Eb_c

  !Calculate relative particle velocity
  u_x = va_x - vb_x
  u_y = va_y - vb_y
  u_z = va_z - vb_z

  u = SQRT(u_x**2 + u_y**2 + u_z**2)
  uperp = SQRT(u_x**2 + u_y**2)
  phi = ATAN(u_y/u_x)
  theta = ATAN(uperp/u_z)

  !Calculate variance of delta, which is used to calculate scattering angle
  CALL calc_variance(local_density,tau_rest_frame,sigma_y,u,delta_t,variance)

  CALL RANDOM_NUMBER(rnum)
  azimuth = twopi * rnum
  !Box-Muller returns a normally distributed random number with unit variance and expectation zero.
  CALL box_muller(delta)
  delta = delta*variance
  zenith = 2.0_rp * ATAN(delta)

  !Calculate change in relative velocity
  change_ux = u_x/uperp*u_z*SIN(zenith)*COS(azimuth) - u_y/uperp*u*SIN(zenith)*SIN(azimuth) - u_x*(1.0_rp-COS(zenith))
  change_uy = u_y/uperp*u_z*SIN(zenith)*COS(azimuth) + u_x/uperp*u*SIN(zenith)*SIN(azimuth) - u_y*(1.0_rp-COS(zenith))
  change_uz = -uperp*SIN(zenith)*COS(azimuth) - u_z*(1.0_rp-COS(zenith))

  !Convert to change in relative momentum in eV/c units
  change_px = change_ux / conversion_factor
  change_py = change_uy / conversion_factor
  change_pz = change_uz / conversion_factor
END SUBROUTINE collide

!-----------------------------------------------------------------------------

SUBROUTINE box_muller(z1,z2)
  IMPLICIT NONE

  REAL(rp) z1
  REAL(rp), OPTIONAL :: z2
  REAL(rp) u1, u2

  CALL RANDOM_NUMBER(u1)
  CALL RANDOM_NUMBER(u2)

  z1 = SQRT(-2.0_rp*LOG(u1))*COS(twopi*u2)
  IF( PRESENT(z2) ) THEN
    z2 = SQRT(-2.0_rp*LOG(u1))*SIN(twopi*u2)
  ENDIF

END SUBROUTINE box_muller

!------------------------------------------------------------------------------

SUBROUTINE subtract_center(canonical_dist, adjusted_dist, dump_vals)
  USE bmad

  IMPLICIT NONE

  TYPE(coord_struct) canonical_dist(:)
  TYPE(coord_struct) adjusted_dist(:) 

  LOGICAL, OPTIONAL :: dump_vals

  REAL(rp) avg_x, avg_xp, avg_y, avg_yp, avg_z, avg_zp
  INTEGER Nmp

  Nmp = size(canonical_dist)
  
  ! move distribution of particles back to the origin for fitting
  avg_x  = SUM(canonical_dist(:)%vec(1))/Nmp
  avg_xp = SUM(canonical_dist(:)%vec(2))/Nmp
  avg_y  = SUM(canonical_dist(:)%vec(3))/Nmp
  avg_yp = SUM(canonical_dist(:)%vec(4))/Nmp
  avg_z  = SUM(canonical_dist(:)%vec(5))/Nmp
  avg_zp = SUM(canonical_dist(:)%vec(6))/Nmp

!  IF( PRESENT(dump_vals) ) THEN 
!    WRITE(9876,'(6ES14.4)') avg_x, avg_xp, avg_y, avg_yp, avg_z, avg_zp
!  ENDIF

  adjusted_dist(:)%vec(1) = canonical_dist(:)%vec(1) - avg_x
  adjusted_dist(:)%vec(2) = canonical_dist(:)%vec(2) - avg_xp
  adjusted_dist(:)%vec(3) = canonical_dist(:)%vec(3) - avg_y
  adjusted_dist(:)%vec(4) = canonical_dist(:)%vec(4) - avg_yp
  adjusted_dist(:)%vec(5) = canonical_dist(:)%vec(5) - avg_z
  adjusted_dist(:)%vec(6) = canonical_dist(:)%vec(6) - avg_zp
END SUBROUTINE

!------------------------------------------------------------------------------

END MODULE mc_ibs_mod


