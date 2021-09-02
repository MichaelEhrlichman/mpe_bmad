program bayes

  use fgsl

implicit none

real, parameter :: pi = 3.141592653

real, parameter :: a = 2.0
real, parameter :: s = 5.0
real, parameter :: a0 = 0.0
real, parameter :: s0 = 3.0

real, parameter :: va_x1 = 0.0
real, parameter :: va_x2 = 5.0
real, parameter :: vs_x1 = 0.01
real, parameter :: vs_x2 = 80.0
integer, parameter :: nmodel = 100

real, parameter :: meas_x1 = 0.0
real, parameter :: meas_x2 = 10.0
integer, parameter :: nmeas = 20.0

integer, parameter :: niter = 10000

integer i,j

real t(nmeas)
real meas(nmeas)
real posterior(nmodel,nmodel)
real lik(nmodel)
real va(nmodel), vs(nmodel)
real pr(nmodel)
real po(nmodel)
real a_s(niter), v_s(niter)
real dbug

type(fgsl_rng_type) :: trng
type(fgsl_rng) :: rng
real gammar

trng = fgsl_rng_default
rng = fgsl_rng_alloc(trng)

call random_seed()

call make_regressor()
call take_meas()
call make_va()
call make_vs()
call make_posterior()
call likelihood()
call prior()
call make_posterior2()
call gibbs_sampling()

!do i=1,100000
!  dbug = gammarnd(5.0,1.0)
!  write(1010,*) dbug
!enddo
!stop !FOO

open(35,name='meas.dat')
do i=1,nmeas
  write(35,*) t(i), meas(i)
enddo
close(35)

open(45,name='lik.dat')
do i=1,nmodel
  write(45,*) va(i), lik(i)
enddo
close(45)

open(55,name='pr.dat')
do i=1,nmodel
  write(55,*) va(i), pr(i)
enddo
close(55)

open(65,name='po.dat')
do i=1,nmodel
  write(65,*) va(i), po(i)
enddo
close(65)

open(75,name='posterior.dat')
do i=1,nmodel
  do j=1,nmodel
    write(75,'(2es14.4,es18.5e3)') va(j), vs(i), posterior(i,j)
  enddo
enddo
close(75)

open(85,name='gibbs.dat')
do i=1,niter
    write(85,'(2es14.4)') a_s(i), v_s(i)
enddo
close(85)

contains

subroutine gibbs_sampling()
  real n1, n2
  real mu, sigma
  integer nnn

  a_s(1) = 2.0d0
  v_s(1) = 10.0d0

  nnn = 1

  do i=2,niter
    call box_muller(n1,n2)
    mu = sum(t*meas)/sum(t**2)
    sigma = sqrt(v_s(i-1)/sum(t**2))
    a_s(i) = n1*sigma+mu
    gammar = fgsl_ran_gamma(rng, nmeas/2.0d0, 2.0d0/sum( (meas-a_s(i-1.0d0)*t)**2))
    v_s(i) = 1.0d0/gammar
  enddo
end subroutine

subroutine make_regressor()
  real dx
  integer i

  !generate regressor
  dx = (meas_x2-meas_x1)/(nmeas-1.0d0)
  do i=1,nmeas
    t(i) = meas_x1 + dx*(i-1.0d0)
  enddo
end subroutine

subroutine take_meas()
  real dx
  real n1, n2
  integer i

  !generate noisy data
  do i=1,nmeas,2
    call box_muller(n1,n2)
    meas(i) = a*t(i) + s*n1
    meas(i+1) = a*t(i+1) + s*n2
  enddo
end subroutine

subroutine make_posterior()
  integer i, j
  real S(nmeas)
  real li1, pr1

  do i=1,nmodel
    do j=1,nmodel
      S = va(j)*t
      li1 = exp( -nmeas/2.0 * log(vs(i)) - sum( (meas-S)**2/2.0/vs(i)) )
      pr1 = 1.0/vs(i)* normpdf(va(i),0.0,1000.0)
      posterior(i,j) = li1*pr1
    enddo
  enddo
  posterior = posterior / sum(posterior)
end subroutine

subroutine likelihood()
  lik = exp(-0.5*sum((spread(meas,2,nmodel)-outer_product(t,va))**2,1)/(s**2))
  lik = lik/sum(lik)
end subroutine

subroutine make_va()
  real dx
  integer i

  dx = (va_x2-va_x1)/(nmodel-1)

  do i=1,nmodel
    va(i) = va_x1 + dx*(i-1)
  enddo
end subroutine

subroutine make_vs()
  real dx
  integer i

  dx = (vs_x2-vs_x1)/(nmodel-1)

  do i=1,nmodel
    vs(i) = vs_x1 + dx*(i-1)
  enddo
end subroutine

subroutine box_muller(x1,x2)
  real x1, x2
  real u1, u2

  call random_number(u1)
  call random_number(u2)

  x1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
  x2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
end subroutine

subroutine prior()
  pr = exp((-0.5*(va-a0)**2)/(s0**2))
  pr = pr / sum(pr)
end subroutine

subroutine make_posterior2()
  real beta, beta0
  real sp, ap

  beta = 1.0/(s**2)
  beta0 = 1.0/(s0**2)

  sp = 1.0/sqrt(beta0+beta*dot_product(t,t))
  ap = ( beta0*a0 + beta*dot_product(t,meas) ) / (beta0 + beta*dot_product(t,t))
  po = normpdf(va,ap,sp)
  po = po/sum(po)
end subroutine

function mean(y) result(m)
  real y(:), m
  m = sum(y)/size(y)
end function

elemental real function normpdf(x,mu,sig)
  real, intent(in) :: x, mu, sig
  normpdf = exp(-((x-mu)**2)/2.0/sig/sig)/sig/sqrt(2.0*pi)
end function

function outer_product(a,b) result(c)
  real a(:), b(:)
  real c(size(a),size(b))
  integer i,j
  forall (i=1:size(a))
      forall(j=1:size(b)) c(i,j) = a(i)*b(j)
  end forall
end function

end program
