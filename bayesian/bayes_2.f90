program bayes

implicit none

real, parameter :: pi = 3.141592653

real, parameter :: a = 5.0
real, parameter :: s = 3.0
real, parameter :: a0 = 0.0
real, parameter :: s0 = 3.0

integer, parameter :: nmeas = 20
integer, parameter :: nmodel = 100
real, parameter :: model_x1 = -15.0
real, parameter :: model_x2 = 15.0
real, parameter :: meas_x1 = -1.0
real, parameter :: meas_x2 = 1.0

integer i

real x(nmeas)
real meas(nmeas)
real lik(nmodel)
real va(nmodel)
real pr(nmodel)
real po(nmodel)

call take_meas()
call eval_model()
call likelihood()
call prior()
call posterior()

open(35,name='meas.dat')
do i=1,nmeas
  write(35,*) meas(i), 0
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

contains

subroutine likelihood()
  lik = exp(-0.5*sum((spread(meas,2,nmodel)-outer_product(x,va))**2,1)/(s**2))
  lik = lik/sum(lik)
end subroutine

subroutine eval_model()
  real dx
  integer i

  dx = (model_x2-model_x1)/(nmodel-1)

  do i=1,nmodel
    va(i) = model_x1 + dx*(i-1)
  enddo
end subroutine

subroutine take_meas()
  real dx
  real n1, n2
  integer i

  dx = (meas_x2-meas_x2)/(nmeas-1)
  do i=1,nmeas
    x(i) = meas_x1 + dx*(i-1)
  enddo

  do i=1,nmeas,2
    call box_muller(n1,n2)
    meas(i) = a*x(i) + s*n1
    meas(i+1) = a*x(i+1) + s*n2
  enddo
end subroutine

subroutine box_muller(x1,x2)
  real x1, x2
  real u1, u2

  call random_seed()
  call random_number(u1)
  call random_number(u2)

  x1 = sqrt(-2.0d0*log(u1))*cos(2.0d0*pi*u2)
  x2 = sqrt(-2.0d0*log(u1))*sin(2.0d0*pi*u2)
end subroutine

subroutine prior()
  pr = exp((-0.5*(va-a0)**2)/(s0**2))
  pr = pr / sum(pr)
end subroutine

subroutine posterior()
  real beta, beta0
  real sp, ap

  beta = 1.0/(s**2)
  beta0 = 1.0/(s0**2)

  sp = 1.0/sqrt(beta0+beta*dot_product(x,x))
  ap = ( beta0*a0 + beta*dot_product(x,meas) ) / (beta0 + beta*dot_product(x,x))
  po = normpdf(va,ap,sp)
  po = po/sum(po)
end subroutine

function mean(y) result(m)
  real y(:), m
  m = sum(y)/size(y)
end function

function normpdf(x,mu,sig) result(y)
  real x(:), mu, sig
  real y(size(x))

  y = exp(-((x-mu)**2)/2.0/sig/sig)/sig/sqrt(2.0*pi)
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
