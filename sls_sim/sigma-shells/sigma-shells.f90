program sigma_shells

implicit none

character(100) filename

integer, parameter :: rp = 8

integer ia,ja
integer f
integer Jparts
integer N, Nsig

real(rp) Cx, Cy, Cz
real(rp) Jx, Jy, Jz
real(rp) emit_x, emit_y, emit_z
real(rp) p(8,6)
real(rp) long_beta, long_alpha, long_gamma
real(rp) px,pxp,py,pyp,pz,pzp
real(rp) beta_x, alpha_x, beta_y, alpha_y
real(rp) eta_x, etap_x, eta_y, etap_y
real(rp) co_x, co_xp, co_y, co_yp, co_z, co_zp

Jparts = 21  !odd
Nsig=7

beta_x = 14.81938088
alpha_x = -0.13890670
beta_y = 5.70138064
alpha_y = -0.41533232
eta_x = 0.11628635
etap_x = 0.0
eta_y = 0.0
etap_y = 0.0
long_beta = 6.0
long_alpha = 1.3111E-03
long_gamma = 1.6666E-01

emit_z =  4.36669636842e-6  ! bunch_length * energy_spread
emit_x =  1.80598E-09
emit_y =  1.80598E-10

!Produce surface distributions
do N=1, Nsig
  write(filename,'(a,i1,a)') "sigma_shell_",N,".0sig.raw"
  open(55,file=filename)
  do ia=1,Jparts
    do ja=1,ia
      Cx = 1 - (ia-1.0)/(Jparts-1.0)
      Cy = 1 - Cx - (ja-1.0)/(Jparts-1.0)
      Cz = 1 - Cx - Cy

      call chop(Cx,1.0d-6)
      call chop(Cy,1.0d-6)
      call chop(Cz,1.0d-6)

      Jx = Cx * emit_x * N
      Jy = Cy * emit_y * N
      Jz = Cz * emit_z * N

      px =  sqrt(2*Jx*beta_x)
      pxp = -alpha_x*sqrt(2*Jx/beta_x)
      py =  sqrt(2*Jy*beta_y)
      pyp = -alpha_y*sqrt(2*Jy/beta_y)
      pz  = sqrt(2*Jz*long_beta)
      pzp = -long_alpha*sqrt(2*Jz/long_beta)

      p(1,1)=  px; p(1,2)=  pxp; p(1,3)=  py; p(1,4)=  pyp; p(1,5)=  pz; p(1,6)=  pzp;
      p(2,1)=  px; p(2,2)=  pxp; p(2,3)= -py; p(2,4)= -pyp; p(2,5)=  pz; p(2,6)=  pzp;
      p(3,1)= -px; p(3,2)= -pxp; p(3,3)=  py; p(3,4)=  pyp; p(3,5)=  pz; p(3,6)=  pzp;
      p(4,1)= -px; p(4,2)= -pxp; p(4,3)= -py; p(4,4)= -pyp; p(4,5)=  pz; p(4,6)=  pzp;

      p(5,1)=  px; p(5,2)=  pxp; p(5,3)=  py; p(5,4)=  pyp; p(5,5)= -pz; p(5,6)= -pzp;
      p(6,1)=  px; p(6,2)=  pxp; p(6,3)= -py; p(6,4)= -pyp; p(6,5)= -pz; p(6,6)= -pzp;
      p(7,1)= -px; p(7,2)= -pxp; p(7,3)=  py; p(7,4)=  pyp; p(7,5)= -pz; p(7,6)= -pzp;
      p(8,1)= -px; p(8,2)= -pxp; p(8,3)= -py; p(8,4)= -pyp; p(8,5)= -pz; p(8,6)= -pzp;

      do f=1,8
        write(55,'(6es12.3,3f13.4)') co_x  + p(f,6)*eta_x  + p(f,1), &
                                     co_xp + p(f,6)*etap_x + p(f,2), &
                                     co_y  + p(f,6)*eta_y  + p(f,3), &
                                     co_yp + p(f,6)*etap_y + p(f,4), &
                                     co_z  + p(f,5), &
                                     co_zp + p(f,6), Cx, Cy, Cz
      enddo
    enddo
  enddo
  close(55)
enddo

contains
  subroutine chop(x,eps)
    implicit none
    real(rp) x,eps

    if(abs(x) < eps) x = 0
  end subroutine
end program


