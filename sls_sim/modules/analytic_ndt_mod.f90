module analytic_ndt_mod

use bmad

implicit none

type ndt_struct
  complex(rp) h20001
  complex(rp) h00201
  complex(rp) h10002
  complex(rp) h21000
  complex(rp) h30000
  complex(rp) h10110
  complex(rp) h10020
  complex(rp) h10200

  complex(rp) h31000
  complex(rp) h40000
  complex(rp) h20110
  complex(rp) h11200
  complex(rp) h20020
  complex(rp) h20200
  complex(rp) h00310
  complex(rp) h00400
end type

contains

subroutine print_ndts(ndts)

implicit none

type(ndt_struct) ndts

write(*,*) "h20001 = ", abs(ndts%h20001)
write(*,*) "h00201 = ", abs(ndts%h00201)
write(*,*) "h10002 = ", abs(ndts%h10002)
write(*,*) "h21000 = ", abs(ndts%h21000)
write(*,*) "h30000 = ", abs(ndts%h30000)
write(*,*) "h10110 = ", abs(ndts%h10110)
write(*,*) "h10020 = ", abs(ndts%h10020)
write(*,*) "h10200 = ", abs(ndts%h10200)

write(*,*) "h31000 = ", abs(ndts%h31000)
write(*,*) "h40000 = ", abs(ndts%h40000)
write(*,*) "h20110 = ", abs(ndts%h20110)
write(*,*) "h11200 = ", abs(ndts%h11200)
write(*,*) "h20020 = ", abs(ndts%h20020)
write(*,*) "h20200 = ", abs(ndts%h20200)
write(*,*) "h00310 = ", abs(ndts%h00310)
write(*,*) "h00400 = ", abs(ndts%h00400)

end subroutine

subroutine analytic_ndt(lat,n)

use bmad

implicit none

type(lat_struct) lat
type(ndt_struct) n
integer i,j
real(rp) pxi, pyi, pxj, pyj
real(rp) bxi, byi, bxj, byj

n%h20001 = 0
n%h00201 = 0
n%h10002 = 0
n%h21000 = 0
n%h30000 = 0
n%h10110 = 0
n%h10020 = 0
n%h10200 = 0

n%h31000 = 0
n%h40000 = 0
n%h20110 = 0
n%h11200 = 0
n%h20020 = 0
n%h20200 = 0
n%h00310 = 0
n%h00400 = 0
DO i=1,lat%n_ele_track
  pxi = lat%ele(i)%a%phi
  pyi = lat%ele(i)%b%phi
  bxi = lat%ele(i)%a%beta
  byi = lat%ele(i)%b%beta
  n%h20001   = n%h20001   + lat%ele(i)%value(l$)*(lat%ele(i)%value(k1$) - 2.0*lat%ele(i)%value(k2$)*lat%ele(i)%a%eta)*bxi * EXP( CMPLX(0.0d0,2.0d0*pxi) )
  n%h00201   = n%h00201   + lat%ele(i)%value(l$)*(lat%ele(i)%value(k1$) - 2.0*lat%ele(i)%value(k2$)*lat%ele(i)%a%eta)*byi * EXP( CMPLX(0.0d0,2.0d0*pyi) )
  n%h10002   = n%h10002   + lat%ele(i)%value(l$)*(lat%ele(i)%value(k1$) - lat%ele(i)%value(k2$)*lat%ele(i)%a%eta)*lat%ele(i)%a%eta*sqrt(bxi) * EXP( CMPLX(0.0d0,1.0d0*pxi) )
  n%h21000   = n%h21000   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*bxi**(3./2.) * EXP( CMPLX(0.0d0,1.0d0*pxi) )
  n%h30000   = n%h30000   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*bxi**(3./2.) * EXP( CMPLX(0.0d0,3.0d0*pxi) )
  n%h10110   = n%h10110   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*bxi**(1./2.)*byi * EXP( CMPLX(0.0d0,1.0d0*pxi) )
  n%h10020   = n%h10020   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*bxi**(1./2.)*byi * EXP( CMPLX(0.0d0,pxi-2.0*pyi) )
  n%h10200   = n%h10200   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*bxi**(1./2.)*byi * EXP( CMPLX(0.0d0,pxi+2.0*pyi) )
  DO j=1,lat%n_ele_track
    pxj = lat%ele(j)%a%phi
    pyj = lat%ele(j)%b%phi
    bxj = lat%ele(j)%a%beta
    byj = lat%ele(j)%b%beta
    IF (j .ne. i) THEN
      n%h31000 = n%h31000 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        bxi**(3./2.)*bxj**(3./2.) * &
                        exp( cmplx(0.0d0,3.0*pxi-pxj) )
      n%h40000 = n%h40000 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        bxi**(3./2.)*bxj**(3./2.) * &
                        exp( cmplx(0.0d0,3.0*pxi+pxj) )
      n%h20110 = n%h20110 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        bxi**(1./2.)*bxj**(1./2.)*byi * &
                        ( bxj*exp( -cmplx(0.0d0,pxi-3.0*pxj) ) - exp( cmplx(0.0d0,pxi+pxj) ) + &
                        2.0*byj*exp( cmplx(0.0d0,pxi+pxj+2.0*pyi-2.0*pyj )) )
      n%h11200 = n%h11200 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        bxi**(1./2.)*bxj**(1./2.)*byi * &
                        ( bxj*(exp( -cmplx(0.0d0,pxi-pxj-2.0*pyi) ) - exp( cmplx(0.0d0,pxi-pxj+2.0*pyi)) ) + &
                        2.0*byj*(exp( cmplx(0.0d0,pxi-pxj+2.0*pyi)) + exp( -cmplx(0.0d0,pxi-pxj-2.0*pyi)) ) )
      n%h20020 = n%h20020 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        sqrt(bxi*bxj)*byi * &
                        ( bxj*exp(-cmplx(0.0d0,pxi-3.0*pxj+2.0*pyi)) - (bxj+4.0*byj)*exp(cmplx(0.0d0,pxi+pxj-2.0*pyi)) )
      n%h20200 = n%h20200 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        sqrt(bxi*bxj)*byi * &
                        ( bxj*exp(-cmplx(0.0d0,pxi-3.0*pxj-2.0*pyi)) - (bxj-4.0*byj)*exp(cmplx(0.0d0,pxi+pxj+2.0*pyi)) )
      n%h00310 = n%h00310 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        sqrt(bxi*bxj)*byi*byj * &
                        (exp(cmplx(0.0d0,pxi-pxj+2.0*pyi))-exp(-cmplx(0.0d0,pxi-pxj-2.0d0*pyi)))
      n%h00400 = n%h00400 + sign(1.0,1.0*(j-i)) * &
                        lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(j)%value(k2$)*lat%ele(j)%value(l$) * &
                        sqrt(bxi*bxj)*byi*byj * &
                        exp(cmplx(0.0d0,pxi-pxj+2.0*pyi+2.0*pyj))
    ENDIF
  ENDDO
ENDDO
n%h20001 = n%h20001 / 8.0
n%h00201 = -1.0 * n%h00201 / 8.0
n%h10002 = n%h10002 / 2.0
n%h21000 = -1.0 * n%h21000 / 8.0
n%h30000 = -1.0 * n%h30000 / 24.0
n%h10110 =  n%h10110 / 4.0
n%h10020 =  n%h10020 / 8.0
n%h10200 =  n%h10200 / 8.0

n%h31000 =  n%h31000 / 32.0_rp * CMPLX(0.0,1.0_rp)
n%h40000 =  n%h40000 / 64.0_rp * CMPLX(0.0,1.0_rp)
n%h20110 =  n%h20110 / 32.0_rp * CMPLX(0.0,1.0_rp)
n%h11200 =  n%h11200 / 32.0_rp * CMPLX(0.0,1.0_rp)
n%h20020 =  n%h20020 / 64.0_rp * CMPLX(0.0,1.0_rp)
n%h20200 =  n%h20200 / 64.0_rp * CMPLX(0.0,1.0_rp)
n%h00310 =  n%h00310 / 32.0_rp * CMPLX(0.0,1.0_rp)
n%h00400 =  n%h00400 / 32.0_rp * CMPLX(0.0,1.0_rp)

end subroutine

end module












