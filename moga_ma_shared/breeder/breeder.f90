program breeder

use pisa_mod

implicit none

type(breeder_params_struct) breeder_params
type(pool_struct), allocatable :: parents(:), offspring(:)

integer i

logical err_flag

namelist / breeder_nl / breeder_params

open(unit = 10, file = 'breeder.in')
read(10, nml = breeder_nl)
close(10)

allocate(parents(800))
do i=1,size(parents)
  allocate(parents(i)%x(19))
enddo
allocate(offspring(3200))
do i=1,size(offspring)
  allocate(offspring(i)%x(19))
enddo

!call read_initial_population(pool, n_dep, alpha, n_vars, seed_file, err_flag)
call read_initial_population(parents, 0, 800, 19, 'moga_seeds.dat.uniq', err_flag)
if(err_flag) then
  write(*,*) "Error reading initial population.  aborting."
  stop
endif

call kangal_breeder_bare(parents, offspring, breeder_params)

call write_pool(offspring, 'moga_seeds.bred', 11)

end program






