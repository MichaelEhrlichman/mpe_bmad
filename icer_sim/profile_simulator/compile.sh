#nvfortran -cuda profile_simulator_gpu.f90 -o profile_simulator_gpu 
nvfortran -cuda test_gpu_code.f90 -o test_gpu_code -cudalib=curand
