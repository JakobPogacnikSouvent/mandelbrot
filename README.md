# About
This is a project for the class High Performance Computing, of the M2 Data Science program at Université de Lille.
The project contains c programs that compute the mandelbrot set image in four different ways
1. `mandel.c` - sequentially
2. `mpi_mandel.c` - parallel computing using MPI and naive static workload allocation
3. `mpi_master.c` - parallel computing using MPI and dynamic master-worker workload allocation
4. `mandel.cu` - parellel computing using CUDA

The files `mpi_master_no_count.c` and `mandel_no_count.cu` are similar to `mpi_master.c` and `mandel.cu` but do not
calculate the average number of iteration per pixel performed. Withouth the added workload they output they output more
precise timings of the algorithms and are used to calculate perfomance improvements.

The programs return the image `mandel.ras` in the raster image format.

# Quickstart
To download and run the project follow the following instructions.

## Download the files
```
git clone https://github.com/JakobPogacnikSouvent/mandelbrot.git
```

## Compile the code
For `mandel.c`
```
gcc -o mandel mandel.c
```

For `mpi_mandel.c`, `mpi_master.c` and `mpi_master_no_count.c`
```
mpicc -o mpi_master mpi_master.c -lm
```

For `mandel.cu` and `mandel_no_count.cu`
```
nvcc  -o mandel_cuda mandel.cu --generate-code arch=compute_60,code=sm_60 -O3
```

## Run the code
For `mandel`
```
./mandel 800 800 -2 -2 2 2 10000
```

For `mpi_mandel`, `mpi_master` and `mpi_master_no_count`
```
oarsub -q default -l cluster=1/core=4 'mpirun --mca pml ^ucx -machinefile $OAR_NODEFILE ./mpi_master 800 800 -2 -2 2 2 10000'
```

For `mandel_cuda` and `mandel_cuda_no_count`
```
./mandel_cuda 800 800 -2 -2 2 2 10000 16 16
```

## Input arguments
The programs expect the following inputs
```
mandel dimx dimy xmin ymin xmax ymax depth block_x block_y
```

| Argument           | Description                                   |
|-------------------|-----------------------------------------------|
| `dimx, dimy`       | Width and height of the generated image      |
| `xmin, ymin`       | Minimum values of the computation domain in the complex plane |
| `xmax, ymax`       | Maximum values of the computation domain in the complex plane |
| `depth`            | Maximum number of iterations for the Mandelbrot computation |
| `block_x, block_y`  | Kernel block size (CUDA only) |

