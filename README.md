# About
This is a project for the class High Performance Computing, of the M2 Data Science program at Université de Lille.
The project contains c programs that compute the mandelbrot set image in three different ways
1. `mandel.c` - sequentially
2. `mpi_mandel.c` - parallel computing using MPI and naive static workload allocation
3. `mpi_master.c` - parallel computing using MPI and dynamic master-worker workload allocation

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

For `mpi_mandel.c` and `mpi_master.c`
```
mpicc -o mpi_master mpi_master.c -lm
```

## Run the code
For `mandel`
```
./mandel 800 800 -2 -2 2 2 10000
```

For `mpi_mandel` and `mpi_master`
```
oarsub -q default -l cluster=1/core=4 'mpirun --mca pml ^ucx -machinefile $OAR_NODEFILE ./mpi_master 800 800 -2 -2 2 2 10000'
```

## Input arguments
The programs expect the following inputs
```
mandel dimx dimy xmin ymin xmax ymax depth
```

| Argument           | Description                                   |
|-------------------|-----------------------------------------------|
| `dimx, dimy`       | Width and height of the generated image      |
| `xmin, ymin`       | Minimum values of the computation domain in the complex plane |
| `xmax, ymax`       | Maximum values of the computation domain in the complex plane |
| `depth`            | Maximum number of iterations for the Mandelbrot computation |

