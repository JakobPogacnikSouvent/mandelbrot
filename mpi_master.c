/*
 * M2 Data Science 
 * Computation of the Mandelbrot set; sequential version
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>	/* timings */
#include <string.h>     /* for memset */
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#include "rasterfile.h"



char info[] = "\
Usage:\n\
      mandel dimx dimy xmin ymin xmax ymax depth\n\
\n\
      dimx,dimy: sizes of the image to be generated\n\
      xmin,ymin,xmax,ymax: computation domain in the complex plan\n\
      depth: maximum number of iterations\n\
\n\
Some execution examples:\n\
      mandel 800 800 0.35 0.355 0.353 0.358 200\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 500\n\
      mandel 800 800 -0.736 -0.184 -0.735 -0.183 300\n\
      mandel 800 800 -1.48478 0.00006 -1.48440 0.00044 100\n\
      mandel 800 800 -1.5 -0.1 -1.3 0.1 10000\n\
";



double my_gettimeofday(){
  struct timeval tmp_time;
  gettimeofday(&tmp_time, NULL);
  return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}




/**
 * Conversion from one LINUX integer (4 bytes) into one SUN integer
 * @param i integer to convert
 * @return converted integer
 */

int swap(int i) {
  int init = i; 
  int conv;
  unsigned char *o, *d;
	  
  o = ( (unsigned char *) &init) + 3; 
  d = (unsigned char *) &conv;
  
  *d++ = *o--;
  *d++ = *o--;
  *d++ = *o--;
  *d++ = *o--;
  
  return conv;
}


/*** 
 * By Francois-Xavier MOREL (M2 SAR, oct2009): 
 */

unsigned char power_component(int i, int p) {
  unsigned char o;
  double iD=(double) i;

  iD/=255.0;
  iD=pow(iD,p);
  iD*=255;
  o=(unsigned char) iD;
  return o;
}

unsigned char cos_component(int i, double freq) {
  unsigned char o;
  double iD=(double) i;
  iD=cos(iD/255.0*2*M_PI*freq);
  iD+=1;
  iD*=128;
  
  o=(unsigned char) iD;
  return o;
}

/*** 
 * Coloring choice: define (only) one of the constants below: 
 */
//#define ORIGINAL_COLOR
#define COS_COLOR 

#ifdef ORIGINAL_COLOR
#define RED_COMPONENT(i)    ((i)/2)
#define GREEN_COMPONENT(i)     ((i)%190)
#define BLUE_COMPONENT(i)     (((i)%120) * 2)
#endif /* #ifdef ORIGINAL_COLOR */
#ifdef COS_COLOR
#define RED_COMPONENT(i)    cos_component(i,13.0)
#define GREEN_COMPONENT(i)     cos_component(i,5.0)
#define BLUE_COMPONENT(i)     cos_component(i+10,7.0)
#endif /* #ifdef COS_COLOR */


/**
 *  Save the data array into the 8-bit Rasterfile format 
 *  with a palette of 256 gray levels, from white (0 value)
 *  to black (255)
 *    @param name Image name 
 *    @param width Image width
 *    @param height Image height
 *    @param p pointer the the buffer containing the image
 */

void save_rasterfile( char *name, int width, int height, unsigned char *p) {
  FILE *fd;
  struct rasterfile file;
  int i;
  unsigned char o;

  if ( (fd=fopen(name, "w")) == NULL ) {
	printf("Error when creating the file %s \n",name);
	exit(1);
  }

  file.ras_magic  = swap(RAS_MAGIC);	
  file.ras_width  = swap(width);	  /* image width (in pixels) */
  file.ras_height = swap(height);         /* image height (in pixels) */
  file.ras_depth  = swap(8);	          /* depth of each pixel (1, 8 or 24 )   */
  file.ras_length = swap(width*height);   /* image size (number of bytes)	*/
  file.ras_type    = swap(RT_STANDARD);	  /* file type */
  file.ras_maptype = swap(RMT_EQUAL_RGB);
  file.ras_maplength = swap(256*3);

  fwrite(&file, sizeof(struct rasterfile), 1, fd); 
  
  /* Color palette: red component */
  i = 256;
  while( i--) {
    o = RED_COMPONENT(i);
    fwrite( &o, sizeof(unsigned char), 1, fd);
  }

  /* Color palette: green component */
  i = 256;
  while( i--) {
    o = GREEN_COMPONENT(i);
    fwrite( &o, sizeof(unsigned char), 1, fd);
  }

  /* Color palette: blue component */
  i = 256;
  while( i--) {
    o = BLUE_COMPONENT(i);
    fwrite( &o, sizeof(unsigned char), 1, fd);
  }

  // To check the line order in the image: 
  //fwrite( p, width*height/3, sizeof(unsigned char), fd);
  
  // To see the '0' color:
  // memset (p, 0, width*height);
  
  fwrite( p, width*height, sizeof(unsigned char), fd);
  fclose( fd);
}

/**
 * Given the coordinates of a point \f$c=a+ib\f$ in the complex plan,
 * this function returns the corresponding color by estimating 
 * the distance of the point to the Mandelbrot set.
 * Let the following complex sequence: 
 * \f[
 * \left\{\begin{array}{l}
 * z_0 = 0 \\
 * z_{n+1} = z_n^2 + c
 * \end{array}\right.
 * \f]
 * the number of iterations required for the sequence to diverge is the 
 * number \f$ n \f$ for which \f$ |z_n| > 2 \f$. 
 * This number is reduced to a value between 0 and 255, hence corresponding 
 * to a color in the palette. 
 */

unsigned char xy2color(double a, double b, int prof) {
  double x, y, temp, x2, y2;
  int i;

  x = y = 0.;
  for( i=0; i<prof; i++) {
    /* saving the previous value of x (which will be overwritten) */
    temp = x;
    /* new values for x and y */
    x2 = x*x;
    y2 = y*y;
    x = x2 - y2 + a;
    y = 2*temp*y + b;
    if( x2 + y2 > 4.0) break;
  }
  return (i==prof)?255:(int)((i%255)); 
}

// Generate k-th row of the Mandelbrot image
void generate_row(int k, int w, double xmin, double ymin, double xinc, double yinc, int prof, unsigned char *row) {
    double x, y;

    // Compute the y-coordinate of the k-th row
    y = ymin + k * yinc;

    for (int j = 0; j < w; j++) {
        x = xmin + j * xinc;
        row[j] = xy2color(x, y, prof);
    }
}

/* 
 * Main part: for each grid point, run xy2color()
 */

int main(int argc, char *argv[]) {
  /* Computation domain in the complex plan */
  double xmin, ymin;
  double xmax, ymax;
  /* Image sizes */
  int w,h;
  /* Pixel sizes */
  double xinc, yinc;
  /* Iteration depth */
  int prof;
  /* Intermediate variables */
  double x, y;
  /* Timing: */
  double debut, fin;
  
  int rank, p;
  MPI_Status status;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Timing start */
  debut = my_gettimeofday();

  if(argc == 1) fprintf( stderr, "%s\n", info);
  if (p < 2) {
    fprintf(stderr, "Need at least 2 processes.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // TODO: Print error if values don't divide evenly
  // TODO: Print error if 1 process, warning if 2
  
  /* Default values */
  xmin = -2; ymin = -2;
  xmax =  2; ymax =  2;
  w = h = 800;
  prof = 10000;
  
  /* Get command line parameters */
  if(argc > 1) w    = atoi(argv[1]);
  if(argc > 2) h    = atoi(argv[2]);
  if(argc > 3) xmin = atof(argv[3]);
  if(argc > 4) ymin = atof(argv[4]);
  if(argc > 5) xmax = atof(argv[5]);
  if(argc > 6) ymax = atof(argv[6]);
  if(argc > 7) prof = atoi(argv[7]);

  /* Computing the pixel sizes */
  xinc = (xmax - xmin) / (w-1);
  yinc = (ymax - ymin) / (h-1);
  
  /* Parameter display */
  fprintf( stderr, "Domain: {[%lg,%lg]x[%lg,%lg]}\n", xmin, ymin, xmax, ymax);
  fprintf( stderr, "Pixel sizes: %lg %lg\n", xinc, yinc);
  fprintf( stderr, "Depth: %d\n",  prof);
  fprintf( stderr, "Image sizes: %dx%d\n", w, h);
  
  if (rank == p-1) {
    // Master = p-1

    unsigned char *mandel = malloc(w*h*sizeof(unsigned char));

    if(mandel == NULL) {
      fprintf( stderr, "Array memory allocation error \n");
    }

    // Process p1 is working on row worker_row[p1]
    int worker_row[p-1];

    for (int worker_ind = 0; worker_ind < p-1; worker_ind++) {
      worker_row[worker_ind] = worker_ind;
    }

    int next_row = p-1;
    int running_procceses = p-1;
    while (running_procceses > 0) {
      // MPI_probe
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,  &status);
      int sender = status.MPI_SOURCE;

      int receiving_row = worker_row[sender];

      // MPI_receive
      MPI_Recv(&mandel[w * receiving_row], w, MPI_UNSIGNED_CHAR, sender, 0, MPI_COMM_WORLD, &status);

      // MPI_send next row
      if (next_row < h) {
        MPI_Send(&next_row, 1, MPI_INT, sender, 0, MPI_COMM_WORLD);
        worker_row[sender] = next_row;
        next_row++;
      } else {
        int stop = -1;
        MPI_Send(&stop, 1, MPI_INT, sender, 0, MPI_COMM_WORLD);
        running_procceses--;

        fprintf( stderr, "Stoping process: %d, %d processes remaining\n", sender, running_procceses);
      }

    }

    save_rasterfile( "mandel.ras", w, h, mandel);
    free(mandel);

  } else {
    // Workers
    unsigned char *row = malloc( w*sizeof(unsigned char));

    if (row == NULL) {
        fprintf(stderr, "Memory allocation error \n");
    }
    
    // First row is always the same so hardcode it into workers
    generate_row(rank, w, xmin, ymin, xinc, yinc, prof, row);
    MPI_Send(row, w, MPI_UNSIGNED_CHAR, p-1, 0, MPI_COMM_WORLD);

    while (1) {
      int row_index;

      MPI_Recv(&row_index, 1, MPI_INT, p-1, 0, MPI_COMM_WORLD, &status);

      // If -1 then no more rows
      if (row_index == -1){
        fprintf( stderr, "Process: %d, stopping\n", rank);
        break;
      }

      generate_row(row_index, w, xmin, ymin, xinc, yinc, prof, row);

      MPI_Send(row, w, MPI_UNSIGNED_CHAR, p-1, 0, MPI_COMM_WORLD);
    }    
    free(row);
  }
  
  /* Timing stop */
  fin = my_gettimeofday();
  fprintf( stderr, "Total computation time for process %d: %g sec\n", rank, fin - debut);
  fprintf( stdout, "%g\n", fin - debut);

  MPI_Finalize();
  return 0;
}

