/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 1
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and CNRS / INRIA / University of Pau
 ! (c) Princeton University / California Institute of Technology and CNRS / INRIA / University of Pau
 !                             July 2012
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <string.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#define MAX(a, b) (((a) > (b)) ? (a) : (b))


void save_to_max_surface_file_(float* maxval) {
  int rank;
  char filename[BUFSIZ];
  FILE* fp;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else
  rank = 0;
#endif
  sprintf(filename,"maxval_surface_proc_%03d.dat",rank);
  fp = fopen(filename,"a+");
  fprintf(fp,"%e\n",*maxval);
  fclose(fp);
}


void save_fvector_(float* vector, int* size, int* id, int* cpu_or_gpu) {
  FILE* fp;
  char filename[BUFSIZ];
  if(*cpu_or_gpu == 0) {
    sprintf(filename, "debug_output_cpu_%d.dat",*id);
  }
  else {
    sprintf(filename, "debug_output_gpu_%d.dat",*id);
  }
  fp = fopen(filename, "wb");
  printf("writing vector, vector[0]=%e\n",vector[0]);
  fwrite(vector, sizeof(float), *size, fp);
  fclose(fp);

}

void save_ivector_(int* vector, int* size, int* id, int* cpu_or_gpu) {
  FILE* fp;
  char filename[BUFSIZ];
  if(*cpu_or_gpu == 0) {
    sprintf(filename, "debug_output_cpu_%d.dat",*id);
  }
  else {
    sprintf(filename, "debug_output_gpu_%d.dat",*id);
  }
  fp = fopen(filename, "wb");
  fwrite(vector, sizeof(int), *size, fp);
  fclose(fp);

}


void get_max_from_surface_file_(int* nodes_per_iterationf,int* NSTEP) {
  int nodes_per_iteration = *nodes_per_iterationf;
  char filename[BUFSIZ];
  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  sprintf(filename,"/scratch/eiger/rietmann/SPECFEM3D_AIGLE/OUTPUT_FILES/DATABASES_MPI/proc%06d_surface_movie",procid);

  FILE* fp; int it;
  printf("Opening %s for analysis\n",filename);
  fp = fopen(filename,"rb");
  //char* errorstr;
  if(fp == 0) {
    //errorstr = (char*) strerror(errno);
    printf("FILE ERROR:%s\n",(char*) strerror(errno));
    perror("file error\n");
    exit(1);
  }

  float* vector = (float*)malloc(nodes_per_iteration*sizeof(float));
  float max_val;
  int i;
  max_val = 0.0;
  for(it=0;it<*NSTEP;it++) {
    int pos = (sizeof(float)*nodes_per_iteration)*(it);
    fseek(fp,pos,SEEK_SET);
    fread(vector,sizeof(float),nodes_per_iteration,fp);
    for(i=0;i<nodes_per_iteration;i++) {
      max_val = MAX(max_val,vector[i]);
    }
    if(it % 500 == 0) {
      printf("scanning it=%d\n",it);
    }
  }
  printf("max_val=%e\n",max_val);
}

void compare_two_vectors_exact_(int* sizef,float* vector1,float* vector2,int* num_errors) {

  int size = *sizef;
  int i;
  int error_count = 0;

  for(i=0;i<size;++i) {
    if(vector1[i] != vector2[i]) {
      error_count++;
      if(error_count < 10) {
  printf("err[%d]: %e != %e\n",i,vector1[i],vector2[i]);
      }
    }
  }
  printf("**** Error Count: %d ****\n",error_count);
  *num_errors = error_count;
}

void compare_two_vectors_(int* sizef,float* vector1,float* vector2,int* num_errors) {

  int size = *sizef;
  int i;
  int error_count = 0;
  for(i=0;i<size;++i) {
    if(vector1[i] != 0) {
      if( fabsf(vector1[i]-vector2[i])/vector1[i] > 0.01) {
  if(fabsf(vector1[i]-vector2[i]) > 1e-20) {
  error_count++;
      if(error_count<10) {
        printf("err[%d]: %e != %e\n",i,vector1[i],vector2[i]);
      }
      }
      }
    }
    /* if(vector1[i] != vector2[i]) { */
    /*   if(fabsf(vector1[i]-vector2[i]) > 1e-25) { */
    /*  error_count++; */
    /*  if(error_count<50) { */
    /*    printf("err[%d]: %e != %e\n",i,vector1[i],vector2[i]); */
    /*  } */
    /*   } */
    /* } */
  }
  printf("**** Error Count: %d ****\n",error_count);
  *num_errors = error_count;
}

void compare_surface_files_(int* bytes_per_iteration, int* number_of_iterations) {

  char* cpu_file = "/scratch/eiger/rietmann/SPECFEM3D/OUTPUT_FILES/DATABASES_MPI/cpu_proc000001_surface_movie";
  char* gpu_file = "/scratch/eiger/rietmann/SPECFEM3D/OUTPUT_FILES/DATABASES_MPI/cpu_v2_proc000001_surface_movie";

  FILE* fp_cpu;
  fp_cpu = fopen(cpu_file,"rb");
  //char* errorstr;
  if(fp_cpu == 0) {
    //errorstr = (char*) strerror(errno);
    //printf("CPU FILE ERROR:%s\n",errorstr);
    printf("CPU FILE ERROR:%s\n",(char*) strerror(errno));
    perror("cpu file error\n");
  }
  FILE* fp_gpu;
  fp_gpu = fopen(gpu_file,"rb");

  if(fp_gpu == NULL) {
    //errorstr = (char*) strerror(errno);
    //printf("GPU FILE ERROR:%s\n",errorstr);
    printf("GPU FILE ERROR:%s\n",(char*) strerror(errno));
    perror("gpu file error\n");
  }

  /* pause_for_debug(); */

  float* gpu_vector = (float*)malloc(*bytes_per_iteration);
  float* cpu_vector = (float*)malloc(*bytes_per_iteration);
  int i,it,error_count=0;
  for(it=0;it<*number_of_iterations;it++) {
    int pos = (*bytes_per_iteration)*(it);

    fseek(fp_cpu,pos,SEEK_SET);
    fseek(fp_gpu,pos,SEEK_SET);

    int number_of_nodes = *bytes_per_iteration/sizeof(float);
    fread(cpu_vector,sizeof(float),number_of_nodes,fp_cpu);
    fread(gpu_vector,sizeof(float),number_of_nodes,fp_gpu);
    int size = number_of_nodes;
    float gpu_min_val=10;
    float gpu_max_val=0;
    float cpu_min_val=10;
    float cpu_max_val=0;
    if(it<100) {
      for(i=0;i<size;i++) {
  if((fabs(cpu_vector[i] - gpu_vector[i])/(fabs(cpu_vector[i])+1e-31) > 0.01)) {
    if(error_count < 30) printf("ERROR[%d]: %g != %g\n",i,cpu_vector[i], gpu_vector[i]);
    if(cpu_vector[i] > 1e-30) error_count++;
  }
  if(gpu_vector[i]>gpu_max_val) gpu_max_val = gpu_vector[i];
  if(gpu_vector[i]<gpu_min_val) gpu_min_val = gpu_vector[i];
  if(cpu_vector[i]>cpu_max_val) cpu_max_val = cpu_vector[i];
  if(cpu_vector[i]<cpu_min_val) cpu_min_val = cpu_vector[i];
      }
      printf("%d Total Errors\n",error_count);
      printf("size:%d\n",size);
      printf("GPU:[min/max]=%e/%e\n",gpu_min_val,gpu_max_val);
      printf("CPU:[min/max]=%e/%e\n",cpu_min_val,cpu_max_val);
    }
  }
  printf("End of Surface Compare\n");
  exit(1);
}


void compare_fvector_(float* vector, int* size, int* id, int* cpu_or_gpu) {
  FILE* fp;
  char cmp_filename[BUFSIZ];
  float* compare_vector = (float*)malloc(*size*sizeof(float));
  if(*cpu_or_gpu == 0) { //swap gpu/cpu for compare
    sprintf(cmp_filename, "debug_output_gpu_%d.dat",*id);
  }
  else {
    sprintf(cmp_filename, "debug_output_cpu_%d.dat",*id);
  }
  fopen(cmp_filename, "rb");
  /* read the values */
  if((fp=fopen(cmp_filename, "rb"))==NULL) {
    printf("Cannot open comparison file %s.\n",cmp_filename);
    exit(1);
  }
  if(fread(compare_vector, sizeof(float), *size, fp) != *size) {
    if(feof(fp))
       printf("Premature end of file.");
    else
       printf("File read error.");
  }

  fclose(fp);

  int i;
  int error_count=0;
  for(i=0;i<*size;i++) {
    if((fabs(vector[i] - compare_vector[i])/vector[i] > 0.0001)) {
      if(error_count < 30) {
        printf("ERROR[%d]: %f != %f\n",i,compare_vector[i], vector[i]);
      }
      error_count++;
      /* if(compare_vector[i] > 1e-30) error_count++; */
    }
  }
  printf("%d Total Errors\n",error_count);
  printf("size:%d\n",*size);
  /* for(i=0;i<30;i++) { */
  /*   printf("val[%d]: %g != %g\n",i,compare_vector[i], vector[i]); */
  /*   /\* printf("error_check[%d]= %g\n",abs(vector[i] - compare_vector[i])/vector[i]); *\/ */
  /* } */
}

void compare_ivector_(int* vector, int* size, int* id, int* cpu_or_gpu) {
  FILE* fp;
  char cmp_filename[BUFSIZ];
  int* compare_vector = (int*)malloc(*size*sizeof(int));
  if(*cpu_or_gpu == 0) { //swap gpu/cpu for compare
    sprintf(cmp_filename, "debug_output_gpu_%d.dat",*id);
  }
  else {
    sprintf(cmp_filename, "debug_output_cpu_%d.dat",*id);
  }
  fopen(cmp_filename, "rb");
  /* read the values */
  if((fp=fopen(cmp_filename, "rb"))==NULL) {
    printf("Cannot open comparison file %s.\n",cmp_filename);
    exit(1);
  }
  if(fread(compare_vector, sizeof(int), *size, fp) != *size) {
    if(feof(fp))
       printf("Premature end of file.");
    else
       printf("File read error.");
  }

  fclose(fp);

  int i;
  int error_count=0;
  for(i=0;i<*size;i++) {
    if((abs(vector[i] - compare_vector[i])/vector[i] > 0.01) && error_count < 30) {
      printf("ERROR[%d]: %d != %d\n",i,compare_vector[i], vector[i]);
      error_count++;
    }
  }
  printf("%d Total Errors\n",error_count);
}
