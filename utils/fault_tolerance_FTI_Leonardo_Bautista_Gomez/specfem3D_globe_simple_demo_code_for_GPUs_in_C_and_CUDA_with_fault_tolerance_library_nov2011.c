
//
// This is the main program
//

// standard includes
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

// include for CUDA
#include <cuda.h>
#include <cuda_runtime.h>

// path to the header file created by the mesher (which contains the number of
// mesh elements, of grid points etc.)
#include "../DATABASES_FOR_SOLVER/values_from_mesher_C.h"

// include the CUDA kernels
#include "specfem3D_kernels.h"

// include our own header
#include "specfem3D.h"

// standard header for MPI (do NOT move above because USE_MPI is set or not in "specfem3D.h")
#ifdef USE_MPI
  #include <mpi.h>
#endif

#include "fti.h"

#define DIV 4


#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
// read 2-D addressing for summation between slices: indirection arrays used to fill the MPI buffers
static void read_arrays_buffers_solver(char* prname, int* iboolleft_xi, int* iboolright_xi, int* iboolleft_eta, int* iboolright_eta,int* p_npoin2D_xi,int* p_npoin2D_eta);

// get indirection from faces to edges of mesh slices
static void get_indirect_addressing_1D_buffers(int* iboolleft_xi, int* iboolright_xi, int* iboolleft_eta, int* iboolright_eta,int npoin2D_xi,int npoin2D_eta, int* ibool1D_leftxi_lefteta, int* ibool1D_leftxi_righteta, int* ibool1D_rightxi_lefteta, int* ibool1D_rightxi_righteta);

// subroutine to assemble the elastic forces at each time step for the mesh slices based on non-blocking MPI.
//
// It returns to the calling routine if no communications have arrived yet.
//
// the assembling process is finished when *phase == 4
//
// we keep launching that subroutine on the host in parallel to the execution of the
// CUDA kernels that compute the inner elements located inside each mesh slice
// until the assembly process is finished.
//
// Launched in a polling way by the host after computation of all the outer elements of a given mesh slice
// while the GPU is executing kernel 2 on inner elements.
static void assemble_MPI(int iproc_xi, int iproc_eta, int addressing[NPROC_XI][NPROC_ETA], int npoin2D_xi, int npoin2D_eta, int* phase,
    MPI_Request* request_send_xi_to_left, MPI_Request* request_recv_xi_from_left,
    MPI_Request* request_send_xi_to_right, MPI_Request* request_recv_xi_from_right,
    MPI_Request* request_send_eta_to_up, MPI_Request* request_recv_eta_from_up,
    MPI_Request* request_send_eta_to_down, MPI_Request* request_recv_eta_from_down,
    float* buffer_send_xi_to_left, float* buffer_recv_xi_from_left,
    float* buffer_send_xi_to_right, float* buffer_recv_xi_from_right,
    float* buffer_send_eta_to_up, float* buffer_recv_eta_from_up,
    float* buffer_send_eta_to_down, float* buffer_recv_eta_from_down,
    int* ibool1D_leftxi_lefteta, int* ibool1D_leftxi_righteta, int* ibool1D_rightxi_lefteta, int* ibool1D_rightxi_righteta, float* accel_all_MPI_buffers);

// mass matrix assembly routine, using blocking MPI:
// subroutine to assemble the mass matrix at each time step for the mesh slices based on blocking MPI.
// this is done once and for all before the time loop therefore it is
// sufficient to use blocking MPI
static void assemble_MPI_blocking_scalar(int* iboolleft_xi, int* iboolright_xi, int* iboolleft_eta, int* iboolright_eta, int iproc_xi, int iproc_eta, int addressing[NPROC_XI][NPROC_ETA], int npoin2D_xi, int npoin2D_eta, float* buffer_send_xi, float* buffer_recv_xi, float* buffer_send_eta, float* buffer_recv_eta, float rmass[NGLOB]);
#endif

static void print_CUDA_error_if_any(cudaError_t,int);

// returns a handle to a compute device
int dev;

// call a Fortran routine to read the unformatted binary data files created by the Fortran mesher
//// DK DK 33333333333333 now in Fortran
  extern void read_arrays_solver_(float xix[NSPEC*NGLL3_ALIGN],float xiy[NSPEC*NGLL3_ALIGN],float xiz[NSPEC*NGLL3_ALIGN],float etax[NSPEC*NGLL3_ALIGN],float etay[NSPEC*NGLL3_ALIGN],float etaz[NSPEC*NGLL3_ALIGN],float gammax[NSPEC*NGLL3_ALIGN],float gammay[NSPEC*NGLL3_ALIGN],float gammaz[NSPEC*NGLL3_ALIGN],float kappav[NSPEC*NGLL3_ALIGN],float muv[NSPEC*NGLL3_ALIGN],int ibool[NSPEC*NGLL3_ALIGN],float rmass_inverse[NGLOB], int* myrank,
#ifdef DISPLAY_POSITION_SOURCE_RECEIVER
  double xstore[NSPEC][NGLLZ][NGLLY][NGLLX], double ystore[NSPEC][NGLLZ][NGLLY][NGLLX], double zstore[NSPEC][NGLLZ][NGLLY][NGLLX],
#else
  double xstore[1][1][1][1], double ystore[1][1][1][1], double zstore[1][1][1][1],
#endif
  int* read_mesh_coordinate_arrays);

int main(int argc, char *argv[])
{

        int nb_color_tot;
        int NPROCTOT;

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
        int we_are_done_computing_outer_elements_first;
#endif
        int nb_color[2];
        char prname[200];
        char filename[300];

        int myrank;
#ifdef USE_MPI
        int nbprocs;
#endif

        // For checkpointing and performance counters
        double t1, t2, t3, t4, t5, t6, t7, t8;
        t1 = MPI_Wtime();
        FILE *fd;
        char ckptfile[50];

// initialize MPI or start a run on only one GPU with no MPI
#ifdef USE_MPI
        MPI_Init(&argc,&argv);

        // FTI initialization
        FTI_Init("config.fti");

        MPI_Comm_size(FTI_COMM_WORLD,&nbprocs);
        MPI_Comm_rank(FTI_COMM_WORLD,&myrank);
#else
        myrank = 0;
#endif

// detection of the number of GPU cards
        int deviceCount;
        cudaGetDeviceCount(&deviceCount);
        if (myrank == 0) {
          printf("\nYou have %d GPU card(s) per PCIe bus (or similar)\n",deviceCount);
          printf("(note that there can be several PCIe buses on a TESLA node)\n\n"); }

// assign the right GPU card
        cudaSetDevice(myrank % deviceCount);
//#ifdef USE_MPI
//        MPI_Barrier(FTI_COMM_WORLD);
//#endif
//        int running_on_what_device;
//        cudaGetDevice(&running_on_what_device);
//        printf("hello, this is proc %d using GPU card %d\n\n",myrank,running_on_what_device);
//        fflush(stdout);
//#ifdef USE_MPI
//        MPI_Barrier(FTI_COMM_WORLD);
//#endif

// returns a handle to a compute device
        cuDeviceGet(&dev, 0);

// get device properties
        struct cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp,dev);

// exit if the machine has no CUDA-enabled device
        if (deviceCount <= 0 || (deviceProp.major == 9999 && deviceProp.minor == 9999))
              { fprintf(stderr,"No CUDA-enabled device found, exiting...\n\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1); }

        if (myrank == 0) {
/* from driver_types.h , here is what can be displayed:
  char   name[256];                 ///< ASCII string identifying device
  size_t totalGlobalMem;            ///< Global memory available on device in bytes
  size_t sharedMemPerBlock;         ///< Shared memory available per block in bytes
  int    regsPerBlock;              ///< 32-bit registers available per block
  int    warpSize;                  ///< Warp size in threads
  size_t memPitch;                  ///< Maximum pitch in bytes allowed by memory copies
  int    maxThreadsPerBlock;        ///< Maximum number of threads per block
  int    maxThreadsDim[3];          ///< Maximum size of each dimension of a block
  int    maxGridSize[3];            ///< Maximum size of each dimension of a grid
  int    clockRate;                 ///< Clock frequency in kilohertz
  size_t totalConstMem;             ///< Constant memory available on device in bytes
  int    major;                     ///< Major compute capability
  int    minor;                     ///< Minor compute capability
  size_t textureAlignment;          ///< Alignment requirement for textures
  int    deviceOverlap;             ///< Device can concurrently copy memory and execute a kernel
  int    multiProcessorCount;       ///< Number of multiprocessors on device
  int    kernelExecTimeoutEnabled;  ///< Specified whether there is a run time limit on kernels
  int    integrated;                ///< Device is integrated as opposed to discrete
  int    canMapHostMemory;          ///< Device can map host memory with cudaHostAlloc/cudaHostGetDevicePointer
  int    computeMode;               ///< Compute mode (See ::cudaComputeMode)  */


// display device properties
          printf("Device Name = %s\n",deviceProp.name);
          printf("multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
          printf("totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
          printf("totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
          printf("sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
// This will display for instance:
//        Maximum number of threads per block:           512
//        Maximum sizes of each dimension of a block:    512 x 512 x 64
//        Maximum sizes of each dimension of a grid:     65535 x 65535 x 1
//
//        Thererefore the largest possible number of threads that can be launched on a device is:
//        65K*65K  * 512*512*64  * 512, which is a HUGE number of threads.
          printf("Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
          printf("Maximum size of each dimension of a block: %d x %d x %d\n",deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
          printf("Maximum sizes of each dimension of a grid: %d x %d x %d\n",deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);

          printf("Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
#ifdef PRINT_canMapHostMemory
          if(deviceProp.canMapHostMemory)
            { printf("canMapHostMemory: TRUE\n"); }
          else
            { printf("canMapHostMemory: FALSE\n"); }

          if(deviceProp.deviceOverlap)
            { printf("deviceOverlap: TRUE\n"); }
          else
            { printf("deviceOverlap: FALSE\n"); }
#endif

// make sure that the device has compute capability >= 1.3
          if (deviceProp.major < 1)
                { fprintf(stderr,"Compute capability major number should be at least 1, exiting...\n\n");
#ifdef USE_MPI
                  MPI_Abort(FTI_COMM_WORLD,1);
#endif
                  exit(1); }
#ifdef STOP_IF_NOT_CAPABILITY_13
          if (deviceProp.major == 1 && deviceProp.minor < 3)
                { fprintf(stderr,"Compute capability should be at least 1.3, exiting...\n");
                  fprintf(stderr,"You can comment out #define STOP_IF_NOT_CAPABILITY_13 in specfem3D.h and recompile with something lower than sm_13 in the Makefile...\n\n");
#ifdef USE_MPI
                  MPI_Abort(FTI_COMM_WORLD,1);
#endif
                  exit(1); }
#endif

          if (deviceProp.major == 1) {
            if (deviceProp.minor <= 1) {
              printf("The number of registers per multiprocessor is 8192.\nThe maximum number of active threads per multiprocessor is 768.\n\n");
            }
            else if (deviceProp.minor <= 3) {
              printf("The number of registers per multiprocessor is 16384.\nThe maximum number of active threads per multiprocessor is 1024.\n\n");
            }
          }
          }

// see if we can really test the CUDA 2.2 'zero copy' feature or not
#ifdef USE_ZERO_COPY
#ifndef ACTUALLY_ASSEMBLE_MPI_SLICES
        fprintf(stderr,"cannot test USE_ZERO_COPY when ACTUALLY_ASSEMBLE_MPI_SLICES is off, exiting...\n\n");
#ifdef USE_MPI
        MPI_Abort(FTI_COMM_WORLD,1);
#endif
        exit(1);
#endif
#endif

// prefix for the mesh files that need to be read
        sprintf(prname,"../DATABASES_FOR_SOLVER/proc%06d_reg1_",myrank);

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
        int threads_buf;
        int grid_buf_xi,grid_buf_eta,grid_buf_xi_eta;

/////////////        clock_t waiting_beg, waiting_end;
/////////////        float waiting_tot = 0.f;

        int assembling_phase;  //////////////////////  ,wait_flag;
        int addressing[NPROC_XI][NPROC_ETA];
        int iproc_xi_slice[NPROC_XI*NPROC_ETA];
        int iproc_eta_slice[NPROC_XI*NPROC_ETA];
        int iproc_num,iproc_xi,iproc_eta,npoin2D_xi = 0,npoin2D_eta = 0;

        int* ibool1D_leftxi_lefteta;
        int* ibool1D_leftxi_righteta;
        int* ibool1D_rightxi_lefteta;
        int* ibool1D_rightxi_righteta;

        float *buffer_send_xi_to_left, *buffer_recv_xi_from_left;
        float *buffer_send_xi_to_right, *buffer_recv_xi_from_right;
        float *buffer_send_eta_to_up, *buffer_recv_eta_from_up;
        float *buffer_send_eta_to_down, *buffer_recv_eta_from_down;

        MPI_Request request_send_xi_to_left, request_recv_xi_from_left;
        MPI_Request request_send_xi_to_right, request_recv_xi_from_right;
        MPI_Request request_send_eta_to_up, request_recv_eta_from_up;
        MPI_Request request_send_eta_to_down, request_recv_eta_from_down;

        float* accel_all_MPI_buffers;

        float* d_accel_all_MPI_buffers;

        int iboolleft_xi[NGLOB2DMAX_XMIN_XMAX];
        int iboolright_xi[NGLOB2DMAX_XMIN_XMAX];
        int iboolleft_eta[NGLOB2DMAX_YMIN_YMAX];
        int iboolright_eta[NGLOB2DMAX_YMIN_YMAX];

        int* d_iboolleft_xi;
        int* d_iboolright_xi;
        int* d_iboolleft_eta;
        int* d_iboolright_eta;
#endif

// approximate density of the geophysical medium in which the source is located
// this value is only a constant scaling factor therefore it does not really matter
#define rho 4500.f

        static float rmass_inverse[NGLOB];

        static float hprime_xx[NGLL2];
        static float hprimewgll_xx[NGLL2];
        static float wgllwgll_xy[NGLL2];
        static float wgllwgll_xz[NGLL2];
        static float wgllwgll_yz[NGLL2];

#ifdef WRITE_SEISMOGRAM
// record a seismogram to check that the simulation went well
        static float seismogram_x[NSTEP];
        static float seismogram_y[NSTEP];
        static float seismogram_z[NSTEP];
#endif

        static int ispec,iglob,i,j,k,it;
        int *number_of_elements_in_this_color;
        static int nb_elem_color;

// to compute the maximum of the norm of the displacement vector, to make sure that
// the simulation remains stable
        static float Usolidnorm_local,Usolidnorm_global,current_value;
        static float displ[NDIM][NGLOB];
        static float veloc[NDIM][NGLOB];

        static float time;
        clock_t timeloop_begin;
/////////////        clock_t k2_beg, k2_end;
/////////////        float k2_tot;
///////////////////////////////////////        float k2_prev = 0.f;
        float timeloop_total;
        static float source_accel;

// to read external files
  FILE *IIN;

  float memory_size;

#ifdef USE_MPI
#ifdef USE_SERIAL_CASCADE_FOR_IOs
  MPI_Status status;
  int you_can_go_ahead_and_start_reading;
#endif
#endif

  if (myrank == 0) {
    printf("\nNSPEC = %d\n",NSPEC);
    printf("NGLOB = %d\n\n",NGLOB);
    printf("NSTEP = %d\n",NSTEP/DIV);
    printf("deltat = %f\n\n",deltat);

// estimate total memory size (the size of a real number is 4 bytes)
// we perform the calculation in single precision rather than integer
// to avoid integer overflow in the case of very large meshes
    memory_size = 4.f * ((3.f*NDIM + 1.f) * NGLOB + 12.f * (float)(NGLLX*NGLLY*NGLLZ)*(float)(NSPEC));
    printf("approximate total memory size used on each GPU = %f Mb\n\n",memory_size/1024.f/1024.f);
    }

// make sure that we can use the Deville et al. (2002) routines
  if(NGLLX != 5 || NGLLY != 5 || NGLLZ != 5) {
         fprintf(stderr,"we must have NGLLX = NGLLY = NGLLZ = 5 to be able to use the CUDA kernels, exiting...\n");
#ifdef USE_MPI
         MPI_Abort(FTI_COMM_WORLD,1);
#endif
         exit(1);
       }

// make sure that we can use MPI
#ifdef USE_MPI
  if(NPROC_XI < 2 || NPROC_ETA < 2) {
         fprintf(stderr,"we must have NPROC_XI >= 2 and NPROC_ETA >= 2 to be able to use MPI, exiting...\n");
         MPI_Abort(FTI_COMM_WORLD,1);
         exit(1);
       }
#endif

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
#ifdef USE_ZERO_COPY
// allocation of pinned memory needs to be done first because often if it is done
// after other large standard memory allocations it fails by lack of space

// set the flag to allow access to pinned host memory by the device.
        print_CUDA_error_if_any(cudaSetDeviceFlags(cudaDeviceMapHost),1122);

// allocate buffer written by the device to pinned host memory
        print_CUDA_error_if_any(cudaHostAlloc((void**)&accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float), cudaHostAllocMapped),23);

// retrieve the device pointers corresponding to the mapped, pinned host buffers allocated by cudaHostAlloc()
        print_CUDA_error_if_any(cudaHostGetDevicePointer((void**)&d_accel_all_MPI_buffers,(void*)accel_all_MPI_buffers,0),27);
#else
// from the CUDA manual:
// cudaMallocHost() allocates host memory that is page-locked and accessible to the device.
// The driver tracks the virtual memory ranges allocated with this function and automatically
// accelerates calls to functions such as cudaMemcpy*(). Since the memory can be accessed directly
// by the device, it can be read or written with much higher bandwidth than pageable memory obtained
// with functions such as malloc(). Allocating excessive amounts of memory with cudaMallocHost() may
// degrade system performance, since it reduces the amount of memory available to the system for paging.
// As a result, this function is best used sparingly to allocate staging areas for data exchange between
// host and device.
        print_CUDA_error_if_any(cudaMallocHost((void**)&accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float)),1723);
        print_CUDA_error_if_any(cudaMalloc((void**)&d_accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float)),1727);
#endif
#endif

//
// device arrays
//

// global arrays
        float* d_displ;
        float* d_veloc;
        float* d_accel;

// local arrays
        int* d_ibool;

        float* d_rmass_inverse;

        float* d_xix;
        float* d_xiy;
        float* d_xiz;
        float* d_etax;
        float* d_etay;
        float* d_etaz;
        float* d_gammax;
        float* d_gammay;
        float* d_gammaz;

        float* d_kappav;
        float* d_muv;

//
// device memory allocation
//

// slices of mesher arrays
        print_CUDA_error_if_any(cudaMalloc((void**) &d_xix, NSPEC*NGLL3_ALIGN*sizeof(float)),1);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_xiy, NSPEC*NGLL3_ALIGN*sizeof(float)),2);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_xiz, NSPEC*NGLL3_ALIGN*sizeof(float)),3);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_etax, NSPEC*NGLL3_ALIGN*sizeof(float)),4);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_etay, NSPEC*NGLL3_ALIGN*sizeof(float)),5);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_etaz, NSPEC*NGLL3_ALIGN*sizeof(float)),6);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_gammax, NSPEC*NGLL3_ALIGN*sizeof(float)),7);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_gammay, NSPEC*NGLL3_ALIGN*sizeof(float)),8);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_gammaz, NSPEC*NGLL3_ALIGN*sizeof(float)),9);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_kappav, NSPEC*NGLL3_ALIGN*sizeof(float)),10);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_muv, NSPEC*NGLL3_ALIGN*sizeof(float)),11);

        print_CUDA_error_if_any(cudaMalloc((void**) &d_ibool, NSPEC*NGLL3_ALIGN*sizeof(int)),12);

// global arrays on the device
        print_CUDA_error_if_any(cudaMalloc((void**) &d_displ, NDIM*NGLOB*sizeof(float)),13);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_veloc, NDIM*NGLOB*sizeof(float)),14);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_accel, NDIM*NGLOB*sizeof(float)),15);
        print_CUDA_error_if_any(cudaMalloc((void**) &d_rmass_inverse, NGLOB*sizeof(float)),22);

#ifdef USE_TEXTURES
        bindTexturesDispl(d_displ);
        bindTexturesAccel(d_accel);
#endif

//
// host local arrays
//
        float* xix;
        float* xiy;
        float* xiz;
        float* etax;
        float* etay;
        float* etaz;
        float* gammax;
        float* gammay;
        float* gammaz;
        float* kappav;
        float* muv;

        int* ibool;

#ifdef DISPLAY_POSITION_SOURCE_RECEIVER
    double xstore[NSPEC][NGLLZ][NGLLY][NGLLX];
    double ystore[NSPEC][NGLLZ][NGLLY][NGLLX];
    double zstore[NSPEC][NGLLZ][NGLLY][NGLLX];
    int read_mesh_coordinate_arrays = 1;
#else
    double xstore[1][1][1][1];
    double ystore[1][1][1][1];
    double zstore[1][1][1][1];
    int read_mesh_coordinate_arrays = 0;
#endif

//
// host memory allocation
//
        xix = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        xiy = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        xiz = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        etax = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        etay = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        etaz = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        gammax = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        gammay = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        gammaz = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        kappav = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));
        muv = (float*)malloc(NSPEC*NGLL3_ALIGN*sizeof(float));

        ibool = (int*)malloc(NSPEC*NGLL3_ALIGN*sizeof(int));

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
////////////// DK DK implement a cascade to read the files serially otherwise there are crashes at CEA/CCRT in Paris
#ifdef USE_MPI
#ifdef USE_SERIAL_CASCADE_FOR_IOs
        if (myrank > 0) { MPI_Recv(&you_can_go_ahead_and_start_reading, 1, MPI_INT, myrank-1, 0, FTI_COMM_WORLD, &status); }
#endif
#endif
// read 2-D addressing for summation between slices: indirection arrays used to fill the MPI buffers
        read_arrays_buffers_solver(prname,iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta,&npoin2D_xi,&npoin2D_eta);

        if (npoin2D_xi != NGLOB2DMAX_XMIN_XMAX)
              { fprintf(stderr,"We have npoin2D_xi = %d\n",npoin2D_xi);
                fprintf(stderr,"and NGLOB2DMAX_XMIN_XMAX = %d\n",NGLOB2DMAX_XMIN_XMAX);
                fprintf(stderr,"We should have npoin2D_xi == NGLOB2DMAX_XMIN_XMAX, exiting...\n\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1); }
        if (npoin2D_eta != NGLOB2DMAX_YMIN_YMAX)
              { fprintf(stderr,"We have npoin2D_eta = %d\n",npoin2D_eta);
                fprintf(stderr,"and NGLOB2DMAX_YMIN_YMAX = %d\n",NGLOB2DMAX_YMIN_YMAX);
                fprintf(stderr,"We should have npoin2D_eta == NGLOB2DMAX_YMIN_YMAX, exiting...\n\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1); }

        threads_buf = BLOCK_SIZE_K_BUF;
        grid_buf_xi = (int)(ceilf((float)npoin2D_xi/(float)BLOCK_SIZE_K_BUF));
        grid_buf_eta = (int)(ceilf((float)npoin2D_eta/(float)BLOCK_SIZE_K_BUF));
        grid_buf_xi_eta = (int)(ceilf((float)NGLOB2DMAX_ALL/(float)BLOCK_SIZE_K_BUF));

////////////// DK DK implement a cascade to read the files serially otherwise there are crashes at CEA/CCRT in Paris
#ifdef USE_MPI
#ifdef USE_SERIAL_CASCADE_FOR_IOs
        you_can_go_ahead_and_start_reading = TRUE;
        if (myrank < nbprocs-1) { MPI_Send(&you_can_go_ahead_and_start_reading, 1, MPI_INT, myrank+1, 0, FTI_COMM_WORLD); }
#endif
#endif

        print_CUDA_error_if_any(cudaMalloc((void**)&d_iboolright_xi,npoin2D_xi*sizeof(int)),270);
        print_CUDA_error_if_any(cudaMalloc((void**)&d_iboolleft_xi,npoin2D_xi*sizeof(int)),280);
        print_CUDA_error_if_any(cudaMalloc((void**)&d_iboolright_eta,npoin2D_eta*sizeof(int)),290);
        print_CUDA_error_if_any(cudaMalloc((void**)&d_iboolleft_eta,npoin2D_eta*sizeof(int)),300);

        print_CUDA_error_if_any(cudaMemcpy(d_iboolright_xi,iboolright_xi,npoin2D_xi*sizeof(int),cudaMemcpyHostToDevice),91);
        print_CUDA_error_if_any(cudaMemcpy(d_iboolleft_xi,iboolleft_xi,npoin2D_xi*sizeof(int),cudaMemcpyHostToDevice),92);
        print_CUDA_error_if_any(cudaMemcpy(d_iboolright_eta,iboolright_eta,npoin2D_eta*sizeof(int),cudaMemcpyHostToDevice),93);
        print_CUDA_error_if_any(cudaMemcpy(d_iboolleft_eta,iboolleft_eta,npoin2D_eta*sizeof(int),cudaMemcpyHostToDevice),94);

        ibool1D_leftxi_lefteta = (int*)malloc(npoin2D_eta*sizeof(int));
        ibool1D_leftxi_righteta = (int*)malloc(npoin2D_eta*sizeof(int));
        ibool1D_rightxi_lefteta = (int*)malloc(npoin2D_eta*sizeof(int));
        ibool1D_rightxi_righteta = (int*)malloc(npoin2D_eta*sizeof(int));

        for (j=0;j<npoin2D_eta;j++) {
                ibool1D_leftxi_lefteta[j] = -1;
                ibool1D_leftxi_righteta[j] = -1;
                ibool1D_rightxi_lefteta[j] = -1;
                ibool1D_rightxi_righteta[j] = -1;
        }

        buffer_send_xi_to_left = (float*)malloc(NDIM*npoin2D_xi*sizeof(float));
        buffer_recv_xi_from_left = (float*)malloc(NDIM*npoin2D_xi*sizeof(float));
        buffer_send_xi_to_right = (float*)malloc(NDIM*npoin2D_xi*sizeof(float));
        buffer_recv_xi_from_right = (float*)malloc(NDIM*npoin2D_xi*sizeof(float));

        buffer_send_eta_to_up = (float*)malloc(NDIM*npoin2D_eta*sizeof(float));
        buffer_recv_eta_from_up = (float*)malloc(NDIM*npoin2D_eta*sizeof(float));
        buffer_send_eta_to_down = (float*)malloc(NDIM*npoin2D_eta*sizeof(float));
        buffer_recv_eta_from_down = (float*)malloc(NDIM*npoin2D_eta*sizeof(float));

// get indirection from faces to edges of mesh slices
        get_indirect_addressing_1D_buffers(iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta,npoin2D_xi,npoin2D_eta,ibool1D_leftxi_lefteta,ibool1D_leftxi_righteta,ibool1D_rightxi_lefteta,ibool1D_rightxi_righteta);

// read slice addressing (i.e., the MPI topology of the mesh slices)
        if (myrank == 0) {
//      printf("proc %d reads file addressing.txt\n",myrank);
        if((IIN=fopen("../DATABASES_FOR_SOLVER/addressing.txt","r"))==NULL) {
                fprintf(stderr,"Cannot open file addressing.txt, exiting...\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
        for (i=0;i<nbprocs;i++) {
                fscanf(IIN,"%d\n",&iproc_num);
                fscanf(IIN,"%d\n",&iproc_xi);
                fscanf(IIN,"%d\n",&iproc_eta);
                addressing[iproc_xi][iproc_eta] = iproc_num;
                iproc_xi_slice[iproc_num] = iproc_xi;
                iproc_eta_slice[iproc_num] = iproc_eta;
        }
        fclose(IIN);
        }

// broadcast the information read on the master to the nodes
        MPI_Bcast(addressing,NPROC_XI*NPROC_ETA,MPI_INT,0,FTI_COMM_WORLD);
        MPI_Bcast(iproc_xi_slice,NPROC_XI*NPROC_ETA,MPI_INT,0,FTI_COMM_WORLD);
        MPI_Bcast(iproc_eta_slice,NPROC_XI*NPROC_ETA,MPI_INT,0,FTI_COMM_WORLD);

        iproc_xi = iproc_xi_slice[myrank];
        iproc_eta = iproc_eta_slice[myrank];
#endif

////////////// DK DK implement a cascade to read the files serially otherwise there are crashes at CEA/CCRT in Paris
#ifdef USE_MPI
#ifdef USE_SERIAL_CASCADE_FOR_IOs
        if (myrank > 0) { MPI_Recv(&you_can_go_ahead_and_start_reading, 1, MPI_INT, myrank-1, 0, FTI_COMM_WORLD, &status); }
#endif
#endif

// read the mesh from an external file
//// DK DK 33333333333333 now in Fortran
//// DK DK 33333333333333 but still open and close the file just to check that it exists on the disk and exit if not
        if (myrank == 0) { printf("proc %d reads file database.dat\n",myrank);}
        strcpy(filename,prname);
        strcat(filename,"database.dat");
        if((IIN=fopen(filename,"r"))==NULL) {
                fprintf(stderr,"Cannot open file database.dat, exiting...\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
        fclose(IIN);

//// DK DK 33333333333333 now in Fortran
        read_arrays_solver_(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,kappav,muv,ibool,rmass_inverse,&myrank,xstore,ystore,zstore,&read_mesh_coordinate_arrays);

        for (ispec=0;ispec<NSPEC;ispec++) {
                for (k=0;k<NGLLZ;k++) {
                        for (j=0;j<NGLLY;j++) {
                                for (i=0;i<NGLLX;i++) {
// read real numbers here
//                                      fscanf(IIN, "%e\n", &xix[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &xiy[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &xiz[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &etax[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &etay[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &etaz[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &gammax[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &gammay[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &gammaz[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &kappav[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
//                                      fscanf(IIN, "%e\n", &muv[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
// read an integer here
//                                      fscanf(IIN, "%d\n", &ibool[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]);
// subtract one because indices start at zero in C but this array was created by a Fortran
// program and therefore starts at one in the file stored on the disk
                                        ibool[ispec*NGLL3_ALIGN+k*NGLL2+j*NGLLX+i]--;
                                }
                        }
                }
        }
//      for (i=0;i<NGLOB;i++) {
// at this point we read the unassembled mass matrix, not its inverse,
// because we first need to assemble it with MPI before being able to invert it
//              fscanf(IIN,"%e\n", &rmass_inverse[i]);
//      }
//// DK DK 33333333333333 now in Fortran
//      fclose(IIN);

// read pointers to the several colored parts of the mesh
//      printf("proc %d reads file number_of_elements_in_this_color.dat\n",myrank);
        strcpy(filename,prname);
        strcat(filename,"number_of_elements_in_this_color.dat");
        if((IIN=fopen(filename,"r"))==NULL) {
                fprintf(stderr,"Cannot open file number_of_elements_in_this_color.dat, exiting...\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
// nb_color[0]: number of colors for the outer elements
// nb_color[1]: number of colors for the inner elements
        fscanf(IIN,"%d\n", &nb_color[0]); // outer elements
        fscanf(IIN,"%d\n", &nb_color[1]); // inner elements

        nb_color_tot = nb_color[0] + nb_color[1];

// array of indices that indicate the beginning of each color
        number_of_elements_in_this_color = (int*) malloc(nb_color_tot*sizeof(int));

        for(i = 0; i < nb_color_tot; i++) {
                fscanf(IIN,"%d\n",&number_of_elements_in_this_color[i]);
        }
        fclose(IIN);

////////////// DK DK implement a cascade to read the files serially otherwise there are crashes at CEA/CCRT in Paris
#ifdef USE_MPI
#ifdef USE_SERIAL_CASCADE_FOR_IOs
        you_can_go_ahead_and_start_reading = TRUE;
        if (myrank < nbprocs-1) { MPI_Send(&you_can_go_ahead_and_start_reading, 1, MPI_INT, myrank+1, 0, FTI_COMM_WORLD); }
#endif
#endif

// read the derivative matrices from an external file
        if (myrank == 0) {
//       printf("proc %d reads file matrices.dat\n",myrank);
        if((IIN=fopen("../DATABASES_FOR_SOLVER/matrices.dat","r"))==NULL) {
                fprintf(stderr,"Cannot open file matrices.dat, exiting...\n");
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }

        for (j=0;j<NGLLY;j++) {
                for (i=0;i<NGLLX;i++) {
                        fscanf(IIN, "%e\n", &hprime_xx[j*NGLLX+i]);
                        fscanf(IIN, "%e\n", &hprimewgll_xx[j*NGLLX+i]);
                        fscanf(IIN, "%e\n", &wgllwgll_yz[j*NGLLX+i]);
                        fscanf(IIN, "%e\n", &wgllwgll_xz[j*NGLLX+i]);
                        fscanf(IIN, "%e\n", &wgllwgll_xy[j*NGLLX+i]);
                }
        }
        fclose(IIN);
        }

// broadcast the information read on the master to the nodes
#ifdef USE_MPI
        MPI_Bcast(hprime_xx,NGLL2,MPI_FLOAT,0,FTI_COMM_WORLD);
        MPI_Bcast(hprimewgll_xx,NGLL2,MPI_FLOAT,0,FTI_COMM_WORLD);
        MPI_Bcast(wgllwgll_yz,NGLL2,MPI_FLOAT,0,FTI_COMM_WORLD);
        MPI_Bcast(wgllwgll_xz,NGLL2,MPI_FLOAT,0,FTI_COMM_WORLD);
        MPI_Bcast(wgllwgll_xy,NGLL2,MPI_FLOAT,0,FTI_COMM_WORLD);
#endif

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
// assemble the mass matrix once and for all before the time loop using non-blocking MPI
//    printf("proc %d assembles its mass matrix with blocking MPI\n",myrank);
      assemble_MPI_blocking_scalar(iboolleft_xi, iboolright_xi, iboolleft_eta, iboolright_eta, iproc_xi, iproc_eta, addressing, npoin2D_xi, npoin2D_eta, buffer_send_xi_to_right, buffer_recv_xi_from_left, buffer_send_eta_to_up, buffer_recv_eta_from_down,rmass_inverse);
#endif

// invert the diagonal mass matrix once and for all after assembling it with MPI above
      for (i = 0; i < NGLOB; i++) { rmass_inverse[i] = 1.f/rmass_inverse[i]; }

// copy these arrays to the device
// we could/should store the mass matrix "rmass_inverse" in texture memory in the future
//      printf("proc %d exports its data to the GPU\n",myrank);
        print_CUDA_error_if_any(cudaMemcpy(d_rmass_inverse,rmass_inverse,NGLOB*sizeof(float),cudaMemcpyHostToDevice),31);
        print_CUDA_error_if_any(cudaMemcpy(d_xix,xix,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),32);
        print_CUDA_error_if_any(cudaMemcpy(d_xiy,xiy,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),33);
        print_CUDA_error_if_any(cudaMemcpy(d_xiz,xiz,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),34);
        print_CUDA_error_if_any(cudaMemcpy(d_etax,etax,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),35);
        print_CUDA_error_if_any(cudaMemcpy(d_etay,etay,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),36);
        print_CUDA_error_if_any(cudaMemcpy(d_etaz,etaz,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),37);
        print_CUDA_error_if_any(cudaMemcpy(d_gammax,gammax,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),38);
        print_CUDA_error_if_any(cudaMemcpy(d_gammay,gammay,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),39);
        print_CUDA_error_if_any(cudaMemcpy(d_gammaz,gammaz,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),40);
        print_CUDA_error_if_any(cudaMemcpy(d_kappav,kappav,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),41);
        print_CUDA_error_if_any(cudaMemcpy(d_muv,muv,NSPEC*NGLL3_ALIGN*sizeof(float),cudaMemcpyHostToDevice),42);
        print_CUDA_error_if_any(cudaMemcpy(d_ibool,ibool,NSPEC*NGLL3_ALIGN*sizeof(int),cudaMemcpyHostToDevice),43);

// copy these arrays to the device in constant memory
        setConst_hprime_xx(hprime_xx);
        setConst_hprimewgll_xx(hprimewgll_xx);
        setConst_wgllwgll_yz(wgllwgll_yz);
        setConst_wgllwgll_xz(wgllwgll_xz);
        setConst_wgllwgll_xy(wgllwgll_xy);

// then free these arrays on the host
        free(xix);
        free(xiy);
        free(xiz);
        free(etax);
        free(etay);
        free(etaz);
        free(gammax);
        free(gammay);
        free(gammaz);
        free(kappav);
        free(muv);

// clear initial vectors before starting the time loop
        print_CUDA_error_if_any(cudaMemset(d_displ,0,NDIM*NGLOB*sizeof(float)),49);
        print_CUDA_error_if_any(cudaMemset(d_veloc,0,NDIM*NGLOB*sizeof(float)),52);
        print_CUDA_error_if_any(cudaMemset(d_accel,0,NDIM*NGLOB*sizeof(float)),55);

/////////////        k2_tot = 0.f;

// definition of grid sizes for kernels 1 3 4
    int grid_134_x, grid_134_y;
    int nb_blocks = (int)ceil((float)NGLOB/(float)BLOCK_SIZE_K1_K3_K4);
    if (nb_blocks > MAX_ONE_DIMENSION_CUDA_GRID) {
      // 2D grid
      grid_134_x = (int)floor(sqrt((float)nb_blocks));
      grid_134_y = (int)ceil((float)nb_blocks/(float)grid_134_x);
    } else {
      // 1D grid
      grid_134_x = nb_blocks;
      grid_134_y = 1;
    }


// make sure everybody is ready and is done reading the input files before starting the time loop
// we do not care if this barrier is expensive because timing of the run starts in the time loop only
// and therefore this part is ignored. And it is always better for everybody to start the timed
// section of the code at the same time
        cudaThreadSynchronize();
#ifdef USE_MPI
        MPI_Barrier(FTI_COMM_WORLD);
#endif

// compute the total number of processors used
   NPROCTOT = NPROC_XI*NPROC_ETA;

// print the position of the source
#ifdef DISPLAY_POSITION_SOURCE_RECEIVER
#ifdef USE_MPI
                if (myrank == RANK_SOURCE)
#else
                if (myrank == 0)
#endif
                {
                        printf("\nSource X = %.15lf\n",xstore[NSPEC_SOURCE-1][1][1][1]);
                        printf("Source Y = %.15lf\n",ystore[NSPEC_SOURCE-1][1][1][1]);
                        printf("Source Z = %.15lf\n\n",zstore[NSPEC_SOURCE-1][1][1][1]);
                }
#ifdef USE_MPI
        MPI_Barrier(FTI_COMM_WORLD);
#endif

// print the position of the station
#ifdef USE_MPI
                if (myrank == RANK_STATION)
#else
                if (myrank == 0)
#endif
                {
                        printf("Station X = %.15lf\n",xstore[NSPEC_STATION-1][1][1][1]);
                        printf("Station Y = %.15lf\n",ystore[NSPEC_STATION-1][1][1][1]);
                        printf("Station Z = %.15lf\n\n",zstore[NSPEC_STATION-1][1][1][1]);
                }
#ifdef USE_MPI
        MPI_Barrier(FTI_COMM_WORLD);
#endif
#endif

// start the timer
        timeloop_begin = clock();
        t2 = MPI_Wtime();
        t8 = 999999999;

//******************************** LOOP *************************
// beginning of the serial time loop
        for (it = 1; it <= NSTEP/DIV; it++) {

            // compute the maximum of the norm of the displacement vector from time to time and display it
            // in order to monitor the simulation
            // this can remain serial because it is done only every NTSTEP_BETWEEN_OUTPUT_INFO time steps
            if (((it % NTSTEP_BETWEEN_OUTPUT_INFO) == 0 || it == 5 || it == NSTEP) && 0) {

                print_CUDA_error_if_any(cudaMemcpy(displ,d_displ,NDIM*NGLOB*sizeof(float),cudaMemcpyDeviceToHost),9961);

                // compute the maximum of the norm of the displacement vector in this mesh slice
                Usolidnorm_local = -1.f;

                for (iglob = 0; iglob < NGLOB; iglob++) {
                    current_value = sqrtf(displ[X][iglob]*displ[X][iglob] + displ[Y][iglob]*displ[Y][iglob] + displ[Z][iglob]*displ[Z][iglob]);
                    if(current_value > Usolidnorm_local) { Usolidnorm_local = current_value; }
                }

                // broadcast the information read on the master to the nodes
                MPI_Reduce(&Usolidnorm_local,&Usolidnorm_global,1,MPI_FLOAT,MPI_MAX,0,FTI_COMM_WORLD);
                Usolidnorm_global = Usolidnorm_local;

                if (myrank == 0) {

                    // get the timer
                    timeloop_total = ((clock()-timeloop_begin)/(float)CLOCKS_PER_SEC);

                    printf("\nTime step # %d out of %d\n",it,NSTEP);
                    printf("Max norm displacement vector U in the solid (m) = %.8g\n",Usolidnorm_global);
                    printf("Total elapsed time so far for NPROCTOT = %03d : %f\n",NPROCTOT,timeloop_total);
                    if (it >= 100) { printf("Average elapsed time per time step for NPROCTOT = %03d : %f\n",NPROCTOT,timeloop_total/(float)(it-1)); }

                    // write a time stamp file
                    sprintf(prname,"timestamp_%07d_NPROCTOT_%03d_NPROCXI_%02d_NPROCETA_%02d.txt",it,NPROCTOT,NPROC_XI,NPROC_ETA);
                    if((IIN = fopen(prname,"w")) == NULL) {
                        fprintf(stderr,"Cannot create time stamp file, exiting...\n");
                        MPI_Abort(FTI_COMM_WORLD,1);
                        exit(1);
                    }
                    fprintf(IIN,"\nTime step # %d out of %d\n",it,NSTEP);
                    fprintf(IIN,"Max norm displacement vector U in the solid (m) = %.8g\n",Usolidnorm_global);
                    fprintf(IIN,"Total elapsed time so far for NPROCTOT = %03d : %f\n",NPROCTOT,timeloop_total);
                    if (it >= 100) { fprintf(IIN,"Average elapsed time per time step for NPROCTOT = %03d : %f\n",NPROCTOT,timeloop_total/(float)(it-1)); }
                    fprintf(IIN,"\n");
                    fclose(IIN);

                    // check stability of the code, exit if unstable
                    if(Usolidnorm_global > STABILITY_THRESHOLD || Usolidnorm_global < 0) {
                        fprintf(stderr,"code became unstable and blew up\n");
                        MPI_Abort(FTI_COMM_WORLD,1);
                        exit(1);
                    }

                }

            }

/* *********************************************************************************** */
// KERNEL 1
                UpdateDisplVeloc (d_displ, d_veloc, d_accel, grid_134_x, grid_134_y);

// synchronize
// DOM: The calls to cudaMemset below are blocking, so there is no need to sync here
#ifdef ADD_SYNC_BARRIERS
                cudaThreadSynchronize();
#endif

// set acceleration vector to 0
// DOM: This is definitely necessary because of the indirection access to these vectors and the accumulation to them in kernel_2.
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
                print_CUDA_error_if_any(cudaMemset(d_accel,0,NDIM*NGLOB*sizeof(float)),58);
#else
                cudaMemset(d_accel,0,NDIM*NGLOB*sizeof(float));
#endif

                /*******************************************************************************/
                // CHECKPOINT !!!
                if (it % 3000 == 0) {
                    t5 = MPI_Wtime();
                    print_CUDA_error_if_any(cudaMemcpy(displ,d_displ,NDIM*NGLOB*sizeof(float),cudaMemcpyDeviceToHost),81);
                    print_CUDA_error_if_any(cudaMemcpy(veloc,d_veloc,NDIM*NGLOB*sizeof(float),cudaMemcpyDeviceToHost),82);
                    int jt;

                    sprintf(ckptfile,"%s/%s",FTI_Cdir, FTI_File);
                    fd = fopen(ckptfile, "w");
                    fwrite(&it,sizeof(int),1,fd);
                    for(jt = 0; jt < NDIM; jt++) {
                        fwrite(displ[jt],sizeof(float),NGLOB,fd);
                        fwrite(veloc[jt],sizeof(float),NGLOB,fd);
                    }
                    fclose(fd);
                    FTI_Checkpointed();
                    t6 = t6 + (MPI_Wtime() - t5);
                }

                // Recovering in case of failure
                if (FTI_Fail) {
                    int jt;
                    sprintf(ckptfile,"%s/%s",FTI_Cdir, FTI_File);
                    fd = fopen(ckptfile, "r");
                    fread(&it,sizeof(int),1,fd);
                    for(jt = 0; jt < NDIM; jt++) {
                        fread(displ[jt],sizeof(float),NGLOB,fd);
                        fread(veloc[jt],sizeof(float),NGLOB,fd);
                    }
                    fclose(fd);
                    print_CUDA_error_if_any(cudaMemcpy(d_displ,displ,NDIM*NGLOB*sizeof(float),cudaMemcpyHostToDevice),83);
                    print_CUDA_error_if_any(cudaMemcpy(d_veloc,veloc,NDIM*NGLOB*sizeof(float),cudaMemcpyHostToDevice),84);
                    FTI_Restarted();
                }

                /*************************************************************************************/
// KERNEL 2
#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
                assembling_phase = 0;
#endif
                int color_offset = 0;

// loop on color sets
                for(int icolor = 0; icolor < nb_color_tot; icolor++) {

// we first loop on the color sets of the outer elements.
// then we prepare the MPI buffers and we start the non-blocking communications.
// finally, we loop on the color sets of the inner elements while the MPI messages
// are traveling across the network. This way the MPI communications are overlapped
// by calculations on the GPUs.
#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
// nb_color[0]: number of colors for the outer elements
// nb_color[1]: number of colors for the inner elements
                        if (icolor < nb_color[0]) {
                                we_are_done_computing_outer_elements_first = FALSE;
                        } else {
                                we_are_done_computing_outer_elements_first = TRUE;
                        }
#endif

// number of elements in that color
            nb_elem_color = number_of_elements_in_this_color[icolor];

// if necessary, split the color set if its size is bigger than the size of the grid
            int grid_2_x, grid_2_y;
            if (nb_elem_color > MAX_ONE_DIMENSION_CUDA_GRID) {
              // 2D GRID
              grid_2_x = (int)floor(sqrt((float)nb_elem_color));
              grid_2_y = (int)ceil((float)nb_elem_color/(float)grid_2_x);
            } else {
              // 1D GRID
              grid_2_x = nb_elem_color;
              grid_2_y = 1;
            }
/////////////            k2_beg = clock();
            Kernel_2(nb_elem_color, grid_2_x, grid_2_y, d_ibool+color_offset, d_displ, d_accel, d_xix+color_offset, d_xiy+color_offset, d_xiz+color_offset, d_etax+color_offset, d_etay+color_offset, d_etaz+color_offset, d_gammax+color_offset, d_gammay+color_offset, d_gammaz+color_offset, d_kappav+color_offset, d_muv+color_offset);

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
// this part of the code runs on the host in parallel to the execution of kernel 2 on the GPU
// DOM: Indeed it does, and as Kernel_2() is very heavy, I doubt that anything can be gained here
//      by optimising compiler flags for the CPU code.
// get data from the device to fill the MPI communication buffers on the host
// to be able to assemble the contributions of neighboring mesh slices
///////            if (!we_are_done_computing_outer_elements_first && icolor==(nb_color[0]-1)) { // after the last loop on the outer elements
/////////// DK DK the double test was redundant
            if (icolor == (nb_color[0]-1)) { // right after the last loop on the outer elements

// launch the kernel for the buffers along the "xi" and "eta" directions
// wait for kernel 2 to finish
// DOM: I believe this sync is not necessary at all, as there is always an implicit sync between kernel
//      calls on dependent data
#ifdef ADD_SYNC_BARRIERS
                cudaThreadSynchronize();
#endif
                get_MPI_buffers(grid_buf_xi_eta,threads_buf,npoin2D_xi,npoin2D_eta,d_accel_all_MPI_buffers, d_accel, d_iboolright_xi, d_iboolleft_xi,d_iboolright_eta, d_iboolleft_eta);

// transfer the buffers from device to host
#ifndef USE_ZERO_COPY
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
                print_CUDA_error_if_any(cudaMemcpy(accel_all_MPI_buffers,d_accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float),cudaMemcpyDeviceToHost),63);
#else
                cudaMemcpy(accel_all_MPI_buffers,d_accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float),cudaMemcpyDeviceToHost);
#endif
#endif
            }

// assemble the contributions with non-blocking MPI
            if (we_are_done_computing_outer_elements_first) {
// this flag indicates that the assembling phase is not yet finished, because there are 4 phases therefore 4 is the maximum
                if (assembling_phase < 4) {
                    assemble_MPI(iproc_xi, iproc_eta, addressing, npoin2D_xi, npoin2D_eta, &assembling_phase,
                                 &request_send_xi_to_left, &request_recv_xi_from_left,
                                 &request_send_xi_to_right, &request_recv_xi_from_right,
                                 &request_send_eta_to_up, &request_recv_eta_from_up,
                                 &request_send_eta_to_down, &request_recv_eta_from_down,
                                 buffer_send_xi_to_left, buffer_recv_xi_from_left,
                                 buffer_send_xi_to_right, buffer_recv_xi_from_right,
                                 buffer_send_eta_to_up, buffer_recv_eta_from_up,
                                 buffer_send_eta_to_down, buffer_recv_eta_from_down,
                                 ibool1D_leftxi_lefteta, ibool1D_leftxi_righteta, ibool1D_rightxi_lefteta, ibool1D_rightxi_righteta, accel_all_MPI_buffers);

// if last color of the inner elements i.e. last color of the whole mesh slice, we wait for the last communications to arrive and then we start to assemble
/////////////                    wait_flag = 0;
                    if (icolor == (nb_color_tot-1)) {
                      while (assembling_phase < 4) {
/////////////                                      if (!wait_flag) {
// timer for the blocking wait for the last communications
/////////////                                          waiting_beg = clock();
/////////////                                          wait_flag++;
/////////////                                      }
                        assemble_MPI(iproc_xi, iproc_eta, addressing, npoin2D_xi, npoin2D_eta, &assembling_phase,
                                     &request_send_xi_to_left, &request_recv_xi_from_left,
                                     &request_send_xi_to_right, &request_recv_xi_from_right,
                                     &request_send_eta_to_up, &request_recv_eta_from_up,
                                     &request_send_eta_to_down, &request_recv_eta_from_down,
                                     buffer_send_xi_to_left, buffer_recv_xi_from_left,
                                     buffer_send_xi_to_right, buffer_recv_xi_from_right,
                                     buffer_send_eta_to_up, buffer_recv_eta_from_up,
                                     buffer_send_eta_to_down, buffer_recv_eta_from_down,
                                     ibool1D_leftxi_lefteta, ibool1D_leftxi_righteta, ibool1D_rightxi_lefteta, ibool1D_rightxi_righteta, accel_all_MPI_buffers);
                      }
                    }
/////////////                                if (wait_flag) {
/////////////                                    waiting_end = clock();
/////////////                                    waiting_tot += ((waiting_end-waiting_beg)/(float)CLOCKS_PER_SEC);
/////////////                                }
                }
                if (assembling_phase == 4) {
// copy the assembled communication buffers to the card
#ifdef ADD_SYNC_BARRIERS
                    cudaThreadSynchronize();
#endif

#ifndef USE_ZERO_COPY
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
                    print_CUDA_error_if_any(cudaMemcpy(d_accel_all_MPI_buffers,accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float),cudaMemcpyHostToDevice),67);
#else
                    cudaMemcpy(d_accel_all_MPI_buffers,accel_all_MPI_buffers,4*NGLOB2DMAX_ALL*NDIM*sizeof(float),cudaMemcpyHostToDevice);
#endif
#endif

// use a CUDA kernel to update the d_accel vector
// we cannot launch the "xi" and "eta" directions in the same kernel
// because the "xi" and "eta" faces write to the same point in the corners of each mesh slice (with are common to both faces)

// launch the kernel for the buffers along the "xi" direction
#ifdef ADD_SYNC_BARRIERS
                    cudaThreadSynchronize();
#endif
                    update_accel_MPI_xi(grid_buf_xi,threads_buf,npoin2D_xi,d_accel_all_MPI_buffers, d_accel, d_iboolright_xi, d_iboolleft_xi);

// launch the kernel for the buffers along the "eta" direction
#ifdef ADD_SYNC_BARRIERS
                    cudaThreadSynchronize();
#endif
                    update_accel_MPI_eta(grid_buf_eta,threads_buf,npoin2D_eta,d_accel_all_MPI_buffers, d_accel, d_iboolright_eta, d_iboolleft_eta);

// this flag indicates that the assembling phase is finished, because there are 4 phases therefore 4 is the maximum
                    assembling_phase = 5;
                                        }
                                }
#endif  // of ACTUALLY_ASSEMBLE_MPI_SLICES

// synchronize
// DOM: again: sync is not necessary here, I am pretty sure of this. Unless you want k2_end to contain
//      something meaningful.
// DAVID:   Ok, added an ifdef, just let the last one for timing the timeloop without taking into account
//          the code to print basin information on stdout and the code to testing stability
#ifdef ADD_SYNC_BARRIERS
                    cudaThreadSynchronize();
#endif

/////////////                    k2_end = clock();
/////////////                    k2_tot += ((k2_end-k2_beg)/(float)CLOCKS_PER_SEC);
                    color_offset += nb_elem_color*NGLL3_ALIGN;
                } // end of loop on colors

/*************************************************************************************/
// KERNEL 3
//               Kernel_3(d_accel,d_rmass_inverse);
#ifdef ADD_SYNC_BARRIERS
//               cudaThreadSynchronize();
#endif

// add the earthquake source at a given grid point
// this is negligible in terms of computation time and is intrinsically serial because it is done by only
// one grid point out of several millions typically
#ifdef USE_MPI
                if (myrank == RANK_SOURCE)
#else
                if (myrank == 0)
#endif
                {
                        iglob = ibool[(NSPEC_SOURCE-1)*NGLL3_ALIGN+NGLL2+NGLLX+1];  // this is ibool[NSPEC_SOURCE-1][1][1][1]
                        time = (it-1)*deltat;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
// get the vertical component
                        print_CUDA_error_if_any(cudaMemcpy(&source_accel,d_accel+iglob + Z*NGLOB,sizeof(float),cudaMemcpyDeviceToHost),69);
#else
// get the vertical component
                        cudaMemcpy(&source_accel,d_accel+iglob + Z*NGLOB,sizeof(float),cudaMemcpyDeviceToHost);
#endif
// we divide the amplitude of the source by rmass_inverse[iglob_source] here because
// we have merged the calculation of acceleration and velocity below in a single kernel
// and therefore the value of accel[] at that point will be
// multiplied by rmass_inverse[i] in that merged kernel
                        source_accel += 1.e4f * (1.f - 2.f*a*(time-t0)*(time-t0)) * expf(-a*(time-t0)*(time-t0)) / (rho * rmass_inverse[iglob]);
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
// add the source to the vertical component
                        print_CUDA_error_if_any(cudaMemcpy(d_accel+iglob + Z*NGLOB,&source_accel,sizeof(float),cudaMemcpyHostToDevice),70);
#else
// add the source to the vertical component
                        cudaMemcpy(d_accel+iglob + Z*NGLOB,&source_accel,sizeof(float),cudaMemcpyHostToDevice);
#endif
                }

/*************************************************************************************/
// KERNEL 4
                Kernel_3_merged(d_veloc, d_accel,d_rmass_inverse, grid_134_x, grid_134_y);
//              not_used_any_more_Kernel_4(d_veloc,d_accel);


#ifdef WRITE_SEISMOGRAM
// record a seismogram (time evolution of the seismic wave field at a given
// seismic station in a given mesh slice) to check that the simulation went well
#ifdef USE_MPI
                if (myrank == RANK_STATION)
#else
                if (myrank == 0)
#endif
                {
                        iglob = ibool[(NSPEC_STATION-1)*NGLL3_ALIGN+NGLL2+NGLLX+1];  // this is ibool[NSPEC_STATION-1][1][1][1]
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
// get the vertical component
                        print_CUDA_error_if_any(cudaMemcpy(&seismogram_x[it-1],d_displ+iglob + X*NGLOB,sizeof(float),cudaMemcpyDeviceToHost),71);
                        print_CUDA_error_if_any(cudaMemcpy(&seismogram_y[it-1],d_displ+iglob + Y*NGLOB,sizeof(float),cudaMemcpyDeviceToHost),71);
                        print_CUDA_error_if_any(cudaMemcpy(&seismogram_z[it-1],d_displ+iglob + Z*NGLOB,sizeof(float),cudaMemcpyDeviceToHost),71);
#else
// get the vertical component
                        cudaMemcpy(&seismogram_x[it-1],d_displ+iglob + X*NGLOB,sizeof(float),cudaMemcpyDeviceToHost);
                        cudaMemcpy(&seismogram_y[it-1],d_displ+iglob + Y*NGLOB,sizeof(float),cudaMemcpyDeviceToHost);
                        cudaMemcpy(&seismogram_z[it-1],d_displ+iglob + Z*NGLOB,sizeof(float),cudaMemcpyDeviceToHost);
#endif
                }
#endif
                if (it == 1) {
                    t7 = MPI_Wtime();
                } else {
                    t7 = MPI_Wtime() - t7;
                }

                if (t7 < t8) {
                    t8 = t7;
                }

                //
// end of the serial time loop
//
        }
        t3 = MPI_Wtime();

#ifdef WRITE_SEISMOGRAM
// synchronize everything at the end of the time loop, just in case
        cudaThreadSynchronize();
#ifdef USE_MPI
        if (myrank == RANK_STATION)
#else
        if (myrank == 0)
#endif
        {       // save the seismogram file at the end of the run
                if((IIN=fopen("seismogram_CUDA.txt","w"))==NULL) {
                        fprintf(stderr,"Cannot open file seismogram_CUDA.txt, exiting...\n");
#ifdef USE_MPI
                        MPI_Abort(FTI_COMM_WORLD,1);
#endif
                        exit(1);
                }
                for (it=0;it<NSTEP;it++)
                {   fprintf(IIN,"%e  %e  %e  %e\n",it*deltat,seismogram_x[it],seismogram_y[it],seismogram_z[it]);
                }
                fclose(IIN);
        }
#endif

        free(number_of_elements_in_this_color);

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
        free(buffer_send_xi_to_left);
        free(buffer_recv_xi_from_left);
        free(buffer_send_xi_to_right);
        free(buffer_recv_xi_from_right);
        free(buffer_send_eta_to_up);
        free(buffer_recv_eta_from_up);
        free(buffer_send_eta_to_down);
        free(buffer_recv_eta_from_down);

        free(ibool1D_leftxi_lefteta);
        free(ibool1D_leftxi_righteta);
        free(ibool1D_rightxi_lefteta);
        free(ibool1D_rightxi_righteta);
#endif

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES
/////////////        printf("proc %d: of which %f were spent waiting for communications\n",myrank,waiting_tot);
#endif

// close MPI
#ifdef USE_MPI
        MPI_Barrier(FTI_COMM_WORLD);
        FTI_Finalize();

        t4 = MPI_Wtime();
        if (myrank == 0) {
            printf("\n ==> Time elapsed initializing the application : %f\n",t2-t1);
            printf(" ==> Time elapsed in the main loop : %f\n",t3-t2);
            printf(" ==> Time elapsed in checkpointing : %f\n",t6);
            printf(" ==> Time elapsed in computation in main loop : %f\n",(t3-t2)-t6);
            printf(" ==> Average time elapsed per time step : %f\n",(t3-t2)/(NSTEP/DIV));
            printf(" ==> Fastest loop : %f\n",t8);
            printf(" ==> Time elapsed finalizing the application : %f\n",t4-t3);
            printf(" ==> Total time elapsed : %f\n",t4-t1);
            printf("\n === End of the program === \n");
        }

        MPI_Finalize();
#endif
        return 0;
  }

//
// ***************************************************************************
//    below are the subroutines we need to assemble with MPI
// ***************************************************************************
//

#ifdef ACTUALLY_ASSEMBLE_MPI_SLICES

// read 2-D addressing for summation between slices: indirection arrays used to fill the MPI buffers
void read_arrays_buffers_solver(char* prname, int* iboolleft_xi, int* iboolright_xi, int* iboolleft_eta, int* iboolright_eta,int* p_npoin2D_xi,int* p_npoin2D_eta)
{
        char filename[300];
        int i,value_read;
        FILE* IIN;

        (*p_npoin2D_xi) = 0;
        (*p_npoin2D_eta) = 0;

// read 2-D addressing for summation between slices along xi with MPI
        sprintf(filename,"%siboolleft_xi.txt",prname);
        if((IIN=fopen(filename,"r"))==NULL) {
                fprintf(stderr,"Cannot open file %siboolleft_xi.txt, exiting...\n",prname);
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
        i = 0;
        do {
                fscanf(IIN,"%d\n",&value_read);
                if(value_read > 0) {
// subtract 1 because these arrays have been created by a Fortran code, in which indices start at 1 rather than 0
                  iboolleft_xi[i] = value_read - 1;
                  i++;
                  (*p_npoin2D_xi)++;
                }
        } while(value_read > 0);
        fclose(IIN);

        sprintf(filename,"%siboolright_xi.txt",prname);
        if((IIN=fopen(filename,"r"))==NULL) {
                fprintf(stderr,"Cannot open file %siboolright_xi.txt, exiting...\n",prname);
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
        i = 0;
        do {
                fscanf(IIN,"%d\n",&value_read);
                if(value_read > 0) {
// subtract 1 because these arrays have been created by a Fortran code, in which indices start at 1 rather than 0
                  iboolright_xi[i] = value_read - 1;
                  i++;
                }
        } while(value_read > 0);
        fclose(IIN);

// read 2-D addressing for summation between slices along eta with MPI
        sprintf(filename,"%siboolleft_eta.txt",prname);
        if((IIN=fopen(filename,"r"))==NULL) {
                fprintf(stderr,"Cannot open file %siboolleft_eta.txt, exiting...\n",prname);
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
        i = 0;
        do {
                fscanf(IIN,"%d\n",&value_read);
                if(value_read > 0) {
// subtract 1 because these arrays have been created by a Fortran code, in which indices start at 1 rather than 0
                  iboolleft_eta[i] = value_read - 1;
                  i++;
                  (*p_npoin2D_eta)++;
                }
        } while(value_read > 0);
        fclose(IIN);

        sprintf(filename,"%siboolright_eta.txt",prname);
        if((IIN=fopen(filename,"r"))==NULL) {
                fprintf(stderr,"Cannot open file %siboolright_eta.txt, exiting...\n",prname);
#ifdef USE_MPI
                MPI_Abort(FTI_COMM_WORLD,1);
#endif
                exit(1);
        }
        i = 0;
        do {
                fscanf(IIN,"%d\n",&value_read);
                if(value_read > 0) {
// subtract 1 because these arrays have been created by a Fortran code, in which indices start at 1 rather than 0
                  iboolright_eta[i] = value_read - 1;
                  i++;
                }
        } while(value_read > 0);
        fclose(IIN);

        return;
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// get indirection from faces to edges of mesh slices:
// get indirect addressing for the 1D MPI buffers for the vertical edges of mesh slice corners
// i.e. read the point numbers that compose them, which have been created by the mesher
void get_indirect_addressing_1D_buffers(int* iboolleft_xi, int* iboolright_xi, int* iboolleft_eta, int* iboolright_eta,int npoin2D_xi,int npoin2D_eta, int* ibool1D_leftxi_lefteta, int* ibool1D_leftxi_righteta, int* ibool1D_rightxi_lefteta, int* ibool1D_rightxi_righteta) {

        int i,j;

        for(i=0;i<npoin2D_xi;i++) {
                for(j=0;j<npoin2D_eta;j++) {
                        if (iboolleft_xi[i] == iboolleft_eta[j]) {
                                ibool1D_leftxi_lefteta[j] = i;
                        }
                        if (iboolleft_xi[i] == iboolright_eta[j]) {
                                ibool1D_leftxi_righteta[j] = i;
                        }
                        if (iboolright_xi[i] == iboolleft_eta[j]) {
                                ibool1D_rightxi_lefteta[j] = i;
                        }
                        if (iboolright_xi[i] == iboolright_eta[j]) {
                                ibool1D_rightxi_righteta[j] = i;
                        }
                }
        }
        return;
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

void assemble_MPI(int iproc_xi, int iproc_eta, int addressing[NPROC_XI][NPROC_ETA], int npoin2D_xi, int npoin2D_eta, int* phase,
    MPI_Request* request_send_xi_to_left, MPI_Request* request_recv_xi_from_left,
    MPI_Request* request_send_xi_to_right, MPI_Request* request_recv_xi_from_right,
    MPI_Request* request_send_eta_to_up, MPI_Request* request_recv_eta_from_up,
    MPI_Request* request_send_eta_to_down, MPI_Request* request_recv_eta_from_down,
    float* buffer_send_xi_to_left, float* buffer_recv_xi_from_left,
    float* buffer_send_xi_to_right, float* buffer_recv_xi_from_right,
    float* buffer_send_eta_to_up, float* buffer_recv_eta_from_up,
    float* buffer_send_eta_to_down, float* buffer_recv_eta_from_down,
    int* ibool1D_leftxi_lefteta, int* ibool1D_leftxi_righteta, int* ibool1D_rightxi_lefteta, int* ibool1D_rightxi_righteta, float* accel_all_MPI_buffers)
{

// IMPORTANT: after assembling, we first need to overwrite the values of "accel"
// by the values of the "xi" buffers, and then by the values of the "eta" buffers.
// the order ("xi" first, then "eta") does matter because the grid points that
// belong to the common vertical edge are up to date in the "eta" buffers only
// and therefore the corresponding values in the "xi" buffers, which are not correct,
// must be discarded

        int sender,receiver,i;
        MPI_Status status;

        int flag_send_xi_to_left, flag_recv_xi_from_left,
            flag_send_xi_to_right, flag_recv_xi_from_right,
            flag_send_eta_to_up, flag_recv_eta_from_up,
            flag_send_eta_to_down, flag_recv_eta_from_down;

        if (*phase == 0) {
// phase 0: send right_xi to left_xi
//////////// useless because erased when received                        memset(buffer_recv_xi_from_left,0,NDIM*npoin2D_xi*sizeof(float));
// all the slices copy their right face into the buffer
                        memcpy(&buffer_send_xi_to_right[X*npoin2D_xi],&accel_all_MPI_buffers[right_xi*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL],npoin2D_xi*sizeof(float));
                        memcpy(&buffer_send_xi_to_right[Y*npoin2D_xi],&accel_all_MPI_buffers[right_xi*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL],npoin2D_xi*sizeof(float));
                        memcpy(&buffer_send_xi_to_right[Z*npoin2D_xi],&accel_all_MPI_buffers[right_xi*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL],npoin2D_xi*sizeof(float));
// send messages forward along each row
                        if(iproc_xi == 0) {
                                sender = MPI_PROC_NULL;
                        } else {
                                sender = addressing[iproc_xi-1][iproc_eta];
                        }
                        if(iproc_xi == NPROC_XI-1) {
                                receiver = MPI_PROC_NULL;
                        } else {
                                receiver = addressing[iproc_xi+1][iproc_eta];
                        }
                        MPI_Issend(buffer_send_xi_to_right,NDIM*npoin2D_xi, MPI_REAL, receiver,0, FTI_COMM_WORLD, request_send_xi_to_right);
                        MPI_Irecv(buffer_recv_xi_from_left,NDIM*npoin2D_xi, MPI_REAL,sender,0,FTI_COMM_WORLD, request_recv_xi_from_left);

//////////// useless because erased when received                                memset(buffer_recv_xi_from_right,0,NDIM*npoin2D_xi*sizeof(float));
// all the slices copy their left face into the buffer
                                memcpy(&buffer_send_xi_to_left[X*npoin2D_xi],&accel_all_MPI_buffers[left_xi*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL],npoin2D_xi*sizeof(float));
                                memcpy(&buffer_send_xi_to_left[Y*npoin2D_xi],&accel_all_MPI_buffers[left_xi*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL],npoin2D_xi*sizeof(float));
                                memcpy(&buffer_send_xi_to_left[Z*npoin2D_xi],&accel_all_MPI_buffers[left_xi*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL],npoin2D_xi*sizeof(float));
// send messages backward along each row
                                if(iproc_xi == NPROC_XI-1) {
                                        sender = MPI_PROC_NULL;
                                } else {
                                        sender = addressing[iproc_xi+1][iproc_eta];
                                }
                                if(iproc_xi == 0) {
                                        receiver = MPI_PROC_NULL;
                                } else {
                                        receiver = addressing[iproc_xi-1][iproc_eta];
                                }
                                MPI_Issend(buffer_send_xi_to_left,NDIM*npoin2D_xi, MPI_REAL, receiver,0, FTI_COMM_WORLD, request_send_xi_to_left);
                                MPI_Irecv(buffer_recv_xi_from_right,NDIM*npoin2D_xi, MPI_REAL,sender,0,FTI_COMM_WORLD, request_recv_xi_from_right);

                *phase=1;
                return;
        }


        if (*phase == 1) {
// phase 1: perform the sum of the contributions along xi
// if some of the messages have not arrived yet, just return
                        MPI_Test(request_send_xi_to_right, &flag_send_xi_to_right, &status); if (!flag_send_xi_to_right) { return; }
                        MPI_Test(request_recv_xi_from_left, &flag_recv_xi_from_left, &status);; if (!flag_recv_xi_from_left) { return; }
                        MPI_Test(request_send_xi_to_left, &flag_send_xi_to_left, &status);; if (!flag_send_xi_to_left) { return; }
                        MPI_Test(request_recv_xi_from_right, &flag_recv_xi_from_right, &status);; if (!flag_recv_xi_from_right) { return; }

// all the slices add the buffer received to the contributions on the left face
                                if (iproc_xi > 0) {
                                        for (i=0;i<npoin2D_xi;i++) {
                                                accel_all_MPI_buffers[left_xi*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + i] += buffer_recv_xi_from_left[X*npoin2D_xi + i];
                                                accel_all_MPI_buffers[left_xi*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + i] += buffer_recv_xi_from_left[Y*npoin2D_xi + i];
                                                accel_all_MPI_buffers[left_xi*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + i] += buffer_recv_xi_from_left[Z*npoin2D_xi + i];
                                        }
                                        for (int j=0;j<npoin2D_eta;j++) {
// update the common points between buffers
                                                if (ibool1D_leftxi_lefteta[j]>=0) {
                                                        accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_left[X*npoin2D_xi + ibool1D_leftxi_lefteta[j]];
                                                        accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_left[Y*npoin2D_xi + ibool1D_leftxi_lefteta[j]];
                                                        accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_left[Z*npoin2D_xi + ibool1D_leftxi_lefteta[j]];
                                                }
                                                if (ibool1D_leftxi_righteta[j]>=0) {
                                                        accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_left[X*npoin2D_xi + ibool1D_leftxi_righteta[j]];
                                                        accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_left[Y*npoin2D_xi + ibool1D_leftxi_righteta[j]];
                                                        accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_left[Z*npoin2D_xi + ibool1D_leftxi_righteta[j]];
                                                }
                                        }
                                }

// add to right_xi
// all the slices copy the buffer received to the contributions on the right face
                                if(iproc_xi < NPROC_XI-1) {
                                        for(i=0;i<npoin2D_xi;i++) {
                                                accel_all_MPI_buffers[right_xi*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + i] += buffer_recv_xi_from_right[X*npoin2D_xi + i];
                                                accel_all_MPI_buffers[right_xi*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + i] += buffer_recv_xi_from_right[Y*npoin2D_xi + i];
                                                accel_all_MPI_buffers[right_xi*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + i] += buffer_recv_xi_from_right[Z*npoin2D_xi + i];
                                        }
                                        for (int j=0;j<npoin2D_eta;j++) {
                                                if (ibool1D_rightxi_lefteta[j]>=0) {
                                                        accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_right[X*npoin2D_xi + ibool1D_rightxi_lefteta[j]];
                                                        accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_right[Y*npoin2D_xi + ibool1D_rightxi_lefteta[j]];
                                                        accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_right[Z*npoin2D_xi + ibool1D_rightxi_lefteta[j]];
                                                }
                                                if (ibool1D_rightxi_righteta[j]>=0) {
                                                        accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_right[X*npoin2D_xi + ibool1D_rightxi_righteta[j]];
                                                        accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_right[Y*npoin2D_xi + ibool1D_rightxi_righteta[j]];
                                                        accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + j] += buffer_recv_xi_from_right[Z*npoin2D_xi + ibool1D_rightxi_righteta[j]];
                                                }
                                        }
                                }

                                *phase=2;
                                return;
        }


        if (*phase == 2) {
// phase 2: send right_eta to left_eta
//////////// useless because erased when received                        memset(buffer_recv_eta_from_down,0,NDIM*npoin2D_eta*sizeof(float));
// the slices copy their right face into the buffer
                        memcpy(&buffer_send_eta_to_up[X*npoin2D_eta],&accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL],npoin2D_eta*sizeof(float));
                        memcpy(&buffer_send_eta_to_up[Y*npoin2D_eta],&accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL],npoin2D_eta*sizeof(float));
                        memcpy(&buffer_send_eta_to_up[Z*npoin2D_eta],&accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL],npoin2D_eta*sizeof(float));
// send messages forward along each row
                        if(iproc_eta == 0) {
                                sender = MPI_PROC_NULL;
                        } else {
                                sender = addressing[iproc_xi][iproc_eta-1];
                        }
                        if(iproc_eta == NPROC_ETA-1) {
                                receiver = MPI_PROC_NULL;
                        } else {
                                receiver = addressing[iproc_xi][iproc_eta+1];
                        }
                        MPI_Issend(buffer_send_eta_to_up,NDIM*npoin2D_eta, MPI_REAL, receiver,0, FTI_COMM_WORLD, request_send_eta_to_up);
                        MPI_Irecv(buffer_recv_eta_from_down,NDIM*npoin2D_eta, MPI_REAL,sender,0,FTI_COMM_WORLD, request_recv_eta_from_down);

// the contributions are correctly assembled on the left side of each slice
// now we have to send the result back to the sender
// all slices copy the left face into the buffer
                                memcpy(&buffer_send_eta_to_down[X*npoin2D_eta],&accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL],npoin2D_eta*sizeof(float));
                                memcpy(&buffer_send_eta_to_down[Y*npoin2D_eta],&accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL],npoin2D_eta*sizeof(float));
                                memcpy(&buffer_send_eta_to_down[Z*npoin2D_eta],&accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL],npoin2D_eta*sizeof(float));
// send messages backward along each row
                                if(iproc_eta == NPROC_ETA-1) {
                                        sender = MPI_PROC_NULL;
                                } else {
                                        sender = addressing[iproc_xi][iproc_eta+1];
                                }
                                if(iproc_eta == 0) {
                                        receiver = MPI_PROC_NULL;
                                } else {
                                        receiver = addressing[iproc_xi][iproc_eta-1];
                                }
                                MPI_Issend(buffer_send_eta_to_down,NDIM*npoin2D_eta, MPI_REAL, receiver,0, FTI_COMM_WORLD, request_send_eta_to_down);
                                MPI_Irecv(buffer_recv_eta_from_up,NDIM*npoin2D_eta, MPI_REAL,sender,0,FTI_COMM_WORLD, request_recv_eta_from_up);

                *phase = 3;
                return;
        }


        if (*phase == 3) {
// phase 3: perform the sum of the contributions and send left_eta back to right_eta
// if some of the messages have not arrived yet, just return
                        MPI_Test(request_send_eta_to_up, &flag_send_eta_to_up, &status); if (!flag_send_eta_to_up) { return; }
                        MPI_Test(request_recv_eta_from_down, &flag_recv_eta_from_down, &status); if (!flag_recv_eta_from_down) { return; }
                        MPI_Test(request_send_eta_to_down, &flag_send_eta_to_down, &status); if (!flag_send_eta_to_down) { return; }
                        MPI_Test(request_recv_eta_from_up, &flag_recv_eta_from_up, &status); if (!flag_recv_eta_from_up) { return; }

// all slices add the buffer received to the contributions on the left face
                                if (iproc_eta > 0) {
                                        for (i=0;i<npoin2D_eta;i++) {
                                                accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + i] += buffer_recv_eta_from_down[X*npoin2D_eta + i];
                                                accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + i] += buffer_recv_eta_from_down[Y*npoin2D_eta + i];
                                                accel_all_MPI_buffers[left_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + i] += buffer_recv_eta_from_down[Z*npoin2D_eta + i];
                                        }
                                }

// add to right_eta
// all slices copy the buffer received to the contributions on the right face
                                if(iproc_eta < NPROC_ETA-1) {
                                        for(i=0;i<npoin2D_eta;i++) {
                                                accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + X*NGLOB2DMAX_ALL + i] += buffer_recv_eta_from_up[X*npoin2D_eta + i];
                                                accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Y*NGLOB2DMAX_ALL + i] += buffer_recv_eta_from_up[Y*npoin2D_eta + i];
                                                accel_all_MPI_buffers[right_eta*NDIM*NGLOB2DMAX_ALL + Z*NGLOB2DMAX_ALL + i] += buffer_recv_eta_from_up[Z*npoin2D_eta + i];
                                        }
                                }

                                *phase=4;
                                return;
                        }
        return;
}

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

void assemble_MPI_blocking_scalar(int* iboolleft_xi, int* iboolright_xi, int* iboolleft_eta, int* iboolright_eta, int iproc_xi, int iproc_eta, int addressing[NPROC_XI][NPROC_ETA], int npoin2D_xi, int npoin2D_eta, float* buffer_send_xi, float* buffer_recv_xi, float* buffer_send_eta, float* buffer_recv_eta, float rmass[NGLOB])
{

// IMPORTANT: after assembling, we first need to overwrite the values of the mass matrix "rmass"
// by the values of the "xi" buffers, and then by the values of the "eta" buffers.
// the order ("xi" first, then "eta") does matter because the grid points that
// belong to the common vertical edge are up to date in the "eta" buffers only
// and therefore the corresponding values in the "xi" buffers, which are not correct,
// must be discarded

        int sender,receiver,i;
        MPI_Status status;

// phase 0: send right_xi to left_xi
// assemble along xi only if more than one slice
// DK DK not needed any more // assemble along xi only if more than one slice
//      if (NPROC_XI > 1) {
// the slices copy their right face into the buffer
                for (i=0;i<npoin2D_xi;i++) {
                        buffer_send_xi[i] = rmass[iboolright_xi[i]];
                }
// send messages forward along each row
                if(iproc_xi == 0) {
                        sender = MPI_PROC_NULL;
                } else {
                        sender = addressing[iproc_xi-1][iproc_eta];
                }
                if(iproc_xi == NPROC_XI-1) {
                        receiver = MPI_PROC_NULL;
                } else {
                        receiver = addressing[iproc_xi+1][iproc_eta];
                }
                MPI_Sendrecv(buffer_send_xi, npoin2D_xi, MPI_REAL, receiver,0, buffer_recv_xi,npoin2D_xi, MPI_REAL,sender,0,FTI_COMM_WORLD, &status);


// phase 1: perform the sum of the contributions and send left_xi back to right_xi
// all the slices add the buffer received to the contributions on the left face
                if (iproc_xi > 0) {
                        for (i=0;i<npoin2D_xi;i++) {
                                rmass[iboolleft_xi[i]] += buffer_recv_xi[i];
                        }
                }
// the contributions are correctly assembled on the left side of each slice
// now we have to send the result back to the sender
// all the slices copy their left face into the buffer
                for (i=0;i<npoin2D_xi;i++) {
                        buffer_send_xi[i] = rmass[iboolleft_xi[i]];
                }
// send messages backward along each row
                if(iproc_xi == NPROC_XI-1) {
                        sender = MPI_PROC_NULL;
                } else {
                        sender = addressing[iproc_xi+1][iproc_eta];
                }
                if(iproc_xi == 0) {
                        receiver = MPI_PROC_NULL;
                } else {
                        receiver = addressing[iproc_xi-1][iproc_eta];
                }
                MPI_Sendrecv(buffer_send_xi, npoin2D_xi, MPI_REAL, receiver,0, buffer_recv_xi,npoin2D_xi, MPI_REAL,sender,0,FTI_COMM_WORLD,&status);


// phase 2: overwrite right_xi
// all the slices copy the buffer received to the contributions on the right face
                if(iproc_xi < NPROC_XI-1) {
                        for (i=0;i<npoin2D_xi;i++) {
                                rmass[iboolright_xi[i]] = buffer_recv_xi[i];
                        }
                }
//      }


// phase 3: send right_eta to left_eta
// DK DK not needed any more // assemble along xi only if more than one slice
//      if (NPROC_ETA > 1) {
// all the slices copy their right face into the buffer
                for (i=0;i<npoin2D_eta;i++) {
                        buffer_send_eta[i] = rmass[iboolright_eta[i]];
                }
// send messages forward along each row
                if(iproc_eta == 0) {
                        sender = MPI_PROC_NULL;
                } else {
                        sender = addressing[iproc_xi][iproc_eta-1];
                }
                if(iproc_eta == NPROC_ETA-1) {
                        receiver = MPI_PROC_NULL;
                } else {
                        receiver = addressing[iproc_xi][iproc_eta+1];
                }
                MPI_Sendrecv(buffer_send_eta, npoin2D_eta, MPI_REAL, receiver,0,buffer_recv_eta,npoin2D_eta, MPI_REAL,sender,0,FTI_COMM_WORLD,&status);


// phase 4: perform the sum of the contributions and send left_eta back to right_eta
// all the slices add the buffer received to the contributions on the left face
                if (iproc_eta > 0) {
                        for (i=0;i<npoin2D_eta;i++) {
                                rmass[iboolleft_eta[i]] += buffer_recv_eta[i];
                        }
                }
// the contributions are correctly assembled on the left side of each slice
// now we have to send the result back to the sender
// all slices copy the left face into the buffer
                for (i=0;i<npoin2D_eta;i++) {
                        buffer_send_eta[i] = rmass[iboolleft_eta[i]];
                }
// send messages backward along each row
                if(iproc_eta == NPROC_ETA-1) {
                        sender = MPI_PROC_NULL;
                } else {
                        sender = addressing[iproc_xi][iproc_eta+1];
                }
                if(iproc_eta == 0) {
                        receiver = MPI_PROC_NULL;
                } else {
                        receiver = addressing[iproc_xi][iproc_eta-1];
                }
                MPI_Sendrecv(buffer_send_eta, npoin2D_eta, MPI_REAL, receiver,0, buffer_recv_eta,npoin2D_eta, MPI_REAL,sender,0,FTI_COMM_WORLD,&status);


// phase 5: overwrite right_eta
// all the slices copy the buffer received to the contributions on the right face
                if(iproc_eta < NPROC_ETA-1) {
                        for (i=0;i<npoin2D_eta;i++) {
                                rmass[iboolright_eta[i]] = buffer_recv_eta[i];
                        }
                }
//      }
        return;
}

#endif  // of ACTUALLY_ASSEMBLE_MPI_SLICES

// ***************************************************************************
// ***************************************************************************
// ***************************************************************************

// print CUDA errors thrown by CUDA function calls.
// output an error message and exit
void print_CUDA_error_if_any(cudaError_t err, int num)
{
  if (cudaSuccess != err)
  {
    printf("\nCUDA error !!!!! <%s> !!!!! at CUDA call # %d\n",cudaGetErrorString(err),num);
    fflush(stdout);
#ifdef USE_MPI
    MPI_Abort(FTI_COMM_WORLD,1);
#endif
    exit(0);
  }
  return;
}

