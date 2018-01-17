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
! the Free Software Foundation; either version 3 of the License, or
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

#include "config.h"

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <strstream>

using namespace std;

// macros
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) < (y) ? (y) : (x))
#define exit_error(msg) { fputs (msg,stderr); exit(1); }

// debugging
//#define TRACE(x) printf("%s\n",x);
//#define TRACE(x) printf("thread %lu: %s\n",(unsigned long) pthread_self(),x);
#define TRACE(x)

/* ----------------------------------------------------------------------------------------------- */

#ifdef HAVE_VTK

#pragma message ("\nCompiling with: HAVE_VTK enabled\n")

#include <vtkVersion.h>

// VTK version compilation infos
#if VTK_MAJOR_VERSION <= 5
// VTK version 5.x
#pragma message ("\nCompiling with: VTK major version <= 5\n")

#else
// VTK version 6+
// special apple
#if defined (__APPLE__) || defined(MACOSX)
// because main thread must run GUI... this would require a work-around not implemented yet
#pragma message ("\nSorry VTK version 6+ on Apple platform won't work yet..., please use a VTK version 5.x\n")
#endif

// auto-initialization:
// see http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Factories_now_require_defines
#if VTK_MAJOR_VERSION == 6 && VTK_MINOR_VERSION == 0
// VTK version 6.0
// early versions of vtk 6 won't recognize VTK_MODULE_INIT(..)
#pragma message ("\nCompiling with: VTK version 6.0\n")
#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
//#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL2)
//#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#else
// VTK version 6+
// newer versions auto-initialization
#pragma message ("\nCompiling with: VTK major version 6+\n")
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkInteractionStyle); // needed to get vtkInteractionStyle instance
VTK_MODULE_INIT(vtkRenderingFreeType); // needed to get TextRenderer instance
VTK_MODULE_INIT(vtkRenderingFreeTypeOpenGL);
// by default VTK 7 uses OpenGL2 rendering; add -lvtkRenderingOpenGL2
#if VTK_MAJOR_VERSION == 7
#pragma message ("\nCompiling with: VTK major version 7 -> using OpenGL2 as default rendering\n")
VTK_MODULE_INIT(vtkRenderingOpenGL2);
//VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
#else
VTK_MODULE_INIT(vtkRenderingOpenGL);
//VTK_MODULE_INIT(vtkRenderingVolumeOpenGL);
#endif
#endif

#endif

#include <vtkActor.h>
#include <vtkAppendFilter.h>
#include <vtkBoundingBox.h>
#include <vtkBox.h>
#include <vtkBoxClipDataSet.h>
#include <vtkCamera.h>
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkClipDataSet.h>
#include <vtkClipPolyData.h>
#include <vtkClipVolume.h>
#include <vtkCommand.h>
#include <vtkContourFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkDoubleArray.h>
#include <vtkExtractGeometry.h>
#include <vtkFloatArray.h>
#include <vtkGeometryFilter.h>
#include <vtkGlyph2D.h>
#include <vtkGlyphSource2D.h>
#include <vtkGlyph3D.h>
#include <vtkHexahedron.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkJPEGWriter.h>
#include <vtkLabeledDataMapper.h>
#include <vtkLabelHierarchy.h>
#include <vtkLabelPlacer.h>
#include <vtkLookupTable.h>
#include <vtkLegendBoxActor.h>
#include <vtkMath.h>
#include <vtkMaskPoints.h>
#include <vtkMutexLock.h>
#include <vtkPlane.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPointSetToLabelHierarchy.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkQuad.h>
#include <vtkRegularPolygonSource.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkSetGet.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStringArray.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkTimerLog.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLUnstructuredGridWriter.h>


// VTK with parallel support: -lvtkParallel
#ifdef VTK_VIS_PARALLEL
#include <vtkCompositeRenderManager.h>
#include <vtkMPIController.h>
#endif

/* ----------------------------------------------------------------------------------------------- */


/******************************************************************************/

// USER PARAMETERS
int CONTOUR_FREESURFACE = 1;  // creates contour for free surface topography
int CONTOUR_VOLUME = 1;       // creates contour for volumetric data

// color range to color scale plotting
double gcolor_min  = 0.0;     // minimum value
double gcolor_max  = 0.005;   // maximum value
double gcolor_incr = 0.001;   // increments

// if set to 1 it uses an extract filter, otherwise it will use clip filters
static int use_vtk_extract = 0;

// exagerates vertical dimension (1.0 == no exageration, 1.5 == 50% exageration...) for free surface
static float vertical_exageration = 1.5;

// receiver labels
// note: we will create receiver labels although this can slow down rendering quite a bit.
//       for now, we will put visibility off by default. users will have to interact to show it.
// limit number of receiver labels
static bool use_reduced_labeling = false;
// if reduced flag is true, this is the maximum number of labels used
static int max_receiver_labels = 20;

/******************************************************************************/


// global flags to check if routine was actually called
int SHOW_FREESURFACE = 0;
int SHOW_VOLUMEDATA = 0;
int SHOW_GPU_TEXT = 0;

// function definitions
void wait_vtk_thread();
void* init_vtk_thread(void*);
void create_vtk_thread();
void terminate_vtk_thread();

void interactor_usage();
void save_snapshot_vtu();
void save_snapshot_jpg();
void set_color_scale(int);
void abort_simulation();
void init_vtk();
void init_vtk_window_interactor();
void clean_vtk_arrays();
void sync_rendering(int sec = 20);
void sync_halt_simulation();

extern "C" void FC_FUNC_(initialize_vtkwindow,INITIALIZE_VTKWINDOW)(int*);
extern "C" void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float*,float*,float*);
extern "C" void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int*,float*,float*,float*,int*,int*);
extern "C" void FC_FUNC_(prepare_vtkreceivers,PREPARE_VTKRECEIVERS)(int*, float*, float*, float*, int*, char*);
extern "C" void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int*,float*,float*,float*,int*,int*);
extern "C" void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int*,float*,float*);
extern "C" void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)(int*);

/* ----------------------------------------------------------------------------------------------- */


class Visualization {
  public:
    // 2D surface data
    vtkPoints* points2D;
    vtkFloatArray* data_array2D;

    vtkUnstructuredGrid* volume2D;

    vtkDataSetMapper* mapMesh2D;
    vtkActor* actor2D;

    // 3D volume data
    vtkPoints* points3D;
    vtkFloatArray* data_array3D;

    vtkUnstructuredGrid* volume3D;

    vtkDataSetMapper* mapMesh3D;
    vtkActor* actor3D;

    // clipping
    vtkPlane* clipPlane1;
    vtkPlane* clipPlane2;

    vtkTableBasedClipDataSet* clip1;
    vtkTableBasedClipDataSet* clip2;
    vtkTableBasedClipDataSet* clip1m;
    //vtkClipDataSet* clip1;
    //vtkClipDataSet* clip2;
    //vtkClipDataSet* clip1m;
    vtkAppendFilter* merger;

    vtkExtractGeometry* extract;
    vtkDataSetSurfaceFilter* clippersurface;

    // mesh cleaning
    //vtkGeometryFilter* geometryFilter;
    //vtkPolyData* polydata;
    //vtkCleanPolyData* cleanFilter;

    // colors
    vtkLookupTable* lut2D;
    vtkLookupTable* lut;
    int icolor;
    vtkScalarBarActor* legendcolor;

    // window text
    vtkTextActor* text;
    vtkTextActor* textGPU;
    vtkTextActor* help;

    // user interaction
    volatile int jpegImageOn;
    volatile int colorAdjustOn;
    volatile int bgBlackOn;
    volatile int haltOn;
    volatile int do_restart;

    // rendering
    vtkRenderer* ren;
    vtkRenderWindow* renWin;
    vtkRenderWindowInteractor* iren;

    vtkMutexLock *mutex;
    volatile int do_render;

    // camera
    vtkCamera* camera;

    double pcam[3]; // camera position
    double rclip[2]; // clipping range
    double pfocal[3]; // focal point

    vtkInteractorStyleTrackballCamera *style;

    // source sphere
    vtkSphereSource* sphere;
    vtkPolyDataMapper* mapperSphere;
    vtkActor* actorSphere;
    double pos_source[3]; // source location

    // receiver glyphs
    //vtkGlyph2D* glyph;
    vtkGlyph3D* glyph;
    vtkPolyDataMapper* glyphMapper;
    vtkActor* glyphActor;

    vtkSelectVisiblePoints* labelVis;
    vtkLabeledDataMapper* labelMapper;
    vtkActor2D* labelActor;

    // countours
    vtkContourFilter* contour;
    vtkPolyDataMapper* contourMapper;
    vtkActor* contourActor;

    vtkContourFilter* contour3D;
    vtkPolyDataMapper* contour3DMapper;
    vtkActor* contour3DActor;
};

class VTKmesh{
  public:
    int np;
    int nspec;
    double bounds[6];
};

class VTKState {
  public:
    // vtk rendering
    Visualization vtk;

    // meshes
    VTKmesh freesurface;
    VTKmesh volume;

    // timing
    vtkTimerLog* timer;
};

// global vtk state variable
static VTKState fs;

static string info_string =
  "interactor usage, press:\n"
  "  a            - automatically adjust color scale value (toggle on/off) \n"
  "  <up>/<down>  - change color scale maximum \n"
  "  b            - change background color (toggle white/black) \n"
  "  c            - change color scheme \n"
  "  h            - help (toggle on/off) \n"
  "  mouse click  - move/rotate view \n"
  "  r            - reset view \n"
  "\n"
  "  v            - save snapshot as vtu file\n"
  "  i            - save snapshot as image (jpg-file) (toggle on/off)\n"
  "\n"
  "  5            - show volume (toggle on/off)\n"
  "  6            - show volume contour (toggle on/off)\n"
  "  7            - show freesurface (toggle on/off)\n"
  "  8            - show freesurface contour (toggle on/off)\n"
  "  9            - show source/receivers glyphs (toggle on/off)\n"
  "  0            - show receiver labels (toggle on/off)\n"
  "\n"
  "  <space>      - halt/continue simulation\n"
  "  <escape>,q,e - abort and exit simulation\n"
  "  y            - yes, restart simulation (after the end of the time loop)\n"
  ;

static string info_string_defaults =
  "\n"
  "default interactions: \n"
  "  Mouse bindings:\n"
  "    camera: Button 1 - rotate\n"
  "            Button 2 - pan\n"
  "            Button 3 - zoom\n"
  "            ctrl-Button 1 - spin\n"
  "    actor:  Button 1 - rotate\n"
  "            Button 2 - pan\n"
  "            Button 3 - uniform scale\n"
  "            ctrl-Button 1 - spin\n"
  "            ctrl-Button 2 - dolly\n"
  "\n"
  "  Keyboard bindings (upper or lower case):\n"
  "    j - joystick like mouse interactions\n"
  "    t - trackball like mouse interactions\n"
  "    o - object/ actor interaction\n"
  "    c - camera interaction\n"
  "    r - reset camera view\n"
  "    w - turn all actors wireframe\n"
  "    s - turn all actors surface\n"
  "    u - execute user defined function\n"
  "    p - pick actor under mouse pointer (if pickable)\n"
  "    3 - toggle in/out of 3D mode (if supported by renderer)\n"
  "    e - exit\n"
  "    q - exit\n"
  ;

// time step it from simulation
int global_it_step = 0;

// array for storing peak values (using a simple float array to avoid memory overhead from vtkFloatArray...)
float* global_data_array_peak;

/* ----------------------------------------------------------------------------------------------- */

// threading version

/* ----------------------------------------------------------------------------------------------- */

#include <pthread.h>

struct threadInfo{
  volatile bool finished;
  volatile bool started;
  volatile bool init;
  pthread_mutex_t mutex;
};

struct threadInfo vtk_thread_info;

pthread_t vtk_thread;
pthread_t main_thread;

// Waits until thread is finished with rendering
void wait_vtk_thread() {
  TRACE("wait_vtk_thread");

  int rc;

  // checks if thread still runs
  if (vtk_thread_info.started == true ) {
    void* status;
    rc = pthread_join(vtk_thread, &status);
    if (rc ) {
      printf("Error; return code from pthread_join() is %d\n", rc);
      exit_error("Error in wait_vtk_thread: thread_join failed");
    }

    // checks finished flag
    if (vtk_thread_info.finished == false){
      // thread has completed, but somehow it isn't finished?
      printf("vtk_thread has not finished, but been cancelled\n");
      pthread_mutex_lock(&vtk_thread_info.mutex);
      vtk_thread_info.finished = true;
      pthread_mutex_unlock(&vtk_thread_info.mutex);
    };

    // reset
    pthread_mutex_lock(&vtk_thread_info.mutex);
    vtk_thread_info.started = false;
    pthread_mutex_unlock(&vtk_thread_info.mutex);
  }
}

// dummy argument is needed to avoid compiler warning
void* init_vtk_thread(void* dummy){
  TRACE("init_vtk_thread");
  int rc;

  rc = pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  if (rc != 0) exit_error("Error in init_vtk_thread: thread_cancel state failed");

  // setup vtk window only
  // note: starting the interactor will not return until it has been terminated via TerminateApp()
  init_vtk_window_interactor();

  // exits thread
  TRACE("interactor done, finishing vtk thread");
  pthread_mutex_lock(&vtk_thread_info.mutex);
  vtk_thread_info.finished = true;
  pthread_mutex_unlock(&vtk_thread_info.mutex);

  pthread_exit(NULL);

  TRACE("init_vtk_thread done");
  return NULL;  // Never used, but remove warnings.
}


// creates thread for reading adjoint sources
void create_vtk_thread() {
  TRACE("create_vtk_thread");

  int rc;

  // gets main thread id
  main_thread = pthread_self();

  // initializes thread info
  pthread_mutex_lock(&vtk_thread_info.mutex);
  vtk_thread_info.started = false;
  vtk_thread_info.finished = false;
  vtk_thread_info.init = true;
  pthread_mutex_unlock(&vtk_thread_info.mutex);

  // prepares the thread
  pthread_attr_t attr;

  rc = pthread_attr_init(&attr);
  if (rc != 0 ) exit_error("thread: initialization of thread failed");

  rc = pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  if (rc != 0 ) exit_error("thread: setting thread state failed");

  // sets new thread info
  pthread_mutex_lock(&vtk_thread_info.mutex);
  vtk_thread_info.started = true;
  vtk_thread_info.finished = false;
  pthread_mutex_unlock(&vtk_thread_info.mutex);

  // create and launch the thread.
  rc = pthread_create(&vtk_thread, &attr, init_vtk_thread, NULL);
  if (rc != 0 ){
    printf("ERROR: return code from pthread_create() is %d\n", rc);
    exit_error("thread: creating thread failed");
  }

  // thread infos
  printf("     main thread id: %lu\n",(unsigned long) main_thread);
  printf("     vtk  thread id: %lu\n\n",(unsigned long) vtk_thread);

  TRACE("create_vtk_thread done");
}

// terminates thread by sending a cancel signal
void terminate_vtk_thread(){
  TRACE("terminate_vtk_thread");

  int rc;

  // check if we are called from same vtk thread, then just stop interactor
  // tries to cancel thread
  if (pthread_equal(vtk_thread,pthread_self())){
    // stop interactor
    fs.vtk.iren->TerminateApp();
  }else{
    // stop interactor
    fs.vtk.mutex->Lock();
    fs.vtk.do_render = 2;
    fs.vtk.mutex->Unlock();
    // waits
    sync_rendering();

    // sends a cancel signal to other thread
    rc = pthread_cancel(vtk_thread);
    if (rc != 0) exit_error("thread: canceling thread failed");

    // waits until thread finished
    wait_vtk_thread();

    // exits thread
    pthread_mutex_lock(&vtk_thread_info.mutex);
    vtk_thread_info.finished = true;
    pthread_mutex_unlock(&vtk_thread_info.mutex);
  }

  TRACE("terminate_vtk_thread done");
}

/* ----------------------------------------------------------------------------------------------- */

// callback functions

/* ----------------------------------------------------------------------------------------------- */


// A class not derived from vtkObjectBase

class MyInteractor{
  public:
    void KeypressCallbackFunction(vtkObject* caller,long unsigned int eventId,void* callData ){
      // note: this callback function will only be executed by the vtk thread which called the interactor
      string key = fs.vtk.iren->GetKeySym();
      //cout << "Pressed: " << key << endl;
      // overrides default behaviour (which would stop the window interactor)
      if (key == "e") key = "Escape";
      if (key == "q") key = "Escape";
      // usage
      if (key == "h"){
        interactor_usage();
        // displays help text
        if( fs.vtk.help->GetVisibility() == 1){
          fs.vtk.help->SetVisibility( 0 );
        }else{
          fs.vtk.help->SetVisibility( 1 );
        }
        // update
        fs.vtk.renWin->Render();
      }
      // changes color scales
      if (key == "c"){
        fs.vtk.mutex->Lock();
        fs.vtk.icolor++;
        fs.vtk.mutex->Unlock();
        set_color_scale(fs.vtk.icolor);
        // update
        fs.vtk.renWin->Render();
      }
      // adjust color scaling
      if (key == "a"){
        fs.vtk.mutex->Lock();
        fs.vtk.colorAdjustOn++;
        fs.vtk.colorAdjustOn = fs.vtk.colorAdjustOn % 2 ;
        fs.vtk.mutex->Unlock();
        if( fs.vtk.colorAdjustOn ){
          printf("\ntoggle on: color adjusting\n");
        }else{
          printf("\ntoggle off: color adjusting off\n");
        }
        set_color_scale(fs.vtk.icolor);
        // update
        fs.vtk.renWin->Render();
      }
      // saves vtu snapshot
      if (key == "v"){
        save_snapshot_vtu();
      }
      // saves jpg snapshot
      if (key == "i"){
        // toggles jpeg output flag
        if( fs.vtk.jpegImageOn == 0 ){
          printf("\ntoggle on: save snapshot as jpeg-image\n");
          fs.vtk.mutex->Lock();
          fs.vtk.jpegImageOn = 1;
          fs.vtk.mutex->Unlock();
        }else{
          printf("\ntoggle off: save snapshot as jpeg-image\n");
          fs.vtk.mutex->Lock();
          fs.vtk.jpegImageOn = 0;
          fs.vtk.mutex->Unlock();
        }
        // update
        fs.vtk.renWin->Render();
        save_snapshot_jpg();
      }
      // stops/continues running simulation
      if (key == "space"){
        // toggles halt/continue flag
        if( fs.vtk.haltOn == 0 ){
          printf("\ntoggle on: halting simulation\n");
          fs.vtk.mutex->Lock();
          fs.vtk.haltOn = 1;
          fs.vtk.mutex->Unlock();
          // updates text
          fs.vtk.text->SetInput( "halting simulation, press <space> to continue ..." );
          TRACE("halting done");
        }else{
          printf("\ntoggle off: continue simulation\n");
          fs.vtk.mutex->Lock();
          fs.vtk.haltOn = 0;
          fs.vtk.mutex->Unlock();
          // updates text
          fs.vtk.text->SetInput( "continue simulation..." );
          TRACE("continue done");
          // turn help off
          if( fs.vtk.help->GetVisibility() == 1){
            fs.vtk.help->SetVisibility( 0 );
          }
        }
        // update
        fs.vtk.renWin->Render();
      }
      // resets view
      if (key == "r"){
        // reposition the camera, so that actor can be fully seen
        // range
        fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
        // focal point
        fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
        // position
        fs.vtk.camera->SetPosition( fs.vtk.pcam );
        // view
        fs.vtk.camera->SetViewAngle( 30.0 );
        fs.vtk.camera->SetViewUp( 0.0, 0.0, 1.0 );
        // reposition the camera, so that actor can be fully seen
        fs.vtk.ren->ResetCamera();
        fs.vtk.ren->ResetCameraClippingRange();
        // update
        fs.vtk.renWin->Render();
      }
      // changes color scale maximum values
      if(string(fs.vtk.iren->GetKeySym()) == "Up") {
        // increases color max by 10%
        fs.vtk.mutex->Lock();
        gcolor_max = gcolor_max + gcolor_incr;
        fs.vtk.mutex->Unlock();
        fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
        fs.vtk.lut->Build();
        // update
        fs.vtk.renWin->Render();
      }
      if(string(fs.vtk.iren->GetKeySym()) == "Down") {
        // decreases color max by 10%
        fs.vtk.mutex->Lock();
        gcolor_max = gcolor_max - gcolor_incr;
        if(gcolor_max < gcolor_min ) gcolor_max = gcolor_min;
        fs.vtk.mutex->Unlock();
        fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
        fs.vtk.lut->Build();
        // update
        fs.vtk.renWin->Render();
      }
      // changes background
      if (key == "b"){
        // toggles background
        if( fs.vtk.bgBlackOn == 0 ){
          printf("\ntoggle on: background black/white\n");
          fs.vtk.mutex->Lock();
          fs.vtk.bgBlackOn = 1;
          fs.vtk.mutex->Unlock();
          fs.vtk.ren->SetBackground(0,0,0); // Background color black
        }else{
          printf("\ntoggle off: background black/white\n");
          fs.vtk.mutex->Lock();
          fs.vtk.bgBlackOn = 0;
          fs.vtk.mutex->Unlock();
          fs.vtk.ren->SetBackground(1,1,1); // Background color white
        }
        //update
        fs.vtk.renWin->Render();
      }

      // toggles volume visibility
      if (key == "5"){
        if( SHOW_VOLUMEDATA == 1 ){
          if( fs.vtk.actor3D->GetVisibility() == 1){
            fs.vtk.actor3D->SetVisibility( 0 );
            printf("\ntoggle off: volume\n");
          }else{
            fs.vtk.actor3D->SetVisibility( 1 );
            printf("\ntoggle on: volume\n");
          }
          // update
          fs.vtk.renWin->Render();
        }
      }
      // toggles volume contour visibility
      if (key == "6"){
        if( SHOW_VOLUMEDATA == 1 && CONTOUR_VOLUME == 1 ){
          if( fs.vtk.contour3DActor->GetVisibility() == 1){
            fs.vtk.contour3DActor->SetVisibility( 0 );
            printf("\ntoggle off: volume contour\n");
          }else{
            fs.vtk.contour3DActor->SetVisibility( 1 );
            printf("\ntoggle on: volume contour\n");
          }
          // update
          fs.vtk.renWin->Render();
        }
      }
      // toggles freesurface visibility
      if (key == "7"){
        if( SHOW_FREESURFACE == 1 ){
          if( fs.vtk.actor2D->GetVisibility() == 1){
            fs.vtk.actor2D->SetVisibility( 0 );
            printf("\ntoggle off: free surface\n");
          }else{
            fs.vtk.actor2D->SetVisibility( 1 );
            printf("\ntoggle on: free surface\n");
          }
          // update
          fs.vtk.renWin->Render();
        }
      }
      // toggles freesurface contour visibility
      if (key == "8"){
        if( SHOW_FREESURFACE == 1 && CONTOUR_FREESURFACE == 1 ){
          if( fs.vtk.contourActor->GetVisibility() == 1){
            fs.vtk.contourActor->SetVisibility( 0 );
            printf("\ntoggle off: free surface contour\n");
          }else{
            fs.vtk.contourActor->SetVisibility( 1 );
            printf("\ntoggle on: free surface contour\n");
          }
          // update
          fs.vtk.renWin->Render();
        }
      }
      // toggles source glyph visibility
      if (key == "9"){
        if( fs.vtk.actorSphere->GetVisibility() == 1){
          fs.vtk.actorSphere->SetVisibility( 0 );
          fs.vtk.glyphActor->VisibilityOff();
          printf("\ntoggle off: source/receivers glyph\n");
        }else{
          fs.vtk.actorSphere->SetVisibility( 1 );
          fs.vtk.glyphActor->VisibilityOn();
          printf("\ntoggle on: source/receivers glyph\n");
        }
        // update
        fs.vtk.renWin->Render();
      }
      // toggles receiver labels visibility
      if (key == "0"){
        if( fs.vtk.labelActor->GetVisibility() == 1){
          fs.vtk.labelActor->VisibilityOff();
          printf("\ntoggle off: receiver labels\n");
        }else{
          fs.vtk.labelActor->VisibilityOn();
          printf("\ntoggle on: receiver labels\n");
        }
        fs.vtk.labelVis->Update();
        std::cout << "There are currently: " << fs.vtk.labelVis->GetOutput()->GetNumberOfPoints()
                  << " stations visible." << std::endl;
        // update
        fs.vtk.renWin->Render();
      }

      // exit by escape
      if (key == "Escape"){
        // updates text
        fs.vtk.text->SetInput( "...exiting simulation " );
        // update
        fs.vtk.renWin->Render();
        // abort simulation
        abort_simulation();
      }
      // changes background
      if (key == "y"){
        // toggles restart flag
        if( fs.vtk.do_restart == 0 ){
          printf("\ntoggle on: restart simulation\n");
          if (fs.vtk.haltOn == 1){
            fs.vtk.text->SetInput( "...restart simulation, press <space> to continue " );
          }else{
            fs.vtk.text->SetInput( "...restart simulation " );
          }
          fs.vtk.mutex->Lock();
          fs.vtk.do_restart = 1;
          fs.vtk.mutex->Unlock();
        }else{
          printf("\ntoggle off: stop, don't restart simulation\n");
          if (fs.vtk.haltOn == 1){
            fs.vtk.text->SetInput( "...don't restart simulation, press <space> to continue " );
          }else{
            fs.vtk.text->SetInput( "...don't restart simulation " );
          }
          fs.vtk.mutex->Lock();
          fs.vtk.do_restart = 0;
          fs.vtk.mutex->Unlock();
        }
        //update
        fs.vtk.renWin->Render();
      }
    }
};

static MyInteractor mykey;

class CommandSubclass2 : public vtkCommand {
  public:
    vtkTypeMacro(CommandSubclass2, vtkCommand);

    // note: if compilation error state something like undefined symbols because of vtable reference in vtkCommand
    //       you are compiling with a different compiler than vtk libraries were created (e.g. using gcc here, but clang for vtk library)
    // problematic references are:
    //  virtual void PrintHeader(ostream &os, vtkIndent indent){};
    //  virtual void PrintTrailer(ostream& os, vtkIndent indent){};
    //  virtual void PrintSelf(ostream& os, vtkIndent indent){};
    //
    // please use the same compiler for both and this error should be resolved.

    static CommandSubclass2 *New(){
      return new CommandSubclass2;
    }

    void Execute(vtkObject *vtkNotUsed(caller), unsigned long vtkNotUsed(eventId), void *vtkNotUsed(callData)){
    //void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId), void *vtkNotUsed(callData)){
      //vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
      //
      // note: this execute function will be executed only by the vtk thread which also created the interactor
      // checks for rendering
      if (fs.vtk.do_render == 1){
        TRACE("execute: rendering");
        // render window
        fs.vtk.renWin->Render();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: rendering done");
      }
      // check for terminating interactive
      if (fs.vtk.do_render == 2){
        TRACE("execute: terminating interactor");
        // stop interactor
        fs.vtk.iren->TerminateApp();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: terminate done");
      }
      // render image
      if (fs.vtk.do_render == 3){
        TRACE("execute: rendering image");
        // update
        fs.vtk.renWin->Render();
        // save image
        save_snapshot_jpg();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: rendering image done");
      }
      // render with resetting camera
      if (fs.vtk.do_render == 4){
        TRACE("execute: rendering with camera reset");
        // reposition the camera, so that actor can be fully seen
        fs.vtk.ren->ResetCamera();
        fs.vtk.ren->ResetCameraClippingRange();
        // update
        fs.vtk.renWin->Render();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: rendering with camera reset done");
      }
    }
};


/* ----------------------------------------------------------------------------------------------- */

// helper function

/* ----------------------------------------------------------------------------------------------- */

void interactor_usage(){
  cout << endl;
  cout << "************************************************************************\n" << endl;
  cout << info_string << endl;
  cout << info_string_defaults << endl;
  cout << "\n************************************************************************" << endl;
  cout << endl;
}

// Write vtu file
void save_snapshot_vtu(){
  TRACE("save_snapshot_vtu");

  char filename[180];
  if( global_it_step > 0 ){
    sprintf(filename,"test_snapshot.%6.6d.vtu",global_it_step);
  }else{
    sprintf(filename,"test_snapshot.vtu");
  }

  // creates writer
  vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(fs.vtk.volume3D);
#else
  writer->SetInputData(fs.vtk.volume3D);
#endif
  writer->SetDataModeToAscii();
  writer->Write();

  //clean up
  writer->Delete();

  printf("snapshot written to file: %s\n\n",filename);
}


// Write jpg file
void save_snapshot_jpg(){
  TRACE("save_snapshot_jpg");

  //std::string filename = "test_snapshot.jpg";

  char filename[180];
  if( global_it_step > 0 ){
    sprintf(filename,"test_snapshot.%6.6d.jpg",global_it_step);
  }else{
    sprintf(filename,"test_snapshot.jpg");
  }

  // window filter
  vtkWindowToImageFilter* w2i = vtkWindowToImageFilter::New();
  w2i->SetInput(fs.vtk.renWin);
  w2i->Update();

  // creates writer
  vtkJPEGWriter* writer = vtkJPEGWriter::New();
  //writer->SetFileName(filename.c_str());
  writer->SetFileName(filename);
  writer->SetInputConnection(w2i->GetOutputPort());
  writer->Write();

  //clean up
  writer->Delete();
  w2i->Delete();

  printf("snapshot written to file: %s\n\n",filename);
}

// changes color scheme
void set_color_scale(int icolor){
  TRACE("set_color_scale");

  if( icolor == 0 ){
    // sets (default) rainbow color scale
    //# Set the hue range: from low to high the
    //# table passes through blue, green, yellow,
    //# orange, and red
    gcolor_min = 0.0;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLinear();

    fs.vtk.lut->SetValueRange( 0.6, 1.0 );
    fs.vtk.lut->SetHueRange( 0.66667, 0.0 );
    fs.vtk.lut->SetSaturationRange( 1.0, 1.0 );
  }else if( icolor == 1 ){
    // red-blue color scale
    // Since the default lut is
    // a rainbow lut, we only have
    // to worry about the hue range.
    gcolor_min = 0.0;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLinear();

    fs.vtk.lut->SetValueRange( 0.6, 1.0 );
    fs.vtk.lut->SetHueRange( 0.0, 0.667 );
    fs.vtk.lut->SetSaturationRange( 1.0, 1.0 );
  }else if( icolor == 2 ){
    // red color scale
    gcolor_min = 1.e-3*gcolor_max;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLog10();

    fs.vtk.lut->SetValueRange( 0.4, 1.0 );
    fs.vtk.lut->SetHueRange( 0.0, 0.4 );
    fs.vtk.lut->SetSaturationRange( 0.5, 0.0 );
  }else if( icolor == 3 ){
    // blue color scale
    gcolor_min = 1.e-3*gcolor_max;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLog10();

    fs.vtk.lut->SetValueRange( 0.4, 1.0 );
    fs.vtk.lut->SetHueRange( 0.6, 0.6 );
    fs.vtk.lut->SetSaturationRange( 0.5, 0.0 );
  }else if( icolor == 4 ){
    // black color scale
    gcolor_min = 1.e-3*gcolor_max;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLog10();

    fs.vtk.lut->SetValueRange( 0.4, 1.0 ); // from black to white
    fs.vtk.lut->SetHueRange( 0.0, 1.0 );
    fs.vtk.lut->SetSaturationRange( 0.0, 0.0 ); // no color saturation
  }else{
    // reset
    fs.vtk.icolor = 0;
    // sets (default) rainbow color scale
    //# Set the hue range: from low to high the
    //# table passes through blue, green, yellow,
    //# orange, and red
    gcolor_min = 0.0;
    fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
    fs.vtk.lut->SetScaleToLinear();

    fs.vtk.lut->SetValueRange( 0.6, 1.0 );
    fs.vtk.lut->SetHueRange( 0.66667, 0.0 );
    fs.vtk.lut->SetSaturationRange( 1.0, 1.0 );
  }
  fs.vtk.lut->Build();
}


void abort_simulation(){
  TRACE("abort_simulation");

  // tries to cancel thread / or stop interactor

  // finish with stop - not working if called within same vtk thread
  //FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)();

  // finish interactor thread
  terminate_vtk_thread();
  sleep(1);

  // clean up
  clean_vtk_arrays();

  // abort program
  exit_error("simulation terminates by escape key\n");
}

void init_vtk(){
  TRACE("init_vtk");

  // initializes vtk setup
  printf("vtk: setting up...\n");

  // initializes flags
  fs.vtk.haltOn = 0;
  fs.vtk.jpegImageOn = 0;
  fs.vtk.do_restart = 0;
  fs.vtk.do_render = 0;
  for (int i = 0; i < 6; i++) fs.freesurface.bounds[i] = 0.0;
  for (int i = 0; i < 6; i++) fs.freesurface.bounds[i] = 0.0;

  // timer
  TRACE("     timer\n");
  fs.timer = vtkTimerLog::New();
  if (fs.timer == NULL){ exit_error("Error: invalid timer, could not create object\n"); }

  // mutex
  TRACE("     mutex\n");
  fs.vtk.mutex = vtkMutexLock::New();
  if (fs.vtk.mutex == NULL){ exit_error("Error: invalid mutex, could not create object\n"); }

  // camera position
  TRACE("     camera\n");
  fs.vtk.camera = vtkCamera::New();
  if (fs.vtk.camera == NULL){ exit_error("Error: invalid camera, could not create object\n"); }

  // initial settings
  fs.vtk.rclip[0] = -1.5;
  fs.vtk.rclip[1] = 1.5;
  fs.vtk.pfocal[0] = 0.0;
  fs.vtk.pfocal[1] = 0.0;
  fs.vtk.pfocal[2] = 0.0;
  fs.vtk.pcam[0] = 1.0;
  fs.vtk.pcam[1] = 1.0;
  fs.vtk.pcam[2] = 0.3;

  // range
  fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
  // camera focal point
  fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
  // camera position
  fs.vtk.camera->SetPosition( fs.vtk.pcam );
  // view
  fs.vtk.camera->SetViewAngle( 30.0 );
  fs.vtk.camera->SetViewUp( 0.0, 0.0, 1.0 );

  // renderer
  TRACE("     renderer\n");
  fs.vtk.ren = vtkRenderer::New();
  if (fs.vtk.ren == NULL){ exit_error("Error: invalid renderer, could not create object\n"); }

  // sets camera
  fs.vtk.ren->SetActiveCamera(fs.vtk.camera);

  // Background color white
  fs.vtk.bgBlackOn = 0;
  fs.vtk.ren->SetBackground(1,1,1);

  // text actors
  TRACE("     text\n");
  // GPU flag
  fs.vtk.textGPU = vtkTextActor::New();
  if (fs.vtk.textGPU == NULL){ exit_error("Error: invalid textActor, could not create object\n"); }

  fs.vtk.textGPU->GetTextProperty()->SetFontSize ( 32 );
  fs.vtk.textGPU->SetPosition( 10, 560 );
  if (SHOW_GPU_TEXT){
    fs.vtk.textGPU->SetInput( "GPU" );
    fs.vtk.textGPU->GetTextProperty()->SetColor( 0.2,0.8,0.2 );
  }else{
    fs.vtk.textGPU->SetInput( "CPU" );
    fs.vtk.textGPU->GetTextProperty()->SetColor( 0.8,0.2,0.2 );
  }
  fs.vtk.ren->AddActor2D( fs.vtk.textGPU );

  // help text
  fs.vtk.help = vtkTextActor::New();
  if (fs.vtk.help == NULL){ exit_error("Error: invalid help textActor, could not create object\n"); }

  fs.vtk.help->GetTextProperty()->SetFontSize ( 16 );
  fs.vtk.help->SetPosition( 10, 80 );
  fs.vtk.help->GetTextProperty()->SetColor( 0.5, 0.5, 0.5 );
  fs.vtk.help->SetInput( (const char*) info_string.c_str() );
  fs.vtk.ren->AddActor2D( fs.vtk.help );

  // progress text
  fs.vtk.text = vtkTextActor::New();
  if (fs.vtk.text == NULL){ exit_error("Error: invalid progress textActor, could not create object\n"); }

  fs.vtk.text->GetTextProperty()->SetFontSize ( 16 );
  fs.vtk.text->SetPosition( 10, 530 );
  fs.vtk.text->SetInput( "...initializing data " );
  fs.vtk.text->GetTextProperty()->SetColor( 0.5,0.5,0.5 );
  fs.vtk.ren->AddActor2D( fs.vtk.text );

  // color table
  TRACE("     color table\n");
  int tableSize = 256;
  fs.vtk.icolor = 0; // from blue to red
  fs.vtk.colorAdjustOn = 1; // automatic adjust

  fs.vtk.lut = vtkLookupTable::New();
  if (fs.vtk.lut == NULL){ exit_error("Error: invalid LookupTable, could not create object\n"); }

  fs.vtk.lut->SetNumberOfColors(tableSize);
  // sets (default) rainbow color scale
  set_color_scale(fs.vtk.icolor);

  // render window
  TRACE("     render window\n");
  fs.vtk.renWin = vtkRenderWindow::New();
  if (fs.vtk.renWin == NULL){ exit_error("Error: invalid renderWindow, could not create object\n"); }

  fs.vtk.renWin->AddRenderer(fs.vtk.ren);
  fs.vtk.renWin->SetPosition(500,0);
  fs.vtk.renWin->SetSize(900,600);
  //fs.vtk.renWin->BordersOn();

  TRACE("init_vtk done successfully\n");
}

void init_vtk_window_interactor(){
  TRACE("init_vtk_window_interactor");

  // initializes vtk window interactor
  // note: this function will be called by separate vtk_thread to enable the interactor during the simulation
  //       otherwise it would halt execution of the main program

  // window interactor
  fs.vtk.iren = vtkRenderWindowInteractor::New();
  fs.vtk.iren->SetRenderWindow(fs.vtk.renWin);

  // by default a joystick style is used, this sets a trackball style (similar to paraview GUI)
  fs.vtk.style = vtkInteractorStyleTrackballCamera::New();
  fs.vtk.iren->SetInteractorStyle(fs.vtk.style);

  //printf("window has never been rendered: %d\n",fs.vtk.renWin->GetNeverRendered());

  // Initialize must be called prior to creating timer events.
  fs.vtk.iren->Initialize();

  // callback function
  fs.vtk.iren->AddObserver(vtkCommand::KeyPressEvent, &mykey, &MyInteractor::KeypressCallbackFunction);

  // below call might crash if window is rendered before interactor is initialized
  // error in X11: GLXBadCurrentWindow
  interactor_usage();

  // rendering will be called within this thread by a timer callback method
  CommandSubclass2* timerCallback = CommandSubclass2::New();
  fs.vtk.iren->AddObserver( vtkCommand::TimerEvent, timerCallback );
  // timer in milliseconds
  fs.vtk.iren->CreateRepeatingTimer(500);

  // indicate that window initialization is done
  pthread_mutex_lock(&vtk_thread_info.mutex);
  vtk_thread_info.init = false;
  pthread_mutex_unlock(&vtk_thread_info.mutex);

  // window features
  printf("\n     render window features:\n");
  if (fs.vtk.renWin->SupportsOpenGL()){
    printf("       supports OpenGL\n");
  }else{
    printf("       has no OpenGL support\n");
  }
  if (fs.vtk.renWin->IsDirect()){
    printf("       uses hardware acceleration\n");
  }else{
    printf("       uses no hardware acceleration\n");
  }
#if VTK_MAJOR_VERSION >= 7
  if (fs.vtk.renWin->IsDrawable()){
    printf("       is drawable\n");
  }else{
    printf("       is NOT drawable\n");
  }
#endif
  //printf("       window has never been rendered: %d\n",fs.vtk.renWin->GetNeverRendered());
  printf("\n");

  sleep(1);

  // starts interactor, listens to updates by callback and key events
  // note: this will not return until the interactor gets terminated...
  fs.vtk.iren->Start();

  TRACE("init_vtk_window_interactor done");
}


void clean_vtk_arrays(){
  TRACE("clean_vtk_arrays");

  // free timer
  fs.timer->StopTimer();
  fs.timer->Delete();

  // free arrays
  if(SHOW_FREESURFACE == 1 ){
    fs.vtk.points2D->Delete();
    fs.vtk.data_array2D->Delete();

    fs.vtk.lut2D->Delete();

    // contour
    fs.vtk.contourActor->Delete();
    fs.vtk.contourMapper->Delete();
    fs.vtk.contour->Delete();

    fs.vtk.actor2D->Delete();
    fs.vtk.mapMesh2D->Delete();
  }

  if(SHOW_VOLUMEDATA == 1 ){
    fs.vtk.points3D->Delete();
    fs.vtk.data_array3D->Delete();

    fs.vtk.legendcolor->Delete();
    fs.vtk.lut->Delete();

    if (use_vtk_extract){
      fs.vtk.extract->Delete();
    }else{
      fs.vtk.clipPlane1->Delete();
      fs.vtk.clipPlane2->Delete();
      fs.vtk.clip1->Delete();
      fs.vtk.clip1m->Delete();
      fs.vtk.clip2->Delete();
      fs.vtk.merger->Delete();
    }
    //fs.vtk.geometryFilter->Delete();
    //fs.vtk.cleanFilter->Delete();
    //fs.vtk.clippersurface->Delete();

    fs.vtk.actor3D->Delete();
    fs.vtk.mapMesh3D->Delete();

    free(global_data_array_peak);
  }

  // receivers
  fs.vtk.glyph->Delete();
  fs.vtk.glyphMapper->Delete();
  fs.vtk.glyphActor->Delete();

  // source sphere
  fs.vtk.actorSphere->Delete();
  fs.vtk.mapperSphere->Delete();
  fs.vtk.sphere->Delete();

  fs.vtk.textGPU->Delete();
  fs.vtk.text->Delete();
  fs.vtk.help->Delete();

  fs.vtk.iren->Delete();
  fs.vtk.renWin->Delete();
  fs.vtk.ren->Delete();
  fs.vtk.camera->Delete();
}


void sync_rendering(int sec){
  TRACE("sync_rendering");

  // waits for vtk rendering to finish
  if (sec < 20) sec = 20;

  // time out
  time_t current = time(NULL);
  time_t seconds = sec;
  time_t endtime = current + seconds;

  //printf("render loop time is    : %s", ctime(&current));
  //printf("render loop endtime is : %s", ctime(&endtime));

  while(fs.vtk.do_render != 0){
    sleep(1);
    // check if thread still up
    if (vtk_thread_info.finished == true) break;
    // in case the render flag never changes, time out
    current = time(NULL);
    if (current > endtime){
      printf("render loop time exceeded for vtk rendering, continuing...\n");
      break;
    }
  }
}


void sync_halt_simulation(){
  TRACE("sync_halt_simulation");

  // checks if simulation halted
  // note: this stalls if user selects another key first...!
  if (fs.vtk.haltOn == 1){
    TRACE("simulation halted");
    printf("simulation halted, waiting for pressing <space> again to continue...\n");
    while(fs.vtk.haltOn == 1) {
      sleep(1);
      // check if thread still up
      if (vtk_thread_info.finished == true) break;
    }
    printf("simulation continues...\n");
  }
}


/* ----------------------------------------------------------------------------------------------- */

// public external (callable within fortran) functions

extern "C"
void FC_FUNC_(initialize_vtkwindow,INITIALIZE_VTKWINDOW)(int* GPU_MODE) {
  TRACE("vtk window initialization");

  // flag determines if computation done on GPU or CPU
  SHOW_GPU_TEXT = *GPU_MODE;

  printf("vtk: initializing VTK window...\n");
  printf("     VTK version is %s\n\n",vtkVersion::GetVTKVersion());

#if VTK_MAJOR_VERSION >= 6
// unfortunately, on Apple OsX the main thread must run the GUI, otherwise it crashes like:
// ..
// xspecfem3D[3796:357270] +[NSUndoManager(NSInternal) _endTopLevelGroupings] is only safe to invoke on the main thread.
// ..
#if defined (__APPLE__) || defined(MACOSX)
  cout << "Sorry, running VTK visualization on Mac using newer VTK versions 6+ doesn't work properly yet...." << endl;
  exit(1);
#endif
#endif

#ifdef VTK_VIS_PARALLEL
  // just a test for now...
  // initializes vtk mpi
  vtkMPIController* controller = vtkMPIController::New();
  controller->Initialize();

  int rank = controller->GetLocalProcessId();
  int procs = controller->GetNumberOfProcesses();

  printf("vtk: mpi rank %d - total processes = %d\n\n",rank,procs);
  controller->Finalize();
#endif

  // initializes vtk
  init_vtk();

  // vtk window and interactor
  // starts interactor within separate thread, otherwise this will halt until interactor is terminated
  TRACE("vtk creating thread");
  create_vtk_thread();

  // wait for vtk creation to finish
  time_t current = time(NULL);
  time_t seconds = 10;
  time_t endtime = current + seconds;

  //printf("loop time is    : %s", ctime(&current));
  //printf("loop endtime is : %s", ctime(&endtime));

  while(vtk_thread_info.init == true){
    sleep(1);
    // in case the render flag never changes, time out
    current = time(NULL);
    //printf("loop time is    : %s", ctime(&current));
    //printf("loop endtime is : %s", ctime(&endtime));
    if (current > endtime){
      printf("loop time exceeded for initializing vtk, continuing...\n");
      break;
    }
  }
  printf("     initialization done\n\n");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float* xs_x,float* xs_y, float* xs_z) {
  TRACE("prepare vtk source glyph ...");

  // sets source location
  fs.vtk.pos_source[0] = *xs_x;
  fs.vtk.pos_source[1] = *xs_y;
  fs.vtk.pos_source[2] = *xs_z;

  // terminal output
  printf("vtk: source sphere\n");
  printf("     sphere location: x/y/z = %f / %f / %f \n",fs.vtk.pos_source[0],fs.vtk.pos_source[1],fs.vtk.pos_source[2]);

  // creates sphere around source location
  fs.vtk.sphere = vtkSphereSource::New();
  if (fs.vtk.sphere == NULL){ exit_error("Error: invalid sphere, could not create object\n"); }

  fs.vtk.sphere->SetCenter(fs.vtk.pos_source);
  fs.vtk.sphere->SetRadius(0.02);

  fs.vtk.mapperSphere = vtkPolyDataMapper::New();
  if (fs.vtk.mapperSphere == NULL){ exit_error("Error: invalid sphere mapper, could not create object\n"); }

  fs.vtk.mapperSphere->SetInputConnection(fs.vtk.sphere->GetOutputPort());

  fs.vtk.actorSphere = vtkActor::New();
  if (fs.vtk.actorSphere == NULL){ exit_error("Error: invalid sphere actor, could not create object\n"); }

  fs.vtk.actorSphere->SetMapper(fs.vtk.mapperSphere);

  // adds actor
  fs.vtk.ren->AddActor(fs.vtk.actorSphere);

  // render window (with camera reset)
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 4;
  fs.vtk.mutex->Unlock();

  // waits until rendered
  sleep(5);
  sync_rendering();
  printf("     done\n\n");
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int* free_np,
                                                             float* free_x,
                                                             float* free_y,
                                                             float* free_z,
                                                             int* free_nspec,
                                                             int* free_conn) {
  TRACE("prepare vtk freesurface...");

  float xyz[3];
  int id1,id2,id3,id4;

  // initializes
  SHOW_FREESURFACE = 1;

  // number of points
  fs.freesurface.np = *free_np;

  // terminal output
  printf("vtk: free surface\n");
  printf("     free surface points    : %d\n", fs.freesurface.np);
  printf("     free surface array size: %f (MB)\n", fs.freesurface.np * 4 / 1024. / 1024.);

  // checks
  if( fs.freesurface.np == 0 ){
    fprintf(stderr,"ERROR: VTK window without 2D freesurface points \n");
    exit(1);
  }

  // creates new points and data arrays
  fs.vtk.points2D = vtkPoints::New();
  fs.vtk.points2D->SetNumberOfPoints(fs.freesurface.np);

  fs.vtk.data_array2D = vtkFloatArray::New();
  fs.vtk.data_array2D->SetNumberOfComponents(1);
  fs.vtk.data_array2D->SetNumberOfValues(fs.freesurface.np);
  fs.vtk.data_array2D->SetName("topo");

  for(int i=0;i<fs.freesurface.np;i++) {
    xyz[0] = free_x[i];
    xyz[1] = free_y[i];
    xyz[2] = free_z[i];
    fs.vtk.points2D->SetPoint(i,xyz);
    // assigns vertical z-value
    fs.vtk.data_array2D->SetValue(i,1.0*xyz[2]);
  }
  // Find min and max
  float bounds[2];
  fs.vtk.data_array2D->GetValueRange(bounds);
  double min = bounds[0];
  double max = bounds[1];
  printf("     topo: min = %f max = %f\n", min,max);

  // Find min and max
  double model_bounds[6];
  fs.vtk.points2D->GetBounds(model_bounds);

  // updates freesurface bounds
  // note: will be used for size of receiver glyphs
  for(int i=0;i<6;i++) fs.freesurface.bounds[i] = model_bounds[i];

  printf("     topo size: x min/max = %f/%f   y min/max = %f/%f   z min/max = %f/%f \n",
         fs.freesurface.bounds[0],fs.freesurface.bounds[1],
         fs.freesurface.bounds[2],fs.freesurface.bounds[3],
         fs.freesurface.bounds[4],fs.freesurface.bounds[5]);
  printf("     topo size: x-range = %f   y-range = %f   z-range = %f \n",
         fs.freesurface.bounds[1] - fs.freesurface.bounds[0],
         fs.freesurface.bounds[3] - fs.freesurface.bounds[2],
         fs.freesurface.bounds[5] - fs.freesurface.bounds[4]);

  // black color scale
  int tableSize = 256;
  fs.vtk.lut2D = vtkLookupTable::New();
  fs.vtk.lut2D->SetNumberOfColors(tableSize);
  fs.vtk.lut2D->SetValueRange( 0.1, 1.0 ); // from black to white
  fs.vtk.lut2D->SetHueRange( 0.0, 1.0 );
  fs.vtk.lut2D->SetSaturationRange( 0.0, 0.0 ); // no color saturation
  fs.vtk.lut2D->SetTableRange( min, max );
  fs.vtk.lut2D->Build();

  // creates cell connectivity
  fs.freesurface.nspec = *free_nspec;

  // cells
  vtkCellArray* cells2D = vtkCellArray::New();
  vtkQuad* quad = vtkQuad::New();

  for(int ispec=0;ispec<fs.freesurface.nspec;ispec++){
    id1 = free_conn[0+ispec*4];
    id2 = free_conn[1+ispec*4];
    id3 = free_conn[2+ispec*4];
    id4 = free_conn[3+ispec*4];
    quad->GetPointIds()->SetId(0,id1);
    quad->GetPointIds()->SetId(1,id2);
    quad->GetPointIds()->SetId(2,id3);
    quad->GetPointIds()->SetId(3,id4);
    cells2D->InsertNextCell(quad);
  }
  quad->Delete();

  fs.vtk.volume2D = vtkUnstructuredGrid::New();
  // points
  fs.vtk.volume2D->SetPoints(fs.vtk.points2D);
  fs.vtk.volume2D->GetPointData()->SetScalars(fs.vtk.data_array2D);
  // cells
  fs.vtk.volume2D->SetCells(VTK_QUAD, cells2D);

  // frees array
  cells2D->Delete();

  // contour iso-surfacing
  if (CONTOUR_FREESURFACE){
    fs.vtk.contour = vtkContourFilter::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.contour->SetInput( fs.vtk.volume2D );
#else
    fs.vtk.contour->SetInputData( fs.vtk.volume2D );
#endif
    fs.vtk.contour->SetNumberOfContours(7);
    if (min < 0.0 && max > 0.0){
      fs.vtk.contour->SetValue(0, 0.0);
      fs.vtk.contour->SetValue(1, min*0.25);
      fs.vtk.contour->SetValue(2, max*0.25);
      fs.vtk.contour->SetValue(3, min*0.5);
      fs.vtk.contour->SetValue(4, max*0.5);
      fs.vtk.contour->SetValue(5, min*0.75);
      fs.vtk.contour->SetValue(6, max*0.75);
    }else{
      fs.vtk.contour->SetValue(0, min + 0.5*(max-min)); // middle
      fs.vtk.contour->SetValue(1, min + 0.125*(max-min));
      fs.vtk.contour->SetValue(2, min + 0.25*(max-min));
      fs.vtk.contour->SetValue(3, min + 0.375*(max-min));
      fs.vtk.contour->SetValue(4, min + 0.625*(max-min));
      fs.vtk.contour->SetValue(5, min + 0.75*(max-min));
      fs.vtk.contour->SetValue(6, min + 0.875*(max-min));
    }
    fs.vtk.contour->ComputeScalarsOff();
    fs.vtk.contour->ComputeGradientsOff();

    fs.vtk.contourMapper = vtkPolyDataMapper::New();
    fs.vtk.contourMapper->SetInputConnection(0, fs.vtk.contour->GetOutputPort() );
    fs.vtk.contourMapper->ScalarVisibilityOff();

    fs.vtk.contourActor = vtkActor::New();
    fs.vtk.contourActor->SetMapper( fs.vtk.contourMapper );

    // properties
    fs.vtk.contourActor->GetProperty()->SetOpacity( 1.0 );
    fs.vtk.contourActor->SetScale( 1.0, 1.0, vertical_exageration );
    fs.vtk.contourActor->GetProperty()->SetColor( 0.5, 0.5, 0.5 );

    fs.vtk.ren->AddActor(fs.vtk.contourActor);
  }

  // mapping
  fs.vtk.mapMesh2D = vtkDataSetMapper::New();
#if VTK_MAJOR_VERSION <= 5
  fs.vtk.mapMesh2D->SetInput( fs.vtk.volume2D );
#else
  fs.vtk.mapMesh2D->SetInputData( fs.vtk.volume2D );
#endif

  fs.vtk.mapMesh2D->SetLookupTable(fs.vtk.lut2D);
  fs.vtk.mapMesh2D->SetColorModeToMapScalars();
  fs.vtk.mapMesh2D->SetScalarModeToUsePointData();
  fs.vtk.mapMesh2D->ScalarVisibilityOn();
  fs.vtk.mapMesh2D->UseLookupTableScalarRangeOn();

  // actor
  fs.vtk.actor2D = vtkActor::New();
  fs.vtk.actor2D->SetMapper(fs.vtk.mapMesh2D);

  // properties
  // opacity
  fs.vtk.actor2D->GetProperty()->SetOpacity( 0.5 );
  // vertical exageration
  fs.vtk.actor2D->SetScale( 1.0, 1.0, vertical_exageration );

  // 2D actor
  fs.vtk.ren->AddActor(fs.vtk.actor2D);

  // render window (with camera reset)
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 4;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();
  printf("     done\n\n");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
// character array doesn't work properly...
//void FC_FUNC_(prepare_vtkreceivers,PREPARE_VTKRECEIVERS)(int* nrec_h, float* xr_x,float* xr_y, float* xr_z, char sta_name[][MAX_LENGTH_STATION_NAME]) {
// or single, long string
void FC_FUNC_(prepare_vtkreceivers,PREPARE_VTKRECEIVERS)(int* nrec_h, float* xr_x,float* xr_y, float* xr_z, int* max_len_station, char* sta_names) {
  TRACE("prepare vtk receivers glyph ...");

  int nrec = *nrec_h;

  // checks if anything to do
  if (nrec <= 0) return;

  // terminal output
  printf("vtk: receivers glyphs\n");
  printf("     number of receivers: %d \n",nrec);

  // sets receiver locations
  vtkPoints* points = vtkPoints::New();
  for (int i=0; i<nrec; i++){
    points->InsertNextPoint(xr_x[i],xr_y[i],xr_z[i]);
  }
  vtkPolyData* polydata = vtkPolyData::New();
  polydata->SetPoints(points);

  // creates sphere glyph around location
  vtkSphereSource* sphere = vtkSphereSource::New();
  // determines sphere radius depending on x/y range (updated in setting freesurface)
  double radius = 100.0;
  double dist_min = 1.e31;
  dist_min = MIN(dist_min,fabs(fs.freesurface.bounds[1]-fs.freesurface.bounds[0]));
  dist_min = MIN(dist_min,fabs(fs.freesurface.bounds[3]-fs.freesurface.bounds[2]));
  // estimates cell size (lenght / number of cells )
  if ( dist_min > 0.0){
    radius = dist_min / sqrt(fs.freesurface.nspec);
  }
  printf("     receiver glyph: estimated element size = %f \n",radius);

  // determines largest cell dimension to adapt sphere radius
  int id_cell = fs.vtk.volume2D->GetMaxCellSize();
  double bounds[6];
  fs.vtk.volume2D->GetCellBounds(id_cell,bounds);
  //printf("     receiver glyph: maxcellsize %d type %d\n\n",id_cell,fs.vtk.volume2D->GetCellType(id_cell));
  //printf("     receiver glyph: bounds %f %f %f %f %f %f\n\n",bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5]);
  double cell_max = 0.0;
  cell_max = MAX(cell_max,fabs(bounds[1]-bounds[0])); // dx
  cell_max = MAX(cell_max,fabs(bounds[3]-bounds[2])); // dy
  printf("     receiver glyph: maximum element size = %f \n",cell_max);
  if (cell_max > 0.0){
    radius = cell_max / 2.0;
  }
  printf("     receiver glyph: setting radius size = %f \n\n",radius);

  // sets glyph shape
  sphere->SetRadius(radius);
  sphere->SetPhiResolution(4);
  sphere->SetThetaResolution(4);

  // note: using 2D glyph is not faster than using 3D glyphs...
  // 2D glyph
  //vtkGlyphSource2D* symbol = vtkGlyphSource2D::New();
  //symbol->SetGlyphType(1); // 5 = triangle, 7 = circle
  //symbol->FilledOff();
  //fs.vtk.glyph = vtkGlyph2D::New();
  //fs.vtk.glyph->SetSource(symbol->GetOutput());

  // 3D glyph
  fs.vtk.glyph = vtkGlyph3D::New();
#if VTK_MAJOR_VERSION <= 5
  fs.vtk.glyph->SetSource(sphere->GetOutput());
  fs.vtk.glyph->SetInput(polydata);
#else
  fs.vtk.glyph->SetSourceData(sphere->GetOutput());
  fs.vtk.glyph->SetInputData(polydata);
#endif
  //fs.vtk.glyph->SetScaleModeToDataScalingOff();
  fs.vtk.glyph->OrientOff();
  fs.vtk.glyph->Update();

  // mapper
  fs.vtk.glyphMapper = vtkPolyDataMapper::New();
  fs.vtk.glyphMapper->SetInputConnection(fs.vtk.glyph->GetOutputPort());

  // adds actor
  fs.vtk.glyphActor = vtkActor::New();
  fs.vtk.glyphActor->SetMapper(fs.vtk.glyphMapper);
  fs.vtk.glyphActor->VisibilityOn();
  fs.vtk.glyphActor->GetProperty()->SetColor(0.5, 0.5, 0.5 );

  // adds actor
  fs.vtk.ren->AddActor(fs.vtk.glyphActor);

  // receiver labels
  int label_max, stride, icount;
  // maximum number of displayed labels
  if (use_reduced_labeling){
    // limits number of labels
    label_max = max_receiver_labels; // for a limited amount of labels
    stride = ceil( (float) nrec / (float) label_max );
  }else{
    // no need for an additional array
    label_max = nrec;  // for example: = 20 for a limited display
    stride = 1;
  }

  // sets receiver labels locations
  vtkPoints* pointsLabels = vtkPoints::New();
  icount = 0;
  for (int i=0; i<nrec; i++){
    if (i % stride == 0){
      icount++;
      pointsLabels->InsertNextPoint(xr_x[i],xr_y[i],xr_z[i]);
    }
  }
  printf("     stations total: %d - stride: %d display: %d\n",nrec,stride,icount);
  //printf("     stations strings: %s - length: %d\n",sta_names, (int) len(sta_names));

  vtkPolyData* polydataLabels = vtkPolyData::New();
  polydataLabels->SetPoints(pointsLabels);

  // labeling array
  vtkStringArray* labels = vtkStringArray::New();
  labels->SetNumberOfValues(icount);
  labels->SetName("labels");

  char* buffer = sta_names;
  int max_l = *max_len_station;
  char label[max_l];
  int istart;

  icount = 0;
  for (int i=0; i<nrec; i++){
    if (i % stride == 0){
      // moves pointer to begining of new station name in string
      istart = i * max_l;
      buffer = sta_names + istart; // moves char pointer forward
      // copies length - 1 to be able to add null termination
      strncpy(label,buffer,max_l - 1);
      // adds null character
      label[max_l - 1] = '\0';
      // trims
      char* blank = strchr(label, ' ');
      if (blank != NULL) {
        label[blank - label] = '\0';
      }
      printf("     label station %d: %s\n",i,label);
      // adds label
      labels->SetValue(icount, label);
      icount++;
    }
  }
  printf("     number of labeled stations = %d \n",icount);

  // adds labels to data
  polydataLabels->GetPointData()->AddArray(labels);

  // only visible points should be labelled
  // note: this visibility filter doesn't work correctly,
  //       somehow only points behind contour lines will disappear...? not sure why
  printf("     adding visibility filter\n");
  fs.vtk.labelVis = vtkSelectVisiblePoints::New();
#if VTK_MAJOR_VERSION <= 5
  fs.vtk.labelVis->SetInput( polydataLabels );
#else
  fs.vtk.labelVis->SetInputData( polydataLabels );
#endif
  fs.vtk.labelVis->SetRenderer( fs.vtk.ren );
  // note: if X window looses focus (e.g. by clicking on terminal window), then the update below will crash...
  //fs.vtk.labelVis->Update();

  // sets text color to grey
  //printf("     mapping labels \n");
  fs.vtk.labelMapper = vtkLabeledDataMapper::New();
  fs.vtk.labelMapper->SetInputConnection( fs.vtk.labelVis->GetOutputPort() );
  fs.vtk.labelMapper->SetLabelModeToLabelFieldData();
  fs.vtk.labelMapper->SetFieldDataName("labels");
  fs.vtk.labelMapper->GetLabelTextProperty()->SetColor( 0.8, 0.8, 0.8 );
  fs.vtk.labelMapper->GetLabelTextProperty()->SetFontFamilyToArial();
  fs.vtk.labelMapper->GetLabelTextProperty()->BoldOff();
  fs.vtk.labelMapper->GetLabelTextProperty()->ItalicOff();
  fs.vtk.labelMapper->GetLabelTextProperty()->ShadowOff();

  // note: this can crash on OsX in libGL when accessing the Z-buffer (maybe only defined within other thread)
  //fs.vtk.labelMapper->Update();

  // adds labels as separate actor
  printf("     hiding labels (by default)\n");
  fs.vtk.labelActor = vtkActor2D::New();
  fs.vtk.labelActor->SetMapper( fs.vtk.labelMapper );
  fs.vtk.labelActor->VisibilityOff();

  // adds actor
  fs.vtk.ren->AddActor(fs.vtk.labelActor);

  // render window (with camera reset)
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 4;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();

  // cleanup: not sure if vtk is doing deep copies or not...
  points->Delete();
  sphere->Delete();
  polydata->Delete();
  labels->Delete();
  pointsLabels->Delete();
  polydataLabels->Delete();

  printf("     done\n\n");
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int* vol_np,
                                                 float* vol_x,
                                                 float* vol_y,
                                                 float* vol_z,
                                                 int* vol_nspec,
                                                 int* vol_conn) {
  TRACE("prepare_vtkfield");

  float xyz[3];
  float data_bounds[2];
  double model_bounds[6];

  // debug vtk file output
  static int debug_file = 0;

  // initializes
  SHOW_VOLUMEDATA = 1;

  // volumetric wavefield
  // window text
  fs.vtk.text->SetInput( "...adding wavefield " );

  // update view
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 1;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();

  // volume mesh
  fs.volume.np = *vol_np;
  fs.volume.nspec = *vol_nspec;

  // terminal output
  printf("vtk: volume\n");
  printf("     volume points    : %d\n", fs.volume.np);
  printf("     volume elements  : %d\n", fs.volume.nspec);
  printf("     volume array size: %f (MB)\n", fs.volume.np * 4 / 1024. / 1024.);

  // checks
  if( fs.volume.np == 0 ){
    fprintf(stderr,"ERROR: VTK window without 3D volume points \n");
    exit(1);
  }

  // peak value array
  global_data_array_peak = (float *) malloc( fs.volume.np * sizeof(float));
  if (global_data_array_peak == NULL){
    printf("vtk Error: failed to allocate peak value buffer\n");
    abort();
  }

  // volume data
  // point locations
  fs.vtk.points3D = vtkPoints::New();
  fs.vtk.points3D->SetNumberOfPoints(fs.volume.np);

  // point data (wavefield)
  fs.vtk.data_array3D = vtkFloatArray::New();
  fs.vtk.data_array3D->SetNumberOfComponents(1);
  fs.vtk.data_array3D->SetName("vnorm");
  fs.vtk.data_array3D->SetNumberOfValues(fs.volume.np);

  for(int i=0;i<fs.volume.np;i++) {
    xyz[0] = vol_x[i];
    xyz[1] = vol_y[i];
    xyz[2] = vol_z[i];
    fs.vtk.points3D->SetPoint(i,xyz);
    fs.vtk.data_array3D->SetValue(i,0.0);
    //fs.vtk.data_array3D->SetValue(i,i*1.0);
    // initializes peak values
    global_data_array_peak[i] = 0.0f;
  }
  //fs.vtk.data_array3D->Update();
  // Find min and max
  fs.vtk.data_array3D->GetValueRange(data_bounds);
  double dmin = data_bounds[0];
  double dmax = data_bounds[1];
  printf("     volume data: min = %f max = %f\n\n",dmin,dmax);

  // Find min and max
  fs.vtk.points3D->GetBounds(model_bounds);

  // updates volume bounds
  for(int i=0;i<6;i++) fs.volume.bounds[i] = model_bounds[i];

  double xmin = model_bounds[0];
  double xmax = model_bounds[1];
  double ymin = model_bounds[2];
  double ymax = model_bounds[3];
  double zmin = model_bounds[4];
  double zmax = model_bounds[5];
  printf("     model dimension: xmin/xmax = %f / %f\n",xmin,xmax);
  printf("                      ymin/ymax = %f / %f\n",ymin,ymax);
  printf("                      zmin/zmax = %f / %f\n\n",zmin,zmax);

  // determines approximate minimum point distance to set a merge tolerance
  // has not effect...?
  double dist_min = 1.e31; // something huge
  dist_min = MIN(dist_min,xmax - xmin);
  dist_min = MIN(dist_min,ymax - ymin);
  dist_min = MIN(dist_min,zmax - zmin);
  dist_min = dist_min / pow( fs.volume.nspec,(1./3.));
  printf("     estimated mesh point minimum distance %f\n",dist_min);

  // adjusts camera settings
  fs.vtk.rclip[0] = xmin - 0.1*(xmax-xmin);
  fs.vtk.rclip[1] = xmax + 0.1*(xmax-xmin);
  fs.vtk.pfocal[0] = fs.vtk.pos_source[0]; //(xmax-xmin)/2.0;
  fs.vtk.pfocal[1] = fs.vtk.pos_source[1]; //(ymax-ymin)/2.0;
  fs.vtk.pfocal[2] = fs.vtk.pos_source[2]; // (zmax-zmin)/2.0;
  fs.vtk.pcam[0] = fs.vtk.pos_source[0] + 0.5*(zmax-zmin); // (xmax-xmin)/2.0 - 0.1*(xmax-xmin);
  fs.vtk.pcam[1] = fs.vtk.pos_source[1] + 0.5*(zmax-zmin); // (ymax-ymin)/2.0 + 0.1*(ymax-ymin);
  fs.vtk.pcam[2] = fs.vtk.pos_source[2] + 0.5*(zmax-zmin); // (zmax-zmin)/2.0 + 0.5*(zmax-zmin);

  // range
  fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
  // camer focal point
  fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
  // camer position
  fs.vtk.camera->SetPosition( fs.vtk.pcam );

  // reposition the camera, so that actor can be fully seen
  fs.vtk.ren->ResetCamera();
  fs.vtk.ren->ResetCameraClippingRange();

  // adjust sphere size
  fs.vtk.sphere->SetRadius(0.02*fabs(zmax-zmin));
  fs.vtk.sphere->Modified();

  //
  // unstructured grid
  //
  // creates cell connectivity
  // cells
  vtkCellArray* cells = vtkCellArray::New();
  vtkHexahedron* hex = vtkHexahedron::New();

  int id;
  double dx1,dx2,dy1,dy2;

  double dist_min_x = 1.e31;
  double dist_min_y = 1.e31;

  for (int ispec=0;ispec<fs.volume.nspec;ispec++){
    // ordering of points can change depending on mesh...
    // for example
    //  (1,1,1),(1,1,NGLLZ),(1,NGLLY,NGLLZ),(1,NGLLY,1),(NGLLX,1,1),(NGLLX,1,NGLLZ),(NGLLX,NGLLY,NGLLZ),(NGLLX,NGLLY,1)
    // or
    //  (1,1,1),(NGLLX,1,1),(NGLLX,NGLLY,1),(1,NGLLY,1),(1,1,NGLLX),(NGLLX,1,NGLLZ),(NGLLX,NGLLY,NGLLZ),(1,NGLLY,NGLLZ)
    //printf("cell: %d\n",ispec);

    // sets cell ids
    for (int j=0; j<8; j++){
      id = vol_conn[j+ispec*8];
      if (id < 0 || id >= fs.volume.np) {exit_error("Error: cell id exceeds array bounds");}
      hex->GetPointIds()->SetId(j,id);
      // debug
      //printf("cell corner %d: id %d - x/y/z = %f / %f / %f \n",j,id,vol_x[id],vol_y[id],vol_z[id]);
    }
    // adds cell
    cells->InsertNextCell(hex);

    // estimates minimum cell size along x/y
    // distances along x-direction: for points (NGLLX,1,1,ispec) - (1,1,1,ispec)
    dx1 = fabs(vol_x[vol_conn[4+ispec*8]] - vol_x[vol_conn[0+ispec*8]]);
    dx2 = fabs(vol_x[vol_conn[1+ispec*8]] - vol_x[vol_conn[0+ispec*8]]);
    if (dx1 > dx2){
      dist_min_x = MIN(dist_min_x,dx1);
    }else{
      dist_min_x = MIN(dist_min_x,dx2);
    }

    // distances along y-direction: for points (1,NGLLY,1,ispec) - (1,1,1,ispec)
    dy1 = fabs(vol_y[vol_conn[3+ispec*8]] - vol_y[vol_conn[0+ispec*8]]);
    dy2 = fabs(vol_y[vol_conn[2+ispec*8]] - vol_y[vol_conn[0+ispec*8]]);
    if (dy1 > dy2){
      dist_min_y = MIN(dist_min_y,dy1);
    }else{
      dist_min_y = MIN(dist_min_y,dy2);
    }
  }
  hex->Delete();

  printf("     estimated cell size minimum distance along x = %f\n",dist_min_x);
  printf("     estimated cell size minimum distance along y = %f\n\n",dist_min_y);

  // uses same minimum for both directions
  dist_min_x = MAX(dist_min_x,dist_min_y);
  dist_min_y = MAX(dist_min_x,dist_min_y);

  fs.vtk.volume3D = vtkUnstructuredGrid::New();
  // points
  fs.vtk.volume3D->SetPoints(fs.vtk.points3D);
  fs.vtk.volume3D->GetPointData()->SetScalars(fs.vtk.data_array3D);
  // cells
  fs.vtk.volume3D->SetCells(VTK_HEXAHEDRON, cells);

  // frees array
  cells->Delete();

  // contour iso-surfacing
  if (CONTOUR_VOLUME){
    dmin = gcolor_min;
    dmax = gcolor_max;
    fs.vtk.contour3D = vtkContourFilter::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.contour3D->SetInput( fs.vtk.volume3D );
#else
    fs.vtk.contour3D->SetInputData( fs.vtk.volume3D );
#endif

    fs.vtk.contour3D->SetNumberOfContours(4);
    fs.vtk.contour3D->SetValue(0, dmin + (dmax-dmin)*0.2);
    fs.vtk.contour3D->SetValue(1, dmin + (dmax-dmin)*0.4);
    fs.vtk.contour3D->SetValue(2, dmin + (dmax-dmin)*0.6);
    fs.vtk.contour3D->SetValue(3, dmin + (dmax-dmin)*0.8);

    fs.vtk.contour3D->ComputeScalarsOn();
    fs.vtk.contour3D->ComputeGradientsOff();

    fs.vtk.contour3DMapper = vtkPolyDataMapper::New();
    fs.vtk.contour3DMapper->SetInputConnection(0, fs.vtk.contour3D->GetOutputPort() );
    fs.vtk.contour3DMapper->ScalarVisibilityOn();

    fs.vtk.contour3DMapper->SetLookupTable(fs.vtk.lut);
    fs.vtk.contour3DMapper->SetColorModeToMapScalars();
    fs.vtk.contour3DMapper->SetScalarModeToUsePointData();
    fs.vtk.contour3DMapper->ScalarVisibilityOn();
    fs.vtk.contour3DMapper->UseLookupTableScalarRangeOn();

    fs.vtk.contour3DActor = vtkActor::New();
    fs.vtk.contour3DActor->SetMapper( fs.vtk.contour3DMapper );
    fs.vtk.contour3DActor->GetProperty()->SetOpacity( 0.3 );

    //fs.vtk.contour3DActor->SetScale( 1.0,1.0,1.02 );
    //fs.vtk.contour3DActor->GetProperty()->SetColor( 0.5, 0.5, 0.5 );

    fs.vtk.ren->AddActor(fs.vtk.contour3DActor);
  }

  // idea: let us cut out a portion of the mesh to better look inside the volume, to see the wavefield
  //       the cut portion is determined where the source location is.

  //
  // clip box
  //

  // vector to source
  double v[3];
  v[0] = fs.vtk.pos_source[0];
  v[1] = fs.vtk.pos_source[1];
  v[2] = fs.vtk.pos_source[2];

  if (use_vtk_extract){
    // box dimensions
    if( fabs(v[0]-model_bounds[0]) < fabs(v[0]-model_bounds[1]) ){
      xmin = model_bounds[0];
      xmax = v[0];
    }else{
      xmin = v[0];
      xmax = model_bounds[1];
    }
    if( fabs(v[1]-model_bounds[2]) < fabs(v[1]-model_bounds[3]) ){
      ymin = model_bounds[2];
      ymax = v[1];
    }else{
      ymin = v[1];
      ymax = model_bounds[3];
    }
    zmin = model_bounds[4];
    zmax = model_bounds[5];

    // clip box
    vtkBox* clipbox1 = vtkBox::New();
    clipbox1->SetBounds(xmin,xmax,ymin,ymax,zmin,zmax);

    // extracting cells
    fs.vtk.extract = vtkExtractGeometry::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.extract->SetInput( fs.vtk.volume3D );
#else
    fs.vtk.extract->SetInputData( fs.vtk.volume3D );
#endif

    fs.vtk.extract->SetImplicitFunction( clipbox1 );
    // extracts cells outside of box
    fs.vtk.extract->SetExtractInside(0);
    fs.vtk.extract->Update();
    printf("     extract: number of points = %d\n",(int)fs.vtk.extract->GetOutput()->GetNumberOfPoints());
    printf("     extract: number of cells  = %d\n\n",(int)fs.vtk.extract->GetOutput()->GetNumberOfCells());

    // free memory
    clipbox1->Delete();
  }else{
    // for table based clipping ..
    // note: we will use the extract filter for now, since the table clipping has problems
    //       when clip functions align with cell boundaries

    // to define cut planes, we need a position and a normal vector
    // source location vector / perpendicular vectors
    //
    // vectorproduct(vector1, vector2, product)
    // calculates the vector product of vector1 and vector2
    //    product(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
    //    product(2) = - vector1(1)*vector2(3) + vector1(3)*vector2(1)
    //    product(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)

    // normals
    double pn1[3],pn2[3];

    // general
    //double p1[3],p2[3];
    //p1[0] = 0.0;
    //p1[1] = 0.0;
    //p1[2] = 1.0;
    //p2[0] = 0.0;
    //p2[1] = 1.0;
    //p2[2] = 0.0;
    //pn1[0] = v[1]*p1[2] - v[2]*p1[1];
    //pn1[1] = -v[0]*p1[2] + v[2]*p1[0];
    //pn1[2] = v[0]*p1[1] - v[1]*p1[0];
    //pn2[0] = v[1]*p2[2] - v[2]*p2[1];
    //pn2[1] = -v[0]*p2[2] + v[2]*p2[0];
    //pn2[2] = v[0]*p2[1] - v[1]*p2[0];

    // normal vector aligned with x/y-direction
    pn1[0] = -1.0;
    pn1[1] = 0.0;
    pn1[2] = 0.0;

    pn2[0] = 0.0;
    pn2[1] = -1.0;
    pn2[2] = 0.0;

    // flips normal depending on location of source
    // (to cut out only smaller portion)
    double sign;
    // x-direction
    xmin = model_bounds[0];
    xmax = model_bounds[1];
    if( fabs(v[0] - xmin) < fabs(v[0] - xmax) ){
      sign = -1.0;
    }else{
      sign = 1.0;
    }
    for (int i=0; i<3; i++){
      pn1[i] *= sign;
    }
    // adds x - offset to cut plane
    v[0] -= sign * dist_min_x * 2.0;

    // y-direction
    ymin = model_bounds[2];
    ymax = model_bounds[3];
    if( fabs(v[1] - ymin) < fabs(v[1] - ymax) ){
      sign = -1.0;
    }else{
      sign = 1.0;
    }
    for (int i=0; i<3; i++){
      pn2[i] *= sign;
    }
    // adds y - offset to cut plane
    v[1] -= sign * dist_min_y * 2.0;

    // plane through source location
    fs.vtk.clipPlane1 = vtkPlane::New();
    fs.vtk.clipPlane1->SetNormal( pn1 );
    fs.vtk.clipPlane1->SetOrigin( v );

    fs.vtk.clipPlane2 = vtkPlane::New();
    fs.vtk.clipPlane2->SetNormal( pn2 );
    fs.vtk.clipPlane2->SetOrigin( v );

    // table based clipping
    // note: table base clipping is faster than using vtkClipDataSet.
    //       however, it might miss surfaces when exactly aligned with geometry
    //       the clipdataset avoids some of these problems...
    //fs.vtk.clip1 = vtkClipDataSet::New();
    fs.vtk.clip1 = vtkTableBasedClipDataSet::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.clip1->SetInput( fs.vtk.volume3D );
#else
    fs.vtk.clip1->SetInputData( fs.vtk.volume3D );
#endif

    fs.vtk.clip1->SetClipFunction( fs.vtk.clipPlane1 );
    fs.vtk.clip1->Update();
    printf("     clip1: number of points = %d\n",(int)fs.vtk.clip1->GetOutput()->GetNumberOfPoints());
    printf("     clip1: number of cells  = %d\n",(int)fs.vtk.clip1->GetOutput()->GetNumberOfCells());
    printf("     clip1: merge tolerance  = %f\n\n",fs.vtk.clip1->GetMergeTolerance());

    //fs.vtk.clip1m = vtkClipDataSet::New();
    fs.vtk.clip1m = vtkTableBasedClipDataSet::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.clip1m->SetInput( fs.vtk.volume3D );
#else
    fs.vtk.clip1m->SetInputData( fs.vtk.volume3D );
#endif

    fs.vtk.clip1m->SetClipFunction( fs.vtk.clipPlane1 );
    fs.vtk.clip1m->InsideOutOn();
    fs.vtk.clip1m->Update();

    //fs.vtk.clip2 = vtkClipDataSet::New();
    fs.vtk.clip2 = vtkTableBasedClipDataSet::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.clip2->SetInput( fs.vtk.clip1m->GetOutput() );
#else
    fs.vtk.clip2->SetInputData( fs.vtk.clip1m->GetOutput() );
#endif
    fs.vtk.clip2->SetClipFunction( fs.vtk.clipPlane2 );
    fs.vtk.clip2->Update();

    // merges mesh parts
    fs.vtk.merger = vtkAppendFilter::New();
#if VTK_MAJOR_VERSION <= 5
    fs.vtk.merger->AddInput( fs.vtk.clip1->GetOutput() );
    fs.vtk.merger->AddInput( fs.vtk.clip2->GetOutput() );
#else
    fs.vtk.merger->AddInputData( fs.vtk.clip1->GetOutput() );
    fs.vtk.merger->AddInputData( fs.vtk.clip2->GetOutput() );
#endif
    fs.vtk.merger->MergePointsOn();

    /*
    // instead of two clips, this would create 3 blocks for merging...
    // bounds
    double bb_all[6] = { model_bounds[0],model_bounds[1],model_bounds[2],model_bounds[3],model_bounds[4],model_bounds[5] };
    double dist_x = bb_all[1] - bb_all[0];
    double dist_y = bb_all[3] - bb_all[2];
    double dist_z = bb_all[5] - bb_all[4];

    // 1. mesh part
    // unstructured grid clipper
    double bb1[6] = { bb_all[0], bb_all[1], bb_all[2], bb_all[3], bb_all[4], bb_all[4]+ 0.5*dist_z};
    fs.vtk.clipper1 = vtkBoxClipDataSet::New();
    fs.vtk.clipper1->SetInput( fs.vtk.volume3D );
    fs.vtk.clipper1->SetOrientation( 0 );
    fs.vtk.clipper1->SetBoxClip( bb1[0],bb1[1],bb1[2],bb1[3],bb1[4],bb1[5] );
    fs.vtk.clipper1->Update();

    // 2. mesh part
    double bb2[6] = { bb_all[0]+0.6*dist_x, bb_all[1], bb_all[2], bb_all[3], bb_all[4]+0.5*dist_z, bb_all[5]};
    fs.vtk.clipper2 = vtkBoxClipDataSet::New();
    fs.vtk.clipper2->SetInput( fs.vtk.volume3D );
    fs.vtk.clipper2->SetOrientation( 0 );
    fs.vtk.clipper2->SetBoxClip( bb2[0],bb2[1],bb2[2],bb2[3],bb2[4],bb2[5] );
    fs.vtk.clipper2->Update();

    // 3. mesh part
    double bb3[6] = { bb_all[0], bb_all[0]+0.6*dist_x, bb_all[2], bb_all[2]+0.5*dist_y, bb_all[4]+0.5*dist_z, bb_all[5]};
    fs.vtk.clipper3 = vtkBoxClipDataSet::New();
    fs.vtk.clipper3->SetInput( fs.vtk.volume3D );
    fs.vtk.clipper3->SetOrientation( 0 );
    fs.vtk.clipper3->SetBoxClip( bb3[0],bb3[1],bb3[2],bb3[3],bb3[4],bb3[5] );
    fs.vtk.clipper3->Update();

    // merges mesh parts
    fs.vtk.merger = vtkAppendFilter::New();
    fs.vtk.merger->MergePointsOn();
    fs.vtk.merger->AddInput( fs.vtk.clipper1->GetOutput() );
    fs.vtk.merger->AddInput( fs.vtk.clipper2->GetOutput() );
    fs.vtk.merger->AddInput( fs.vtk.clipper3->GetOutput() );
    */

    fs.vtk.merger->Update();
    printf("     merger: number of points = %d\n",(int)fs.vtk.merger->GetOutput()->GetNumberOfPoints());
    printf("     merger: number of cells  = %d\n",(int)fs.vtk.merger->GetOutput()->GetNumberOfCells());
  }

  // test file
  if( debug_file == 1){
    vtkUnstructuredGrid* data;
    if (use_vtk_extract){
      data = fs.vtk.extract->GetOutput();
    }else{
      data = fs.vtk.merger->GetOutput();
    }
    cout << "    unstructured grid: cells  = " << data->GetNumberOfCells() << endl;
    cout << "    unstructured grid: points = " << data->GetNumberOfPoints() << endl;
    // Write vtu file
    printf("\nwriting unstructured grid data...\n");
    std::string filename = "test_init_snapshot.vtu";
    // creates writer
    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(data);
#else
    writer->SetInputData(data);
#endif

    writer->SetDataModeToAscii();
    writer->Write();
    //clean up
    writer->Delete();
    printf("snapshot written to file: %s\n\n",filename.c_str());
  }

  /*
  // this would use a cleaning/merging of the mesh, but it needs to be converted to polydata
  // converts unstructured grid to polydata
  fs.vtk.geometryFilter = vtkGeometryFilter::New();
  if (use_vtk_extract){
    fs.vtk.geometryFilter->SetInput( fs.vtk.extract->GetOutput());
  }else{
    fs.vtk.geometryFilter->SetInput( fs.vtk.merger->GetOutput());
  }
  fs.vtk.geometryFilter->Update();

  fs.vtk.polydata = fs.vtk.geometryFilter->GetOutput();

  // Removes any duplicate points (works only with polydata)
  fs.vtk.cleanFilter = vtkCleanPolyData::New();
  fs.vtk.cleanFilter->SetInput(fs.vtk.polydata);

  // determines approximate minimum point distance to set a merge tolerance
  // has not effect...?
  //fs.vtk.cleanFilter->SetAbsoluteTolerance( dist_min );
  fs.vtk.cleanFilter->PointMergingOn();
  fs.vtk.cleanFilter->Update();
  printf("     cleaner: number of points = %d\n",(int)fs.vtk.cleanFilter->GetOutput()->GetNumberOfPoints());
  printf("     cleaner: number of cells  = %d\n\n",(int)fs.vtk.cleanFilter->GetOutput()->GetNumberOfCells());
  */

  /*
  // we're only interested in rendernig the surface of the model
  // clipper surface
  fs.vtk.clippersurface = vtkDataSetSurfaceFilter::New();
  if (use_vtk_extract){
    fs.vtk.clippersurface->SetInputConnection(0, fs.vtk.extract->GetOutputPort());
  }else{
    fs.vtk.clippersurface->SetInputConnection(0, fs.vtk.merger->GetOutputPort());
  }
  fs.vtk.clippersurface->Update();
  printf("     surface: number of points = %d\n",(int)fs.vtk.clippersurface->GetOutput()->GetNumberOfPoints());
  printf("     surface: number of cells  = %d\n\n",(int)fs.vtk.clippersurface->GetOutput()->GetNumberOfCells());
  */

  // cell connectivity mapping
  fs.vtk.mapMesh3D = vtkDataSetMapper::New();

  // selects input
  if (use_vtk_extract){
    fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.extract->GetOutputPort());
  }else{
    //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.clip1->GetOutputPort());
    fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.merger->GetOutputPort());
  }
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.volume3D->GetProducerPort());
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.cleanFilter->GetOutputPort());
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.clippersurface->GetOutputPort());

  // coloring
  fs.vtk.mapMesh3D->SetLookupTable(fs.vtk.lut);
  fs.vtk.mapMesh3D->SetColorModeToMapScalars();
  fs.vtk.mapMesh3D->SetScalarModeToUsePointData();
  fs.vtk.mapMesh3D->ScalarVisibilityOn();
  fs.vtk.mapMesh3D->UseLookupTableScalarRangeOn();

  //actor
  fs.vtk.actor3D = vtkActor::New();
  fs.vtk.actor3D->SetMapper(fs.vtk.mapMesh3D);
  fs.vtk.actor3D->GetProperty()->SetRepresentationToSurface();
  //fs.vtk.actor3D->GetProperty()->SetEdgeVisibility(1);

  // 3D actor
  fs.vtk.ren->AddActor(fs.vtk.actor3D);

  // legend for colors
  fs.vtk.legendcolor = vtkScalarBarActor::New();
  fs.vtk.legendcolor->SetLookupTable(fs.vtk.lut);
  fs.vtk.legendcolor->SetTitle( "vnorm (m/s)" );
  fs.vtk.legendcolor->SetWidth( 0.1 );
  //fs.vtk.legendcolor->GetTitleTextProperty()->SetFontSize( 16 );
  //fs.vtk.legendcolor->GetLabelTextProperty()->SetFontSize( 8 );
  fs.vtk.ren->AddActor(fs.vtk.legendcolor);

  // text
  fs.vtk.text->SetInput( "...update view with mouse/keyboard, then press <space> to continue " );

  // renders window (with camera reset)
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 4;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();
  printf("     done\n\n");

  // halt simulation, wait for user interaction
  fs.vtk.mutex->Lock();
  fs.vtk.haltOn = 1;
  fs.vtk.mutex->Unlock();
  sync_halt_simulation();

  // starts timer
  fs.timer->StartTimer();
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int* it_h,float* time_h, float* data) {
  TRACE("visualize_vtkdata");

  float bounds[2];
  double dmin,dmax;

  float time = *time_h;

  global_it_step = *it_h;

  // terminal output
  //printf("     visual: it = %d - simulation time = %f (s)\n",global_it_step,time);

  // time for calculating new wavefield
  fs.timer->StopTimer();
  double timeInSeconds = fs.timer->GetElapsedTime();
  double timeTotalRender = fs.vtk.ren->GetLastRenderTimeInSeconds();
  fs.timer->StartTimer();

  // waits for rendering process (to finish with old wavefield)
  int wait_sec;
  if (timeTotalRender > 20.0){
    // in case last rendering time was longer than default 20s, then increase the waiting time
    wait_sec = (int) 1.5 * timeTotalRender;
  }else{
    wait_sec = 20;
  }
  sync_rendering(wait_sec);

  // frames per second (based on computing new wavefield
  double fps = 1.0/timeInSeconds;
  // updates time string
  char inputString[180];
  sprintf(inputString,"time step: %d - simulation time: %6.3f (s) / fps: %f",global_it_step,time,fps);
  fs.vtk.text->SetInput(inputString);

  // updates data
  if(SHOW_VOLUMEDATA == 1 ){
    for(int i=0;i<fs.volume.np;i++) {
      fs.vtk.data_array3D->SetValue(i,data[i]);
      // sets peak values
      // note: data is the velocity norm, thus strictly positive
      if (data[i] > global_data_array_peak[i]) global_data_array_peak[i] = data[i];
    }
    // mark as modified to update rendering
    fs.vtk.data_array3D->Modified();

    // Find min and max
    fs.vtk.data_array3D->GetValueRange(bounds);
    dmin = bounds[0];
    dmax = bounds[1];

    // adjusts color maximum
    if( fs.vtk.colorAdjustOn ){
      if( gcolor_max < 0.0 ) gcolor_max = 1.e-10;
      gcolor_max = dmax;
      gcolor_incr = 0.05 * dmax;
      if( fs.vtk.icolor >= 2 ) gcolor_min = 1.e-3*gcolor_max;
      fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
      fs.vtk.lut->Build();
    }

    // adjust contour
    if (CONTOUR_VOLUME){
      if( fs.vtk.contour3DActor->GetVisibility() == 1){
        dmin = gcolor_min;
        dmax = gcolor_max;
        fs.vtk.contour3D->SetValue(0, dmin + (dmax-dmin)*0.1);
        fs.vtk.contour3D->SetValue(1, dmin + (dmax-dmin)*0.3);
        fs.vtk.contour3D->SetValue(2, dmin + (dmax-dmin)*0.4);
        fs.vtk.contour3D->SetValue(3, dmin + (dmax-dmin)*0.5);
        fs.vtk.contour3D->SetValue(4, dmin + (dmax-dmin)*0.6);
        fs.vtk.contour3D->SetValue(5, dmin + (dmax-dmin)*0.7);
        fs.vtk.contour3D->SetValue(6, dmin + (dmax-dmin)*0.9);
      }
    }
    printf("     data  : min = %f max = %f \n",dmin,dmax);
  }

  // updates window
  // note: does not wait for rendering to be finished, but continue with computing new wavefield while
  //       separate thread can work on rendering new frame
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 1;
  fs.vtk.mutex->Unlock();

  // saves snapshot image
  if (fs.vtk.jpegImageOn){
    fs.vtk.mutex->Lock();
    fs.vtk.do_render = 3;
    fs.vtk.mutex->Unlock();
    // waits for render and saves images
    sync_rendering();
  }

  // time for rendering
  fs.timer->StopTimer();
  double time_renderer = fs.timer->GetElapsedTime();
  fs.timer->StartTimer();

  // terminal output
  printf("     visual: %s \n",inputString);
  printf("     timer : time waited for rendering       = %f (s) \n", time_renderer);
  printf("     timer : time for total window rendering = %f (s) \n\n", timeTotalRender);
  //printf("     timer : time for    2Dmapper rendering  = %f (s) \n", fs.vtk.mapMesh2D->GetTimeToDraw());
  //printf("     timer : time for    2Dcontour rendering = %f (s) \n", fs.vtk.contourMapper->GetTimeToDraw());
  //printf("     timer : time for    3Dmapper rendering  = %f (s) \n", fs.vtk.mapMesh3D->GetTimeToDraw());
  //printf("     timer : time for    3Dcontour rendering = %f (s) \n", fs.vtk.contour3DMapper->GetTimeToDraw());

  // checks if simulation halted
  sync_halt_simulation();
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)(int *do_restart_h) {
  TRACE("finish_vtkwindow");

  float bounds[2];
  double min,max;

  // initializes
  *do_restart_h = 0;

  // check if interactor thread still running
  if (vtk_thread_info.finished == false){

    // sets halt
    fs.vtk.mutex->Lock();
    fs.vtk.haltOn = 1;
    fs.vtk.mutex->Unlock();
    sync_halt_simulation();

    // updates data
    if(SHOW_VOLUMEDATA == 1 ){
      // text
      fs.vtk.text->SetInput( "...simulation done, showing peak ground-motion values, press <space> to continue" );

      for(int i=0;i<fs.volume.np;i++) {
        // updates with peak values
        fs.vtk.data_array3D->SetValue(i,global_data_array_peak[i]);
      }
      // mark as modified to update rendering
      fs.vtk.data_array3D->Modified();

      // Find min and max
      fs.vtk.data_array3D->GetValueRange(bounds);
      min = bounds[0];
      max = bounds[1];

      // adjusts color maximum
      if( fs.vtk.colorAdjustOn ){
        if( gcolor_max < 0.0 ) gcolor_max = 1.e-10;
        gcolor_max = max;
        gcolor_incr = 0.05 * max;
        if( fs.vtk.icolor >= 2 ) gcolor_min = 1.e-3*gcolor_max;
        fs.vtk.lut->SetTableRange( gcolor_min, gcolor_max );
        fs.vtk.lut->Build();
      }
      fs.vtk.legendcolor->SetTitle( "PG vnorm (m/s)" );

      // adjust contour
      if (CONTOUR_VOLUME){
        if( fs.vtk.contour3DActor->GetVisibility() == 1){
          min = gcolor_min;
          max = gcolor_max;
          fs.vtk.contour3D->SetValue(0, min + (max-min)*0.1);
          fs.vtk.contour3D->SetValue(1, min + (max-min)*0.3);
          fs.vtk.contour3D->SetValue(2, min + (max-min)*0.4);
          fs.vtk.contour3D->SetValue(3, min + (max-min)*0.5);
          fs.vtk.contour3D->SetValue(4, min + (max-min)*0.6);
          fs.vtk.contour3D->SetValue(5, min + (max-min)*0.7);
          fs.vtk.contour3D->SetValue(6, min + (max-min)*0.9);
        }
      }

      // terminal output
      printf("\nvtk: showing Peak Ground-motion values \n");
      printf("     Peak Ground-motion vnorm (m/s)\n");
      printf("     data  : min = %f max = %f \n\n",min,max);
    }
    // render window
    fs.vtk.mutex->Lock();
    fs.vtk.do_render = 1;
    fs.vtk.mutex->Unlock();
    sync_rendering();

    // sets halt
    fs.vtk.mutex->Lock();
    fs.vtk.haltOn = 1;
    fs.vtk.mutex->Unlock();
    sync_halt_simulation();

    // text
    fs.vtk.text->SetInput( "...finishing, <y> to restart simulation, then press <space> to continue" );
    // render window
    fs.vtk.mutex->Lock();
    fs.vtk.do_render = 1;
    fs.vtk.mutex->Unlock();
    sync_rendering();

    // terminal output
    if (fs.vtk.do_restart){
      printf("\nvtk: restart flag is turned on\n");
      printf("     press <y> to turn restarting off...\n\n");
    }else{
      printf("\nvtk: restart flag is turned off\n");
      printf("     press <y> to turn restarting on...\n\n");
    }

    // sets halt
    fs.vtk.mutex->Lock();
    fs.vtk.haltOn = 1;
    fs.vtk.mutex->Unlock();
    sync_halt_simulation();

    // checks if restart flag set
    if (fs.vtk.do_restart){
      printf("vtk: restart flag set, do restart simulation...\n");
      // reset
      fs.vtk.legendcolor->SetTitle( "vnorm (m/s)" );
      for(int i=0;i<fs.volume.np;i++) {
        fs.vtk.data_array3D->SetValue(i,0.0);
        global_data_array_peak[i] = 0.0f;
      }

      // done, return to restart time loop
      *do_restart_h = 1;
      return;
    }

    // finish thread
    terminate_vtk_thread();
  }

  // clean up
  clean_vtk_arrays();

  // wait period
  sleep(1);
}

#else

/* ----------------------------------------------------------------------------------------------- */

// empty stubs

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(initialize_vtkwindow,INITIALIZE_VTKWINDOW)(int*) {
  fprintf(stderr,"ERROR: VTK_VIS enabled without VTK Support. "
                 "To enable VTK support, reconfigure with --enable-vtk and --with-vtk-version=.. flags.\n");
  exit(1);
};

extern "C"
void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float*, float*, float*) {};

extern "C"
void FC_FUNC_(prepare_vtkreceivers,PREPARE_VTKRECEIVERS)(int*, float*, float*, float*, int*, char*) {};

extern "C"
void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int* ,float* ,float* ,float*, int*, int* ) {};

extern "C"
void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int* ,float* ,float* ,float* ,int* ,int* ) {};

extern "C"
void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int*,float*,float*) {};

extern "C"
void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)(int*) {};

#endif // HAVE_VTK

