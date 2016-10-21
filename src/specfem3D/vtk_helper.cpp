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

#include "config.h"

#ifdef HAVE_VTK
#pragma message ("\nCompiling with: HAVE_VTK enabled\n")

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>

using namespace std;

#include <vtkCellArray.h>
#include <vtkProperty.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolygon.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkMutexLock.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay2D.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLDataSetWriter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkTimerLog.h>
#include <vtkPlane.h>
#include <vtkBox.h>
#include <vtkBoxClipDataSet.h>
#include <vtkBoundingBox.h>
#include <vtkClipPolyData.h>
#include <vtkClipVolume.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkLookupTable.h>
#include <vtkLegendBoxActor.h>
#include <vtkScalarBarActor.h>
#include <vtkLegendBoxActor.h>
#include <vtkScalarBarActor.h>
#include <vtkClipDataSet.h>
#include <vtkAppendFilter.h>
#include <vtkJPEGWriter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSphereSource.h>
#include <vtkGeometryFilter.h>
#include <vtkTableBasedClipDataSet.h>
#include <vtkQuad.h>
#include <vtkContourFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkGeometryFilter.h>

/* ----------------------------------------------------------------------------------------------- */

/******************************************************************************/

// USER PARAMETERS
int CONTOUR_FREESURFACE = 1;  // creates contour for free surface topography
int CONTOUR_VOLUME = 1;       // creates contour for volumetric data

// color range to color scale plotting
double gcolor_min  = 0.0;     // minimum value
double gcolor_max  = 0.005;   // maximum value
double gcolor_incr = 0.001;   // increments

/******************************************************************************/

// global flags
int SHOW_FREESURFACE = 0;
int SHOW_VOLUMEDATA = 0;
int SHOW_GPU_TEXT = 0;

/* ----------------------------------------------------------------------------------------------- */

// macros
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define exit_error(msg) { fputs (msg,stderr); exit(1); }
// debugging
//#define TRACE(x) printf("%s\n",x);
//#define TRACE(x) printf("thread %lu: %s\n",pthread_self(),x);
#define TRACE(x)


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
void sync_rendering();
void sync_halt_simulation();

extern "C" void FC_FUNC_(initialize_vtkwindow,INITIALIZE_VTKWINDOW)(int*);
extern "C" void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float*,float*,float*);
extern "C" void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int*,float*,float*,float*,int*,int*);
extern "C" void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int*,float*,float*,float*,int*,int*);
extern "C" void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int*,float*,float*);
extern "C" void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)(int*);

/* ----------------------------------------------------------------------------------------------- */


class Visualization {
  public:
    // 2D surface data
    vtkPoints* points2D;
    vtkFloatArray* data_array2D;
    vtkCellArray* cells2D;
  
    vtkUnstructuredGrid* volume2D;

    vtkDataSetMapper* mapMesh2D;
    vtkActor* actor2D;
  
    // 3D volume data
    vtkPoints* points3D;
    vtkFloatArray* data_array3D;
    vtkCellArray* cells;

    vtkUnstructuredGrid* volume3D;

    vtkDataSetMapper* mapMesh3D;
    vtkActor* actor3D;

    // clipping
    vtkPlane* clipPlane1;
    vtkPlane* clipPlane2;

    //vtkTableBasedClipDataSet* clip1;
    //vtkTableBasedClipDataSet* clip2;
    //vtkTableBasedClipDataSet* clip1m;
    vtkClipDataSet* clip1;
    vtkClipDataSet* clip2;
    vtkClipDataSet* clip1m;

    // visualizing
    vtkAppendFilter* merger;
    vtkDataSetSurfaceFilter* clippersurface;

    vtkGeometryFilter* geometryFilter;
    vtkPolyData* polydata;
    vtkCleanPolyData* cleanFilter;

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

    vtkInteractorStyleTrackballCamera *style;

    vtkMutexLock *mutex;
    volatile int do_render;

    // camera
    vtkCamera* camera;
    double pcam[3]; // camera position
    double rclip[2]; // clipping range
    double pfocal[3]; // focal point
  
    // source sphere
    vtkSphereSource* sphere;
    vtkPolyDataMapper* mapperSphere;
    vtkActor* actorSphere;  
    double pos_source[3]; // source location

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
  "  6            - show freesurface (toggle on/off)\n"
  "  7            - show freesurface contour (toggle on/off)\n"
  "  8            - show volume (toggle on/off)\n"
  "  9            - show volume contour (toggle on/off)\n"
  "  0            - show source sphere (toggle on/off)\n"
  "\n"
  "  <space>      - halt/continue simulation\n"
  "  <escape>,q,e - abort and exit simulation\n"
  "  y            - yes, restart simulation (after the end of the time loop)\n";

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
  "    q - exit\n";

// time step it from simulation
int global_it_step = 0;
// array for storing peak values (using a simple float array to avoid memory overhead from vtkFloatArray...)
float* global_data_array_peak;

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
        //fs.vtk.ren->ResetCamera(); // only resets zoom effect
        // range
        fs.vtk.camera->SetClippingRange( fs.vtk.rclip );
        // focal point
        fs.vtk.camera->SetFocalPoint( fs.vtk.pfocal );
        // position
        fs.vtk.camera->SetPosition( fs.vtk.pcam );
        // view
        fs.vtk.camera->SetViewAngle( 30.0 );
        fs.vtk.camera->SetViewUp( 0.0, 0.0, 1.0 );
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
      // toggles freesurface visibility
      if (key == "6"){
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
      if (key == "7"){
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
      // toggles volume visibility
      if (key == "8"){
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
      if (key == "9"){
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
      // toggles source glyph visibility
      if (key == "0"){
        if( fs.vtk.actorSphere->GetVisibility() == 1){
          fs.vtk.actorSphere->SetVisibility( 0 );
          printf("\ntoggle off: source glyph\n");
        }else{
          fs.vtk.actorSphere->SetVisibility( 1 );
          printf("\ntoggle on: source glyph\n");
        }
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

    static CommandSubclass2 *New(){
      return new CommandSubclass2;
    }

    void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId), void *vtkNotUsed(callData)){
      vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
      // note: this execute function will be executed only by the vtk thread which also created the interactor
      // checks for rendering
      if (fs.vtk.do_render == 1){
        // render window
        TRACE("execute: rendering");
        // update
        fs.vtk.renWin->Render();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: rendering done");
      }
      // check for terminating interactive
      if (fs.vtk.do_render == 2){
        // stop interactor
        TRACE("execute: terminating interactor");
        // stoping
        fs.vtk.iren->TerminateApp();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: terminate done");
      }
      // render image
      if (fs.vtk.do_render == 3){
        // stop interactor
        TRACE("execute: rendering image");
        // update
        fs.vtk.renWin->Render();
        save_snapshot_jpg();
        fs.vtk.mutex->Lock();
        fs.vtk.do_render = 0;
        fs.vtk.mutex->Unlock();
        TRACE("execute: rendering image done");
      }
    }
};

/* ----------------------------------------------------------------------------------------------- */

// threading version

#include <pthread.h>

struct threadInfo{
  volatile bool finished;
  volatile bool started;
  volatile bool init;
  pthread_mutex_t mutex;
};

struct threadInfo vtk_thread_info;

pthread_t vtk_thread;

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

// helper function

/* ----------------------------------------------------------------------------------------------- */

void interactor_usage(){
  cout << endl;
  cout << info_string << endl;
  cout << info_string_defaults << endl;
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
  writer->SetInput(fs.vtk.volume3D);
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

  // camera position
  fs.vtk.camera = vtkCamera::New();

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
  fs.vtk.ren = vtkRenderer::New();
  fs.vtk.ren->SetActiveCamera(fs.vtk.camera);
  
  // Background color white
  fs.vtk.bgBlackOn = 0;
  fs.vtk.ren->SetBackground(1,1,1);

  // text actors
  // GPU flag
  fs.vtk.textGPU = vtkTextActor::New();
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
  fs.vtk.help->GetTextProperty()->SetFontSize ( 16 );
  fs.vtk.help->SetPosition( 10, 80 );
  fs.vtk.help->GetTextProperty()->SetColor( 0.5, 0.5, 0.5 );
  fs.vtk.help->SetInput( (const char*) info_string.c_str() );
  fs.vtk.ren->AddActor2D( fs.vtk.help );

  // progress text
  fs.vtk.text = vtkTextActor::New();
  fs.vtk.text->GetTextProperty()->SetFontSize ( 16 );
  fs.vtk.text->SetPosition( 10, 530 );
  fs.vtk.text->SetInput( "...initializing data " );
  fs.vtk.text->GetTextProperty()->SetColor( 0.5,0.5,0.5 );
  fs.vtk.ren->AddActor2D( fs.vtk.text );

  // color table
  int tableSize = 256;
  fs.vtk.icolor = 0; // from blue to red
  fs.vtk.colorAdjustOn = 1; // automatic adjust
  
  fs.vtk.lut = vtkLookupTable::New();
  fs.vtk.lut->SetNumberOfColors(tableSize);

  // sets (default) rainbow color scale
  set_color_scale(fs.vtk.icolor);

  // render window
  fs.vtk.renWin = vtkRenderWindow::New();
  fs.vtk.renWin->AddRenderer(fs.vtk.ren);

  fs.vtk.renWin->SetPosition(500,0);
  fs.vtk.renWin->SetSize(900,600);
  
  // timer
  fs.timer = vtkTimerLog::New();

  // mutex
  fs.vtk.mutex = vtkMutexLock::New();

  // initializes flags
  fs.vtk.haltOn = 0;
  fs.vtk.jpegImageOn = 0;
  fs.vtk.do_restart = 0;
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

  // Initialize must be called prior to creating timer events.
  fs.vtk.iren->Initialize();

  // callback function
  fs.vtk.iren->AddObserver(vtkCommand::KeyPressEvent, &mykey, &MyInteractor::KeypressCallbackFunction);

  // below call might crash if window is rendered before interactor is initialized
  // error in X11: GLXBadCurrentWindow
  interactor_usage();

  // rendering will be called within this thread by a timer callback method
  vtkSmartPointer<CommandSubclass2> timerCallback = vtkSmartPointer<CommandSubclass2>::New();
  fs.vtk.iren->AddObserver( vtkCommand::TimerEvent, timerCallback );
  // timer in milliseconds
  fs.vtk.iren->CreateRepeatingTimer(500);

  // indicate that window initialization is done
  pthread_mutex_lock(&vtk_thread_info.mutex);
  vtk_thread_info.init = false;
  pthread_mutex_unlock(&vtk_thread_info.mutex);

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
    fs.vtk.cells2D->Delete();

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
    fs.vtk.cells->Delete();

    fs.vtk.legendcolor->Delete();
    fs.vtk.lut->Delete();

    fs.vtk.clipPlane1->Delete();
    fs.vtk.clipPlane2->Delete();

    fs.vtk.clip1->Delete();
    fs.vtk.clip1m->Delete();
    fs.vtk.clip2->Delete();

    fs.vtk.merger->Delete();

    fs.vtk.clippersurface->Delete();

    fs.vtk.geometryFilter->Delete();
    fs.vtk.cleanFilter->Delete();

    fs.vtk.actor3D->Delete();
    fs.vtk.mapMesh3D->Delete();

    free(global_data_array_peak);
  }

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


void sync_rendering(){
  TRACE("sync_rendering");

  // waits for vtk rendering to finish

  // time out
  time_t current = time(NULL);
  time_t seconds = 20;
  time_t endtime = current + seconds;

  //printf("render loop time is    : %s", ctime(&current));
  //printf("render loop endtime is : %s", ctime(&endtime));

  while(fs.vtk.do_render != 0){
    sleep(0.1);
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
      sleep(0.1);
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

  printf("vtk: initializing window...\n");
  // initializes vtk
  init_vtk();

  //sleep(1);

  // vtk window and interactor
  // starts interactor within separate thread, otherwise this will halt until interactor is terminated
  TRACE("vtk creating thread");
  create_vtk_thread();

  //sleep(1);

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

  TRACE("initialized vtk window successfully");
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float* xs_x,float* xs_y, float* xs_z) {
  TRACE("prepare vtk source glyph ...");

  double xyz[3];
  
  // sets source location
  fs.vtk.pos_source[0] = *xs_x;
  fs.vtk.pos_source[1] = *xs_y;
  fs.vtk.pos_source[2] = *xs_z;

  // terminal output
  printf("vtk: source sphere\n");
  printf("     sphere location: %f %f %f \n\n",fs.vtk.pos_source[0],fs.vtk.pos_source[1],fs.vtk.pos_source[2]);
  
  // creates sphere around source location
  fs.vtk.sphere = vtkSphereSource::New();
  fs.vtk.sphere->SetCenter(fs.vtk.pos_source);
  fs.vtk.sphere->SetRadius(0.02);

  fs.vtk.mapperSphere = vtkPolyDataMapper::New();
  fs.vtk.mapperSphere->SetInputConnection(fs.vtk.sphere->GetOutputPort());
 
  fs.vtk.actorSphere = vtkActor::New();
  fs.vtk.actorSphere->SetMapper(fs.vtk.mapperSphere);

  // adds actor
  fs.vtk.ren->AddActor(fs.vtk.actorSphere);

  // render window
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 1;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();
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

  // exagerates vertical dimension (1.0 == no exageration, 1.5 == 50% exageration...)
  static float VERTICAL_EXAGERATION = 1.5;

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
  printf("     topo: min = %f max = %f\n\n", min,max);
  
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
  fs.vtk.cells2D = vtkCellArray::New();
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
    fs.vtk.cells2D->InsertNextCell(quad);
  }
  quad->Delete();

  fs.vtk.volume2D = vtkUnstructuredGrid::New();
  // points
  fs.vtk.volume2D->SetPoints(fs.vtk.points2D);
  fs.vtk.volume2D->GetPointData()->SetScalars(fs.vtk.data_array2D);
  // cells
  fs.vtk.volume2D->SetCells(VTK_QUAD, fs.vtk.cells2D);

  // contour iso-surfacing
  if (CONTOUR_FREESURFACE){
    fs.vtk.contour = vtkContourFilter::New();
    fs.vtk.contour->SetInput(  fs.vtk.volume2D );
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
    fs.vtk.contourActor->SetScale( 1.0, 1.0, VERTICAL_EXAGERATION );
    fs.vtk.contourActor->GetProperty()->SetColor( 0.5, 0.5, 0.5 );
    
    fs.vtk.ren->AddActor(fs.vtk.contourActor);
  }

  // mapping
  fs.vtk.mapMesh2D = vtkDataSetMapper::New();
  fs.vtk.mapMesh2D->SetInput( fs.vtk.volume2D );
  
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
  fs.vtk.actor2D->SetScale( 1.0, 1.0, VERTICAL_EXAGERATION );

  // 2D actor
  fs.vtk.ren->AddActor(fs.vtk.actor2D);

  // render window
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 1;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();
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
  int id1,id2,id3,id4,id5,id6,id7,id8;

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
  double min = data_bounds[0];
  double max = data_bounds[1];
  printf("     volume data: min = %f max = %f\n\n",min,max);

  // Find min and max
  fs.vtk.points3D->GetBounds(model_bounds);
  double xmin = model_bounds[0];
  double xmax = model_bounds[1];
  double ymin = model_bounds[2];
  double ymax = model_bounds[3];
  double zmin = model_bounds[4];
  double zmax = model_bounds[5];
  printf("     model dimension: xmin/xmax = %f / %f\n",xmin,xmax);
  printf("                      ymin/ymax = %f / %f\n",ymin,ymax);
  printf("                      zmin/zmax = %f / %f\n\n",zmin,zmax);

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

  // adjust sphere size
  fs.vtk.sphere->SetRadius(0.02*fabs(zmax-zmin));
  fs.vtk.sphere->Modified();
  
  //
  // unstructured grid
  //
  // creates cell connectivity
  // cells
  fs.vtk.cells = vtkCellArray::New();
  vtkHexahedron* hex = vtkHexahedron::New();
  for(int ispec=0;ispec<fs.volume.nspec;ispec++){
    id1 = vol_conn[0+ispec*8];
    id2 = vol_conn[1+ispec*8];
    id3 = vol_conn[2+ispec*8];
    id4 = vol_conn[3+ispec*8];
    id5 = vol_conn[4+ispec*8];
    id6 = vol_conn[5+ispec*8];
    id7 = vol_conn[6+ispec*8];
    id8 = vol_conn[7+ispec*8];
    hex->GetPointIds()->SetId(0,id1);
    hex->GetPointIds()->SetId(1,id2);
    hex->GetPointIds()->SetId(2,id3);
    hex->GetPointIds()->SetId(3,id4);
    hex->GetPointIds()->SetId(4,id5);
    hex->GetPointIds()->SetId(5,id6);
    hex->GetPointIds()->SetId(6,id7);
    hex->GetPointIds()->SetId(7,id8);
    fs.vtk.cells->InsertNextCell(hex);      
  }
  hex->Delete();
  
  fs.vtk.volume3D = vtkUnstructuredGrid::New();
  // points
  fs.vtk.volume3D->SetPoints(fs.vtk.points3D);
  fs.vtk.volume3D->GetPointData()->SetScalars(fs.vtk.data_array3D);
  // cells
  fs.vtk.volume3D->SetCells(VTK_HEXAHEDRON, fs.vtk.cells);

  // contour iso-surfacing
  if (CONTOUR_VOLUME){
    min = gcolor_min;
    max = gcolor_max;
    fs.vtk.contour3D = vtkContourFilter::New();
    fs.vtk.contour3D->SetInput(  fs.vtk.volume3D );
    fs.vtk.contour3D->SetNumberOfContours(7);
    fs.vtk.contour3D->SetValue(0, min + (max-min)*0.1);
    fs.vtk.contour3D->SetValue(1, min + (max-min)*0.3);
    fs.vtk.contour3D->SetValue(2, min + (max-min)*0.4);
    fs.vtk.contour3D->SetValue(3, min + (max-min)*0.5);
    fs.vtk.contour3D->SetValue(4, min + (max-min)*0.6);
    fs.vtk.contour3D->SetValue(5, min + (max-min)*0.7);
    fs.vtk.contour3D->SetValue(6, min + (max-min)*0.9);

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

  //
  // clip box
  //
  // bounds
  double bb_all[6] = { xmin,xmax,ymin,ymax,zmin,zmax };
  
  // vector to source
  double v[3];
  v[0] = fs.vtk.pos_source[0];
  v[1] = fs.vtk.pos_source[1];
  v[2] = fs.vtk.pos_source[2];

  // source location vector / perpendicular vectors
  //
  // vectorproduct(vector1, vector2, product)
  // calculates the vector product of vector1 and vector2
  //    product(1) = vector1(2)*vector2(3) - vector1(3)*vector2(2)
  //    product(2) = - vector1(1)*vector2(3) + vector1(3)*vector2(1)
  //    product(3) = vector1(1)*vector2(2) - vector1(2)*vector2(1)
  double p1[3],p2[3];
  p1[0] = 0.0;
  p1[1] = 0.0;
  p1[2] = 1.0;
    
  p2[0] = 0.0;
  p2[1] = 1.0;
  p2[2] = 0.0;

  // normals
  double pn1[3],pn2[3];
  /*
  pn1[0] = v[1]*p1[2] - v[2]*p1[1];
  pn1[1] = -v[0]*p1[2] + v[2]*p1[0];
  pn1[2] = v[0]*p1[1] - v[1]*p1[0];

  pn2[0] = v[1]*p2[2] - v[2]*p2[1];
  pn2[1] = -v[0]*p2[2] + v[2]*p2[0];
  pn2[2] = v[0]*p2[1] - v[1]*p2[0];
  */
  pn1[0] = -1.0;
  pn1[1] = 0.0;
  pn1[2] = 0.0;

  pn2[0] = 0.0;
  pn2[1] = -1.0;
  pn2[2] = 0.0;
  
  // flips normal depending on location of source
  if( fabs(v[0]-xmin) < fabs(v[0]-xmax) ){
    pn1[0] *= -1.0;
    pn1[1] *= -1.0;
    pn1[2] *= -1.0;
  }
  if( fabs(v[1]-ymin) < fabs(v[1]-ymax) ){
    pn2[0] *= -1.0;
    pn2[1] *= -1.0;
    pn2[2] *= -1.0;
  }

  // plane through source location
  fs.vtk.clipPlane1 = vtkPlane::New();
  fs.vtk.clipPlane1->SetNormal( pn1 );
  fs.vtk.clipPlane1->SetOrigin( v );

  fs.vtk.clipPlane2 = vtkPlane::New();
  fs.vtk.clipPlane2->SetNormal( pn2 );
  fs.vtk.clipPlane2->SetOrigin( v );

  // table based clipping
  //fs.vtk.clip1 = vtkTableBasedClipDataSet::New();
  // note: table base clipping might miss surfaces when exactly aligned with geometry
  //       the clipdataset avoids some of these problems...
  fs.vtk.clip1 = vtkClipDataSet::New();
  fs.vtk.clip1->SetInput( fs.vtk.volume3D );
  fs.vtk.clip1->SetClipFunction( fs.vtk.clipPlane1 );
  fs.vtk.clip1->Update();

  //fs.vtk.clip1m = vtkTableBasedClipDataSet::New();
  fs.vtk.clip1m = vtkClipDataSet::New();
  fs.vtk.clip1m->SetInput( fs.vtk.volume3D );
  fs.vtk.clip1m->SetClipFunction( fs.vtk.clipPlane1 );
  fs.vtk.clip1m->InsideOutOn();
  fs.vtk.clip1m->Update();

  //fs.vtk.clip2 = vtkTableBasedClipDataSet::New();
  fs.vtk.clip2 = vtkClipDataSet::New();
  fs.vtk.clip2->SetInput( fs.vtk.clip1m->GetOutput() );
  fs.vtk.clip2->SetClipFunction( fs.vtk.clipPlane2 );
  fs.vtk.clip2->Update();

  // merges mesh parts
  fs.vtk.merger = vtkAppendFilter::New();
  fs.vtk.merger->MergePointsOn();
  fs.vtk.merger->AddInput( fs.vtk.clip1->GetOutput() );
  fs.vtk.merger->AddInput( fs.vtk.clip2->GetOutput() );
  fs.vtk.merger->Update();

  // converts unstructuredgrid to polydata
  fs.vtk.geometryFilter = vtkGeometryFilter::New();
  fs.vtk.geometryFilter->SetInput( fs.vtk.merger->GetOutput());
  fs.vtk.geometryFilter->Update();

  fs.vtk.polydata = fs.vtk.geometryFilter->GetOutput();

  // Removes any duplicate points (works only with polydata)
  fs.vtk.cleanFilter = vtkCleanPolyData::New();
  fs.vtk.cleanFilter->SetInput(fs.vtk.polydata);

  // determines approximate minimum point distance to set a merge tolerance
  // has not effect...?
  //double dist_min = 1.e21;
  //dist_min = MIN(dist_min,xmax - xmin);
  //dist_min = MIN(dist_min,ymax - ymin);
  //dist_min = MIN(dist_min,zmax - zmin);
  //dist_min = dist_min / fs.volume.nspec;
  //printf("estimated minimum distance %f\n",dist_min);
  //fs.vtk.cleanFilter->SetAbsoluteTolerance( dist_min );

  fs.vtk.cleanFilter->PointMergingOn();
  fs.vtk.cleanFilter->Update();

  printf("     mapper: number of points = %d\n",fs.vtk.cleanFilter->GetOutput()->GetNumberOfPoints());
  printf("     mapper: number of cells  = %d\n\n",fs.vtk.cleanFilter->GetOutput()->GetNumberOfCells());

  /*
  // instead of two clips, this would create 3 blocks for merging...
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
  fs.vtk.merger->Update();
  */
  
  // test file
  if( debug_file == 1){
    vtkUnstructuredGrid* data = fs.vtk.merger->GetOutput();
    cout << "    merger: cells  = " << data->GetNumberOfCells() << endl;
    cout << "    merger: points = " << data->GetNumberOfPoints() << endl;
    // Write vtu file
    printf("\nwriting unstructured grid data...\n");
    std::string filename = "test_init_snapshot.vtu";
    // creates writer
    vtkXMLUnstructuredGridWriter* writer = vtkXMLUnstructuredGridWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInput(data);
    writer->SetDataModeToAscii();
    writer->Write();
    //clean up
    writer->Delete();
    printf("snapshot written to file: %s\n\n",filename.c_str());
  }
  
  // clipper surface
  fs.vtk.clippersurface = vtkDataSetSurfaceFilter::New();
  fs.vtk.clippersurface->SetInputConnection(0, fs.vtk.merger->GetOutputPort());
  fs.vtk.clippersurface->Update();
  
  // cell connectivity mapping
  fs.vtk.mapMesh3D = vtkDataSetMapper::New();

  // selects input
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.clippersurface->GetOutputPort());
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.merger->GetOutputPort());
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.volume3D->GetProducerPort());
  //fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.clip1->GetOutputPort());
  fs.vtk.mapMesh3D->SetInputConnection(fs.vtk.cleanFilter->GetOutputPort());

  // coloring
  fs.vtk.mapMesh3D->SetLookupTable(fs.vtk.lut);
  fs.vtk.mapMesh3D->SetColorModeToMapScalars();
  fs.vtk.mapMesh3D->SetScalarModeToUsePointData();
  fs.vtk.mapMesh3D->ScalarVisibilityOn();
  fs.vtk.mapMesh3D->UseLookupTableScalarRangeOn();
  
  //actor
  fs.vtk.actor3D = vtkActor::New();
  fs.vtk.actor3D->SetMapper(fs.vtk.mapMesh3D);
  //fs.vtk.actor3D->GetProperty()->SetRepresentationToSurface();
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

  // reposition the camera, so that actor can be fully seen
  fs.vtk.ren->ResetCamera();

  // text
  fs.vtk.text->SetInput( "...update view with mouse/keyboard, then press <space> to continue " );

  // renders window
  fs.vtk.mutex->Lock();
  fs.vtk.do_render = 1;
  fs.vtk.mutex->Unlock();
  // waits until rendered
  sync_rendering();

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
  double min,max;

  float time = *time_h;

  global_it_step = *it_h;

  // terminal output
  printf("     visual: it = %d - simulation time = %f (s)\n",global_it_step,time);

  // time for calculating new wavefield
  fs.timer->StopTimer();
  double timeInSeconds = fs.timer->GetElapsedTime();
  fs.timer->StartTimer();

  // waits for rendering process (to finish with old wavefield)
  sync_rendering();

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
  printf("     timer : time for rendering = %f (s) \n", time_renderer);
  printf("     data  : min = %f max = %f \n\n",min,max);

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
void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float* xs_x,float* xs_y, float* xs_z) {};

extern "C"
void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int* ,float* ,float* ,float*, int*, int* ) {};

extern "C"
void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int* ,float* ,float* ,float* ,int* ,int* ) {};

extern "C"
void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int*,float*,float*) {};

extern "C"
void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)(int*) {};

#endif // HAVE_VTK

