# This program automates the cross-section cutting for Banana-Doughnut kernels
# Qinya Liu, Caltech, May 2007

#
# We start off by loading some Tcl modules. One is the basic VTK library;
# the second is a package for rendering, and the last includes a set
# of color definitions.
#
package require vtk
package require vtkinteraction
package require vtktesting

# === load the unstructured grid data ===
vtkXMLUnstructuredGridReader kReader
    kReader SetFileName vp.vtu

# data mappers
vtkDataSetMapper kMapper
    kMapper SetInput [kReader GetOutput]
    kMapper SetScalarRange 2000 8000

# Actor
vtkActor kActor
    kActor SetMapper kMapper

# === generate plane cut ====
vtkPlane ewPlane
    ewPlane SetOrigin 385720. 3820000 -30000
    ewPlane SetNormal 0 1 0
vtkCutter ewCut
    ewCut SetInput [kReader GetOutput]
    ewCut SetCutFunction ewPlane
vtkDataSetMapper ewMapper
   ewMapper SetInput [ewCut GetOutput]
   ewMapper InterpolateScalarsBeforeMappingOn
   ewMapper SetScalarRange 2000 8000
vtkActor ewActor
    ewActor SetMapper ewMapper

#=== color bar ====
# NOTE ewMapper variable here
vtkScalarBarActor scalarBar
    scalarBar SetLookupTable [ewMapper GetLookupTable]
    scalarBar SetTitle "P-wave (m/s)"
    [scalarBar GetPositionCoordinate] SetCoordinateSystemToNormalizedViewport
    [scalarBar GetPositionCoordinate] SetValue 0.1 0.05
    scalarBar SetOrientationToVertical
    scalarBar SetWidth 0.1
    scalarBar SetHeight 0.2
    scalarBar SetPosition 0.01 0.2
    scalarBar SetLabelFormat "%-#6.3g"
    [scalarBar GetLabelTextProperty] SetColor 0 0 0
    [scalarBar GetLabelTextProperty] SetFontFamilyToTimes
    [scalarBar GetTitleTextProperty] SetColor 0 0 0
    [scalarBar GetTitleTextProperty] SetFontFamilyToTimes

#=== text title ===
vtkTextActor titleActor
   titleActor SetDisplayPosition 375 100
   titleActor SetInput "Lin model with Harvard model -- EW cross-section at UTM-Y = 3820000 m"
set tprop [titleActor GetTextProperty]
   $tprop SetJustificationToCentered
   $tprop SetColor 0 0 1
   $tprop SetFontFamilyToTimes
   $tprop SetFontSize 16

#=== render window and camera positions ===
vtkRenderer ren1
    ren1 SetBackground 1 1 1
vtkRenderWindow renWin
    renWin AddRenderer ren1
    renWin SetSize 800 500 
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

vtkCamera cam1
    cam1 SetPosition 385720. 200000 -30000. 
    cam1 SetFocalPoint 385720. 3823300. -30000.
    cam1 SetViewUp 0 0 1

#=== scene ===============
ren1 SetActiveCamera cam1

# colorbar only works with --> nsActor
# kActor is the volume

#ren1 AddActor nsActor
ren1 AddActor ewActor
#ren1 AddActor hrActor
#ren1 AddActor kActor
ren1 AddActor2D scalarBar
ren1 AddActor2D titleActor

# I don't understand why, but ResetCamera is very important
ren1 ResetCamera 300000 550000 3572000 4075000 -60000 200000
renWin Render

vtkWindowToImageFilter w2i  
  w2i SetInput renWin
vtkPostScriptWriter writer
  writer SetInputConnection [w2i GetOutputPort]
  writer SetFileName "xc_ew.ps"
  writer Write

wm withdraw .

exit

#set cam2 [ren1 GetActiveCamera]
#puts stdout [ $cam2 GetPosition]
#puts stdout [$cam2 GetFocalPoint]
#puts stdout [$cam2 GetViewUp]
