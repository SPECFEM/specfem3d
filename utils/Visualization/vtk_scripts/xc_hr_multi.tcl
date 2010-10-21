# This program automates the cross-section cutting for Banana-Doughnut kernels
# Qinya Liu and Carl Tape, Caltech, May 2007
#
# We start off by loading some Tcl modules. One is the basic VTK library;
# the second is a package for rendering, and the last includes a set
# of color definitions.
#
package require vtk
package require vtkinteraction
package require vtktesting

# === lookup table ===
vtkLookupTable lut
lut SetNumberOfTableValues 25
lut SetTableValue  0 0.72549 0 0 1
lut SetTableValue  1 0.847059 0 0 1
lut SetTableValue  2 0.964706 0 0 1
lut SetTableValue  3 1 0.0862745 0 1
lut SetTableValue  4 1 0.203922 0 1
lut SetTableValue  5 1 0.32549 0 1
lut SetTableValue  6 1 0.447059 0 1
lut SetTableValue  7 1 0.564706 0 1
lut SetTableValue  8 1 0.686275 0 1
lut SetTableValue  9 1 0.807843 0 1
lut SetTableValue  10 1 0.92549 0 1
lut SetTableValue  11 1 1 0 1
lut SetTableValue  12 1 1 0 1
lut SetTableValue  13 1 1 0 1
lut SetTableValue  14 0.858824 1 0.027451 1
lut SetTableValue  15 0.623529 1 0.0666667 1
lut SetTableValue  16 0.392157 1 0.109804 1
lut SetTableValue  17 0.247059 0.980392 0.211765 1
lut SetTableValue  18 0.117647 0.960784 0.32549 1
lut SetTableValue  19 0 0.92549 0.443137 1
lut SetTableValue  20 0 0.701961 0.65098 1
lut SetTableValue  21 0 0.47451 0.854902 1
lut SetTableValue  22 0 0.282353 0.980392 1
lut SetTableValue  23 0 0.168627 0.909804 1
lut SetTableValue  24 0 0.054902 0.839216 1

# === load the unstructured grid data ===
vtkXMLUnstructuredGridReader kReader1
    kReader1 SetFileName vtu_files/lin_model/vp_low_1.vtu

vtkXMLUnstructuredGridReader kReader2
    kReader2 SetFileName vtu_files/lin_model/vp_low_2.vtu

vtkXMLUnstructuredGridReader kReader3
    kReader3 SetFileName vtu_files/lin_model/vp_low_3.vtu

vtkAppendFilter kReader
  kReader AddInput [kReader1 GetOutput]
  kReader AddInput [kReader2 GetOutput]
  kReader AddInput [kReader3 GetOutput]

# data mappers
vtkDataSetMapper kMapper
    kMapper SetInput [kReader GetOutput]
    kMapper SetScalarRange 4687.47 5729.13

# data Actor
vtkActor kActor
    kActor SetMapper kMapper

# === generate plane cut ====
vtkPlane hrPlane1
    hrPlane1 SetOrigin 385719.8 3823329.2 0.0 
    hrPlane1 SetNormal 0 0 1 
vtkCutter hrCut
    hrCut SetInput [kReader GetOutput]
    hrCut SetCutFunction hrPlane1
vtkDataSetMapper hrMapper
   hrMapper SetInput [hrCut GetOutput]
   hrMapper InterpolateScalarsBeforeMappingOn
   hrMapper SetScalarRange 4687.47 5729.13
   hrMapper SetLookupTable lut
vtkActor hrActor
    hrActor SetMapper hrMapper

# === load the coastline vtk file ===
vtkPolyDataReader cReader1
   cReader1 SetFileName coastfile_mod_utm.vtk

vtkPolyDataMapper cMapper1
  cMapper1 SetInput [cReader1 GetOutput]

vtkActor cActor1
  cActor1 SetMapper cMapper1
  eval [cActor1 GetProperty] SetColor $navy
  eval [cActor1 GetProperty] SetLineWidth 2

# === load the borders vtk file ===
vtkPolyDataReader cReader2
   cReader2 SetFileName borderfile_mod_utm.vtk

vtkPolyDataMapper cMapper2
  cMapper2 SetInput [cReader2 GetOutput]

vtkActor cActor2
  cActor2 SetMapper cMapper2
  eval [cActor2 GetProperty] SetColor $navy
  eval [cActor2 GetProperty] SetLineWidth 1

#=== color bar ====
# NOTE hrMapper variable here
vtkScalarBarActor scalarBar
    scalarBar SetLookupTable lut
    #scalarBar SetLookupTable [hrMapper GetLookupTable]
    scalarBar SetTitle "P-wave (m/s)"
    [scalarBar GetPositionCoordinate] SetCoordinateSystemToNormalizedViewport
    [scalarBar GetPositionCoordinate] SetValue 0.1 0.0
    scalarBar SetOrientationToHorizontal
    scalarBar SetWidth 0.5
    scalarBar SetHeight 0.1
    scalarBar SetPosition 0.2 0.
    scalarBar SetLabelFormat "%-#6.3g"
    [scalarBar GetLabelTextProperty] SetColor 0 0 0
    [scalarBar GetLabelTextProperty] SetFontFamilyToTimes
    [scalarBar GetTitleTextProperty] SetColor 0 0 0
    [scalarBar GetTitleTextProperty] SetFontFamilyToTimes

#=== text title ===
vtkTextActor titleActor
   titleActor SetDisplayPosition 250 60
   titleActor SetInput "Lin model with Harvard model -- Depth = -0.0 km"
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
    renWin SetSize 600 600 
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

vtkCamera cam1
    cam1 SetPosition 385719.8 3823329.2 50000.0 
    cam1 SetFocalPoint 385719.8 3823329.2 0.0 
    cam1 SetViewUp 0 1 0

#=== scene ===============
ren1 SetActiveCamera cam1

# colorbar only works with --> nsActor
# kActor is the volume

ren1 AddActor hrActor
ren1 AddActor cActor1
ren1 AddActor cActor2
ren1 AddActor2D scalarBar
ren1 AddActor2D titleActor

# I don't understand why, but ResetCamera is very important
ren1 ResetCamera 300000 550000 3572000 4075000 -60000 200000
renWin Render

vtkWindowToImageFilter w2i  
  w2i SetInput renWin
vtkPostScriptWriter writer
  writer SetInputConnection [w2i GetOutputPort]
  writer SetFileName "xc_hr_02.ps"
  writer Write

wm withdraw .
