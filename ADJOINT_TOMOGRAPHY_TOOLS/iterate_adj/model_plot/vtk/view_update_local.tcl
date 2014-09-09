#
# This program automates the cross-section cutting for Banana-Doughnut kernels
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
    kReader1 SetFileName /net/sierra/raid1/carltape/socal/socal_3D/RUNS/MODELS/dm15/beta_window/dm_mu_kernel_smooth_p090.vtu

# vtkXMLUnstructuredGridReader kReader2
#     kReader2 SetFileName vtu_files/lin_model/vp_low_2.vtu

# vtkXMLUnstructuredGridReader kReader3
#     kReader3 SetFileName vtu_files/lin_model/vp_low_3.vtu

## When the files are merged, everything crashes.
#
#vtkAppendFilter kReader
#  kReader AddInput [kReader1 GetOutput]
#  kReader AddInput [kReader2 GetOutput]
#  kReader AddInput [kReader3 GetOutput]
#
## data mappers
#vtkDataSetMapper kMapper
#    kMapper SetInput [kReader GetOutput]
#    kMapper SetScalarRange 4687.47 5729.13
#
## data Actor
#vtkActor kActor
#    kActor SetMapper kMapper

# data mappers
vtkDataSetMapper kMapper1
    kMapper1 SetInput [kReader1 GetOutput]
    kMapper1 SetScalarRange -8.0000e-01 8.0000e-01
# vtkDataSetMapper kMapper2
#     kMapper2 SetInput [kReader2 GetOutput]
#     kMapper2 SetScalarRange 4687.47 5729.13
# vtkDataSetMapper kMapper3
#     kMapper3 SetInput [kReader3 GetOutput]
#     kMapper3 SetScalarRange 4687.47 5729.13

# data actors
vtkActor kActor1
    kActor1 SetMapper kMapper1
# vtkActor kActor2
#     kActor2 SetMapper kMapper2
# vtkActor kActor3
#     kActor3 SetMapper kMapper3

# === generate plane cut ====
vtkPlane hrPlane1
    hrPlane1 SetOrigin 385719.8 3823329.2 -40000.0 
    hrPlane1 SetNormal 0 0 1 
vtkCutter hrCut1
    hrCut1 SetInput [kReader1 GetOutput]
    hrCut1 SetCutFunction hrPlane1
vtkDataSetMapper hrMapper1
   hrMapper1 SetInput [hrCut1 GetOutput]
   hrMapper1 InterpolateScalarsBeforeMappingOn
   hrMapper1 SetScalarRange -8.0000e-01 8.0000e-01
   hrMapper1 SetLookupTable lut
vtkActor hrActor1
    hrActor1 SetMapper hrMapper1

# vtkPlane hrPlane2
#     hrPlane2 SetOrigin 385719.8 3823329.2 0.0 
#     hrPlane2 SetNormal 0 0 1 
# vtkCutter hrCut2
#     hrCut2 SetInput [kReader2 GetOutput]
#     hrCut2 SetCutFunction hrPlane2
# vtkDataSetMapper hrMapper2
#    hrMapper2 SetInput [hrCut2 GetOutput]
#    hrMapper2 InterpolateScalarsBeforeMappingOn
#    hrMapper2 SetScalarRange 4687.47 5729.13
#    hrMapper2 SetLookupTable lut
# vtkActor hrActor2
#     hrActor2 SetMapper hrMapper2

# vtkPlane hrPlane3
#     hrPlane3 SetOrigin 385719.8 3823329.2 0.0 
#     hrPlane3 SetNormal 0 0 1 
# vtkCutter hrCut3
#     hrCut3 SetInput [kReader3 GetOutput]
#     hrCut3 SetCutFunction hrPlane3
# vtkDataSetMapper hrMapper3
#    hrMapper3 SetInput [hrCut3 GetOutput]
#    hrMapper3 InterpolateScalarsBeforeMappingOn
#    hrMapper3 SetScalarRange 4687.47 5729.13
#    hrMapper3 SetLookupTable lut
# vtkActor hrActor3
#     hrActor3 SetMapper hrMapper3

# === source (sphere->mapper->actor) ====
# vtkSphereSource sourceSphere
#     sourceSphere SetCenter  427523.33   3752613.39 0.0
#     sourceSphere SetRadius 15000.0
#     sourceSphere SetThetaResolution 20
#     sourceSphere SetPhiResolution 20

# vtkPolyDataMapper sourceMapper
#     sourceMapper SetInput [sourceSphere GetOutput]

# vtkActor sourceActor
#     sourceActor SetMapper sourceMapper
#     eval [sourceActor GetProperty] SetColor $hot_pink
#     [sourceActor GetProperty] SetSpecularColor 1 1 1
#     [sourceActor GetProperty] SetSpecular 0.3
#     [sourceActor GetProperty] SetSpecularPower 20
#     [sourceActor GetProperty] SetAmbient 0.2
#     [sourceActor GetProperty] SetDiffuse 0.8

# # === source from file ====

# vtkSphereSource sourceSphere
#      sourceSphere SetRadius 15000.0
#      sourceSphere SetThetaResolution 20
#      sourceSphere SetPhiResolution 20

# vtkPolyDataReader sourceReader
#     sourceReader SetFileName source.vtk

# vtkGlyph3D sourceGlyph
#     sourceGlyph SetInputConnection [sourceReader GetOutputPort]
#     sourceGlyph SetSource [sourceSphere GetOutput]

# vtkPolyDataMapper sourceMapper
#    sourceMapper SetInput [sourceGlyph GetOutput]

# vtkActor sourceActor
#      sourceActor SetMapper sourceMapper
#      eval [sourceActor GetProperty] SetColor $hot_pink
#      [sourceActor GetProperty] SetSpecularColor 1 1 1
#      [sourceActor GetProperty] SetSpecular 0.3
#      [sourceActor GetProperty] SetSpecularPower 20
#      [sourceActor GetProperty] SetAmbient 0.2
#      [sourceActor GetProperty] SetDiffuse 0.8

# # === receivers from file ===

# vtkSphereSource receiverSphere
#      receiverSphere SetRadius 5000.0
#      receiverSphere SetThetaResolution 20
#      receiverSphere SetPhiResolution 20

# vtkPolyDataReader receiverReader
#     receiverReader SetFileName receiver.vtk

# vtkGlyph3D receiverGlyph
#     receiverGlyph SetInputConnection [receiverReader GetOutputPort]
#     receiverGlyph SetSource [receiverSphere GetOutput]

# vtkPolyDataMapper receiverMapper
#    receiverMapper SetInput [receiverGlyph GetOutput]
#    #receiverMapper SetInputConnection [receiverGlyph GetOutputPort]

# vtkActor receiverActor
#      receiverActor SetMapper receiverMapper
#      eval [receiverActor GetProperty] SetColor $grey
#      [receiverActor GetProperty] SetSpecularColor 1 1 1
#      [receiverActor GetProperty] SetSpecular 0.3
#      [receiverActor GetProperty] SetSpecularPower 20
#      [receiverActor GetProperty] SetAmbient 0.2
#      [receiverActor GetProperty] SetDiffuse 0.8

# === load the coastline vtk file ===
vtkPolyDataReader cReader1
   cReader1 SetFileName vtk_files/coastfile_mod_utm_air.vtk

vtkPolyDataMapper cMapper1
  cMapper1 SetInput [cReader1 GetOutput]

vtkActor cActor1
  cActor1 SetMapper cMapper1
  eval [cActor1 GetProperty] SetColor $navy
  eval [cActor1 GetProperty] SetLineWidth 2

# === load the borders vtk file ===
vtkPolyDataReader cReader2
   cReader2 SetFileName vtk_files/borderfile_mod_utm_air.vtk

vtkPolyDataMapper cMapper2
  cMapper2 SetInput [cReader2 GetOutput]

vtkActor cActor2
  cActor2 SetMapper cMapper2
  eval [cActor2 GetProperty] SetColor $navy
  eval [cActor2 GetProperty] SetLineWidth 1

# === load the plate boundary vtk file ===
vtkPolyDataReader cReader3
   cReader3 SetFileName vtk_files/NA_PA_boundary_utm_air.vtk

vtkPolyDataMapper cMapper3
  cMapper3 SetInput [cReader3 GetOutput]

vtkActor cActor3
  cActor3 SetMapper cMapper3
  eval [cActor3 GetProperty] SetColor $navy
  eval [cActor3 GetProperty] SetLineWidth 2

# === load the faults vtk file ===
vtkPolyDataReader cReader4
   cReader4 SetFileName vtk_files/jennings_more_utm_air.vtk

vtkPolyDataMapper cMapper4
  cMapper4 SetInput [cReader4 GetOutput]

vtkActor cActor4
  cActor4 SetMapper cMapper4
  eval [cActor4 GetProperty] SetColor $navy
  eval [cActor4 GetProperty] SetLineWidth 1.5

#=== color bar ====
# NOTE hrMapper variable here
vtkScalarBarActor scalarBar
    scalarBar SetLookupTable lut
    #scalarBar SetLookupTable [hrMapper GetLookupTable]
scalarBar SetTitle " "
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
   titleActor SetDisplayPosition 280 60
   titleActor SetInput "Subspace update for SHEAR-WAVE-SPEED -- model 90 -- Cut at z = 40.0 km"
set tprop [titleActor GetTextProperty]
   $tprop SetJustificationToCentered
   $tprop SetColor 0 0 1
   $tprop SetFontFamilyToTimes
   $tprop SetFontSize 15

#=== render window and camera positions ===
vtkRenderer ren1
    ren1 SetBackground 1 1 1
vtkRenderWindow renWin
    renWin AddRenderer ren1
    renWin SetSize 600 600 
vtkRenderWindowInteractor iren
    iren SetRenderWindow renWin

vtkCamera cam1
    cam1 SetPosition 385719.8 3823329.2 10000.0 
    cam1 SetFocalPoint 385719.8 3823329.2 -40000.0 
    cam1 SetViewUp 0 1 0

#=== scene ===============
ren1 SetActiveCamera cam1

# colorbar only works with --> nsActor
# kActor is the volume

#ren1 AddActor hrActor
ren1 AddActor hrActor1
#ren1 AddActor hrActor2
#ren1 AddActor hrActor3
ren1 AddActor cActor1
ren1 AddActor cActor2
#ren1 AddActor cActor3
ren1 AddActor cActor4
#ren1 AddActor sourceActor
#ren1 AddActor receiverActor
ren1 AddActor2D scalarBar
ren1 AddActor2D titleActor

# I don't understand why, but ResetCamera is very important
ren1 ResetCamera 300000 550000 3572000 4075000 -60000 200000
renWin Render

vtkWindowToImageFilter w2i  
  w2i SetInput renWin
vtkPostScriptWriter writer
  writer SetInputConnection [w2i GetOutputPort]
  writer SetFileName "dm_mu_kernel_smooth_p090_13.ps"
  writer Write

wm withdraw .

exit

#set cam2 [ren1 GetActiveCamera]
#puts stdout [ $cam2 GetPosition]
#puts stdout [$cam2 GetFocalPoint]
#puts stdout [$cam2 GetViewUp]

