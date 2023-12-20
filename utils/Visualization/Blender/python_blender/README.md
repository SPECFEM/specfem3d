# Blender scripting example

Blender for 3D graphics creation<br>
https://www.blender.org/


## Installation

The python script `plot_with_blender.py` uses Blender's python module. To use the script, we also need some routines from vtk which are not provided in the default python version that Blender internally uses. One possibility to use the systems python frameworks is to set an environment variable `BLENDER_SYSTEM_PATH`. For example, on Mac having python installed through [MacPorts](https://www.macports.org), one can set
```
export BLENDER_SYSTEM_PYTHON='/opt/local/Library/Frameworks/Python.framework/Versions/3.10/'
```
For this to work, the python version must match the internal python version from Blender. In this example, Blender version 3.6 uses a python version 3.10.

Another option is to install vtk into the provided Blender python version. For example, on Mac this can be used:
```
/Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10 -m ensurepip
/Applications/Blender.app/Contents/Resources/3.6/python/bin/python3.10 -m pip install vtk==9.2.6
```
In this latter case to avoid problems with blender loading the matplotlib rendering modules, comment out the lines related to vtkRenderingMatplotlib in the corresponding vtk file. For example, on Mac this would be needed:
```
## correct import problem with vtk in blender
# see: https://devtalk.blender.org/t/python-console-undefined-symbol-py-main-upon-import/21143/6
cd /Applications/Blender.app/Contents/Resources/3.6/python/lib/python3.10/site-packages

sed -i "s:from vtkmodules.vtkRenderingMatplotlib import *:#from vtkmodules.vtkRenderingMatplotlib import *:" ./vtk.py
sed -i "s:from .vtkRenderingMatplotlib import *:#from .vtkRenderingMatplotlib import *:" ./vtkmodules/all.py
```
This should work. You could test it by running blender and importing vtk in its scripting environment, or running blender with a short test script `my_test.py`:
```
#my test script
import sys
for path in sys.path: print(path)
import vtk
print("vtk: ",vtk.__path__)
```
and run blender with:
```
blender -noaudio --background --python my_test.py
```


## Simulation setup

First, run a simulation example, e.g., in `EXAMPLES/applications/simple_model/`, and turn on the shakemap flag in the `DATA/Par_file`:
   ```
   CREATE_SHAKEMAP                 = .true.
   ```

   This will create a file `shakingdata` in the example's `OUTPUT_FILES/` folder.
   Then, run the `xcreate_movie_shakemap_AVS_DX_GMT` binary from the root directory to create an AVS UCD format `.inp` output file
   ```
   ../../../../bin/xcreate_movie_shakemap_AVS_DX_GMT <<<EOF
   2
   -1
   3
   2
   EOF
   ```
   This will plot PGA values. Note that you could also create other visualization output files, such as surface & volume movies (as `.vtk` or `.vtu` files).


### Blender rendering

Once the visualization file (e.g., the `AVS_shaking_map.inp` file above) is created, you can use the corresponding file name as parameter in the python script `plot_with_blender.py`, just type:
 ```
 ./plot_with_blender.py --vtk_file=OUTPUT_FILES/AVS_shaking_map.inp
 ```

This will create a Blender file `out.blend` and a rendered image `out.jpg` with the default scene setup. The Blender file `out.blend` can be used as a starting point to open in Blender and do your vis magic!

For more script options, you can type `./plot_with_blender.py --help`.

Feel free to modify and contribute any improvements to this python script - and have fun with Blender for scientific visualization!
