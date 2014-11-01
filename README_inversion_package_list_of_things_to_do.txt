
Subject: inversion package - proposed talking points
From: Ryan Modrak 
Date: 10/27/2014 02:57 PM
To: komatitsch, peterda, Jeroen Tromp, Matthieu Lefebvre, Yanhua Yuan


Over the past weeks, Jeroen, Matthieu, and myself have been meeting regularly to develop an inversion package capable of interfacing with specfem2d, 3d, and 3d_globe.   We were wondering, would you be interested in joining the discussion?

The idea behind the inversion package is similar to what we were talking about last month on the seismo-dev listserv, except rather than thinking in terms of a subdirectory within the specfem3d repository, we were envisioning a standalone python package.  A list of possible talking points is given below:

 
Postprocessing utilities

 

1. We agree with the idea, expressed by Daniel last month, that it would be good to include Fortran utilities for smoothing, combining, and preconditioning kernels within each of the specfem repositories.  Implementation of the model update scheme and line search, on the other hand, is probably better accomplished within a separate python repository.

 

2. For the most part, postprocessing utilities could be made to share a common command line interface.  For example, sum_kernels could be made like smooth_vol_data by having names of kernels passed from the command line rather than hardwired in the source code.

 

3. It is easy to make changes to the mesher or solver and forget that corresponding changes need to be made to the postprocessing utilities.  In the past, several of the programs in src/shared have been broken this way.   There may be precautions we can take against this.


4. Currently, utilities for smoothing, combining, and preconditioning kernels are provided for Specfem3d only.  Before beginning work on versions for Specfem2d and 3d_globe, additional refactoring and testing of the Specfem3d utilities is probably needed.


4bis (added by Dimitri): parallelizing these tools (very useful for large cases); using MPI in the case of Fortran tools, and pyMPI in the case of Python scripts.


Version tagging

 

5. A tagged release system, in which commits to the master branch are periodically assigned major or minor version numbers, could help with development of the inversion package by providing a set of reference points for testing and collaboration

 

6. Besides signifying important releases, tags could be used to provide a clear indication of when user action is necessary: for example, we might say that an update that breaks the parameter file format must be assigned a major version number, rather than minor version number, so that users would know that new parameter files are needed.



Solver standardization

 

7. There are a lot of inconsistencies between Specfem2d, 3d, and 3d_globe involving filename conventions, model and data formats, and procedures for reading and writing adjoint traces. This might be a good time to address this issue.  By doing so now, we might avoid code duplication in the preprocessing utilities as well as in the inversion toolkit.  

 

8. Without proper communication and documentation, file name and format changes of the kind proposed below could be disruptive.  If we go ahead with any changes, we should provide advance warning through the listserv and consider other steps to minimize disruption.

 

9 In specfem2d, the read_external_model routine could be modified so that models can be read in and kernels can be written out using the same format.

 

10. The ‘proc*.bin’ binary format used by Specfem2d is different than the one used by Specfem3d.  In the former all material parameters share a single file, while in the latter each material parameter has its own file.  For consistency between packages, we could perhaps change the way Specfem2d writes binary files.

 

11. In specfem3d, the model_tomography routine appears to be broken.  The model_gll routine supports only isotropic models.  There may be a way to reduce the degree to which model parameters are hardwired to make it easier to deal with alternate parameterizations, for example, involving Thosmen parameters.

 

12. When writing data in Seismic Unix format, Specfem3d creates multiple files--one for each processor.  This makes the data format dependent on the mesh, requiring data be rewritten whenever the mesh changes and making comparison of observations and synthetics difficult.  To avoid these problems, we should probably write all data to a single file.

 

13. For writing seismic data, perhaps we might adopt something like the following convention

 

Specfem2d and Specfem3d, Seismic Unix format

Ux_$event_name$.su – x component

Uy_$event_name$.su – y component

Uz_$event_name$.su - z component

Up_$event_name$.su - pressure

 

Specfem2d and Specfem3d, single precision binary format

Ux_$event_name$.bin – x component

Uy_$event_name$.bin – y component

Uz_$event_name$.bin - z component

Up_$event_name$.bin - pressure

 

Specfem3d_globe, single precision binary format

Ue_$event_name$.bin – East-West component

Un_$event_name$.bin – North-South component

Uz_$event_name$.bin - verticle component

Up_$event_name$.bin - pressure


