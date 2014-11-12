
Subject: inversion package - cleaned up talking points
From: Ryan Modrak 
Date: Nov 2014
To: komatitsch, peterda, Jeroen Tromp, Matthieu Lefebvre, Yanhua Yuan, Emanuele, Federica, + others

Postprocessing utilities

> 1. We agree with the idea, expressed by Daniel last month, that it would be good to include Fortran utilities for smoothing, combining, and preconditioning kernels within each of the specfem repositories.  Implementation of the model update scheme and line search, on the other hand, is probably better accomplished within a separate python repository.
>

i agree, model updates are basically operations on vectors. thus it could be independent from the solver packages. however, if you want to visualize model updates, gradients and/or model differences, it is specific to each package again (mesh topology dependent). for the line search, i think if you use a subset of events and rely on misfit measurements for those, it becomes again package dependent to run the simulations.

another problem has been noticed by emanuele and federica with SPECFEM3D_Cartesian. their mesh creation took several hours or even days, thus they wrote the model_update tool which directly modifies the databases for the solver to avoid re-running the generate_databases executable. this is very package specific, any similar ideas for inputting a new model into the solver databases would be again very specific for each of the 2d/3d packages.

anyway, having a general library for gradient optimizations would be great to have, any effort in this direction is worth trying. there are similar efforts here at ETH, munich and probably also in Dimitri’s group. i would say, it is still worth trying to create an easy-to-use inversion package. we might end up having a few different ones to try out, which is still okay as there is not just one way to do things.


> 2. For the most part, postprocessing utilities could be made to share a common command line interface.  For example, sum_kernels could be made like smooth_vol_data by having names of kernels passed from the command line rather than hardwired in the source code.

yes, sum_kernels takes a list of event names and then uses a set of hardwired kernel names (e.g. bulk_c and bulk_betav, ..). smooth_vol_data is using a single kernel name as input, but no event name list.

i would have liked if the smoothing would also use a set of kernel names rather than only a single kernel. smoothing can become fairly expensive, e.g. run for several hours. the most expensive calculation therein is the construction of the gaussian kernel, which is purely geometry dependent, i.e. depends on the mesh topology and gaussian widths. smoothing a single gradient (bulk_c or bulk_betav or ..) and then calling it several times for all other kernels involved is thus doing most of the computation redundantly. it would be much more efficient to smooth a set of gradients (bulk_c and bulk_betav and …). i would thus prefer the latter.

the only ugly part is that the set of kernel names is hard coded at the moment and steered by a constant parameter in setup/constants_tomography.h, thus it needs to be set before compilation. this could be made more user-friendly, e.g. use a new Par_file_tomo in which one specifies the set of kernel names used for the inversions and also the type of gradient scheme (steepest descent, cg, …), or we extend the command-line arguments to avoid another parameter file.

at the moment, the tomography tools with their corresponding different comand-line arguments as provided are in use by ebru for the global inversions, by me in zurich, and by at least another group at jamstec ( and i’m not sure about min and hejun). changing these input formats should be done within a new package version. thus, i would prefer validating the packages now, tag them and then modify and change these input formats for a next release.


> 3. It is easy to make changes to the mesher or solver and forget that corresponding changes need to be made to the postprocessing utilities.  In the past, several of the programs in src/shared have been broken this way.   There may be precautions we can take against this.

in the past, e.g. the tomography tools have been in a separate package. i think that’s the main reason we saw these diverge. now, all tomography tools which rely on specific mesh topologies (local or global) are within the same package in the same src/ directory. this simplifies the coding. all we need to do is to setup robust buildbot tests to avoid breaking the tools when modifying the meshers.

also, you will have to consider that it is easier to break post-processing tools when they are put in a separate repository, like the idea with the gradient optimization tools. well, unless the interface between the two is well defined and more or less independent of the mesh then it should still work.



> 4. Currently, utilities for smoothing, combining, and preconditioning kernels are provided for Specfem3d only.  Before beginning work on versions for Specfem2d and 3d_globe, additional refactoring and testing of the Specfem3d utilities is probably needed.

i added summing, summing+preconditioning, smoothing, and model updates to the 3d packages, thus they are all available for both, SPECFEM3D_Cartesian and SPECFEM3D_GLOBE. i think they are only not available in the 2D package.



Version tagging

> 5. A tagged release system, in which commits to the master branch are periodically assigned major or minor version numbers, could help with development of the inversion package by providing a set of reference points for testing and collaboration

tagging would be great, but together with validating the code. unfortunately, nobody has had time to validate the full packages. thus we are a bit behind of tagging the master branch.

this is related to the buildbot system. we just need to add validation, regression and unit testing to buildbot as well. once done, the validation in future becomes easier to handle and tagging an easy and valuable thing to do.


> 6. Besides signifying important releases, tags could be used to provide a clear indication of when user action is necessary: for example, we might say that an update that breaks the parameter file format must be assigned a major version number, rather than minor version number, so that users would know that new parameter files are needed.

yes, tagging should follow some standard software coding conventions. bug fixes increase only the minor, changing Par_file or constants.h files which require a re-configuration should increase a major number. new features and major code additions also increase the major.




Solver standardization

> 7. There are a lot of inconsistencies between Specfem2d, 3d, and 3d_globe involving filename conventions, model and data formats, and procedures for reading and writing adjoint traces. This might be a good time to address this issue.  By doing so now, we might avoid code duplication in the preprocessing utilities as well as in the inversion toolkit.

agree, we should constantly put an effort into increasing the consistency between the packages. for me this doubles a lot of my work. and the 2D package has been somewhat out of my scope. so, please just go ahead whenever you encounter such...



> 8. Without proper communication and documentation, file name and format changes of the kind proposed below could be disruptive.  If we go ahead with any changes, we should provide advance warning through the listserv and consider other steps to minimize disruption.

agree, listserv is a good place to discuss major changes.


> 9 In specfem2d, the read_external_model routine could be modified so that models can be read in and kernels can be written out using the same format.

nice idea.



> 10. The ‘proc*.bin’ binary format used by Specfem2d is different than the one used by Specfem3d.  In the former all material parameters share a single file, while in the latter each material parameter has its own file.  For consistency between packages, we could perhaps change the way Specfem2d writes binary files.

in SPECFEM3D, material parameters are only output to single files when the SAVE_MESH_FILES option was used. these single files are only used for post-processing. the solver still uses the material properties which are contained in the database binary by default, keep this in mind.


> 11. In specfem3d, the model_tomography routine appears to be broken.  The model_gll routine supports only isotropic models.  There may be a way to reduce the degree to which model parameters are hardwired to make it easier to deal with alternate parameterizations, for example, involving Thosmen parameters.

for the model_tomography routine, please check again, i fixed the tomography routine a while ago (see e.g. in commit c074b0db0ee4a8bf0f9ef8b649b92f
0ce6a390c6). if it still fails, please let me know.

the model_gll could easily be extended to use an extension like gll_iso, gll_aniso, etc. to make it more specific. please feel free to add and contribute more models.


12. When writing data in Seismic Unix format, Specfem3d creates multiple files--one for each processor.  This makes the data format dependent on the mesh, requiring data be rewritten whenever the mesh changes and making comparison of observations and synthetics difficult.  To avoid these problems, we should probably write all data to a single file.

good idea. it has been a problem in the past to collect everything, as noticed by yang. there can be lots of receivers (>10,000) in exploration setups. collecting them into a single file by the master process might take a while, one would have to re-evaluate this.


> 13. For writing seismic data, perhaps we might adopt something like the following convention
>
> Specfem2d and Specfem3d, Seismic Unix format
> Ux_$event_name$.su – x component
> Uy_$event_name$.su – y component
> Uz_$event_name$.su - z component
> Up_$event_name$.su - pressure
>
> Specfem2d and Specfem3d, single precision binary format
> Ux_$event_name$.bin – x component
> Uy_$event_name$.bin – y component
> Uz_$event_name$.bin - z component
> Up_$event_name$.bin - pressure
>
>
> Specfem3d_globe, single precision binary format
> Ue_$event_name$.bin – East-West component
> Un_$event_name$.bin – North-South component
> Uz_$event_name$.bin - verticle component
> Up_$event_name$.bin - pressure


good suggestion, this could improve things. i think practically, one would have to see if using the event name in the file name makes scripting a bit tedious, as one would have to change this for every different event. maybe using a generic name would make things easier. anyway, feel free to change this as suggested.


Addendum

14. In preparing the original talking points, I overlooked some revisions made by Daniel a few days prior.  The changes, involving postprocessing utilities for Specfem3d and 3d_globe, included some much needed refactoring and improved consistency between packages.  Many thanks for these extremely useful contributions.



15. For tomorrow's discussion, it might be useful to make the following distinction between postprocessing and optimization routines.  Postprocessing routines take kernels as input and return the gradient direction as output, possibly after smoothing and preconditioning.  Optimization routines take the model and gradient direction as input and return a model update direction and step length as output.


16. Using the above terminology and building on Daniel's point by point responses, postprocessing routines depend on the mesh and are tied to a particular solver repository, while optimization routines involve vector operations and can be implemented through a standalone optimization library or inversion package.  In our experience, a simple utility function that converts between the vector data type used in the optimization routines and the "dictionary" data type used in the postprocessing routines can be used to make the optimization routines entirely independent of the mesh, thus getting around difficulties mentioned in Daniel's response to the first talking point.  The end result, we found, is that it is possible to have a standalone inversion package that is implemented mostly in python, that avoids code duplication across multiple repositories, and that through object oriented programming techniques provides flexibility and extensibility to meet a wide variety of user needs.


17. Over the past few months, we have become excited about our inversion package, as two current grad students, one current postdoc, and one previous undergrad, all working on different projects, have had good experiences using it.  A lot more documentation needs to be written before the package is ready for wider release, but if anyone is interested, one or two specific examples could be used for an early preview.


18. It seems everyone agrees on the need for more routine tagging, but obstacles remain in the way of any sort of long term solution.  Perhaps it would make sense to make a single validated, tagged release of both Specfem 3d and 3d_globe a priority, even while long term validation and testing issues are still being sorted out.


19. With the addition of src/tomography utilities to both Specfem3d and 3g_globe, the work of standardizing the packages appears well underway.  Perhaps tomorrow we could go through the talking points to see if there are any remaining inconsistencies between Specfem3d and 3d_globe worth addressing, and to see if any changes to file names or data formats might be justified.

20. While it seems not very practical to overhaul Specfem2d, limited changes to model/ kernel/trace names and formats could be made fairly easily and have an overall positive impact.  Perhaps tomorrow we could decide on what name and format changes to make to Specfem2d.  Myself and perhaps Yanhua could work on making these modifications, being sure to touch base with Etienne, Alexis and other users.



From Carl

Flexwin appears as a stand-along CIG code on github and also now copied within the SPECFEM3D package. This should probably be fixed as part of any overhaul.


===============================================

Proposed action items:
----------------------

Based on Jeroen's suggestions, below are some action items from Wednesday's Skype call.  If anyone can suggest any improvements, in particular if anyone was assigned a task that is not a good match for them, please feel free to email me with changes.  From the result, I'll go ahead and create a set of issues in a new github repository created just for this purpose.


Best regards,

Ryan


 

1. Decide how to organize repositories

 

            Matthieu, Ryan, Dimitri, Daniel

 

 

A. For each of following software categories, a github repository already exists, or may exist in the future:

 

solver package – specfem2d, specfem3d, specfem3d_globe

preprocessing toolkit – obspy, asdf tools

inversion package – seisflows

optimization toolkit – ?

 

B. Rather than create separate repositories to hold postprocessing routines, we agreed it makes sense to maintain a set of postprocessing routines within each existing solver repository.

 

C. In developing optimization routines, the main the challenge lies not in implementing the algorithms, but in integrating with all the other software elements needed for an inversion.  As it happens, dozens of optimization toolkits already exist (see http://en.wikipedia.org/wiki/List_of_optimization_software).  Rather than creating just one more toolkit, it's worth trying to keep sight a more integrated solution, like a workflow manager or inversion package.

 

D. To repeat, in talking about optimization and postprocessing routines, we made the following distinction:  Postprocessing routines take kernels as input and return the gradient direction as output, possibly after smoothing and preconditioning.  Optimization routines take the model and gradient direction as input and return a model update direction and step length as output.)

 

 

2. Figure out how to deal with dependencies

 

            Mattheiu

 

Possible approaches include (A) manual installation guided by website instructions, (B) automated installation via package manager, (C) git submodules, (D) git subtrees, (E) linked libraries (F) other approaches not discussed during Skype call.  Of all these approaches, we need to figure out what makes the most sense for us.

 

 

3. Create a “getting started” webpage for adjoint tomography, with tutorials, faqs, documentation, and various other links

 

            Matthieu

 

           

4. Create buildbot scripts from Daniel’s tests

 

David (?)

 

 

5. Make it easier for developers to add or modify buildbot scripts

 

David

 

 

6. Document major and minor version tagging conventions

 

            Dimitri (?)

Answer from Dimitri: 
I agree with your suggestions in a previous email, thus we are all set and can consider this point as closed.


 

7. Release a new tagged version of Specfem3d and 3d_globe

 

Dimitri, Daniel, Matthieu, David

 
Answer from Dimitri: We can release Specfem3d now; for 3d_globe the only thing we need to clarify is AK135, I have emailed Elliott and cc'ed you. 

 

8.  Think about ways to eliminate code duplication between Specfem3d and 3d_globe

           

            Dimitri, Daniel, Matthieu

 

 Answer from Dimitri: I will do that in Jan 2015. 

 

9. Continue to improve postprocessing utilities for Specfem3d and 3d_globe

 

            ?

 

 

10. Create buildbot tests for Specfem3d and 3d_globe postprocessing utilities

 

            ?

 

 

11. Develop postprocessing utilities for Specfem2d

 

            ?

 
 Answer from Dimitri: Not sure we need 11; we should instead use the 3D tools from 9/ and just apply them to the 2D code, i.e. make the 2D code conform to the 3D standard and conventions; and then we are all set. It seems you will do that in 12 anyway, thus once you are done with 12 and 9 I think we will automatically be done with 11 as well. 

 

12. Change trace/model/kernel names and formats for Specfem2d

 

            Ryan, Yanhua

 

 

13.  Standardize the way adjoint sources are read

 

            Wenjie

 

 

14.  Refactor Seismic Unix reader/writer in Specfem2d, Specfem3d

 

            Ryan (?)
            
