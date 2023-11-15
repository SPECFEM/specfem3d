----------------------------------------------------------------------
readme - xcreate_movie_GMT:
----------------------------------------------------------------------

xcreate_movie_GMT outputs ascii xyz files, convenient for use with GMT.
This code uses significantly less memory than "create_movie_shakemap_AVS_DX_GMT.f90"
and is therefore useful for high resolution runs.

NOTE: This is an ancient version that still uses the old read/compute_parameter_file subroutines.

      To compile properly, you need to copy over the right subroutines you used
      for your runs and modify them accordingly in the main program.

