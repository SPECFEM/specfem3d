
To create individual movie frames, use:
plot_one_movie_frame.bash

Then use mkmpeg4.py to create an AVI or DIVX movie from the frames:
./mkmpeg4.py -v -f number_of_frames_per_second *jpg -o movie.avi

mkmpeg4.py uses mencoder (from http://www.mplayerhq.hu) and is free.

For a list of animation tools, see e.g.:
http://www.gfdl.noaa.gov/products/vis/animation/

