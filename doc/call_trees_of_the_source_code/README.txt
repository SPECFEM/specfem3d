
To visualize the calling tree of SPECFEM3D:

1/ install "doxygen" on your machine, if not already installed (for instance typing " sudo apt-get install doxygen " or similar on Linux machines)

2/ type:

      doxygen Doxyfile_complete_call_tree

   or

      ./run_doxygen.sh

   for a default configuration of the code documentation.

3/ visualize the calling tree using a Web browser; for instance with firefox:

      firefox html/index.html


note:
for generating call tree graphs, you also need graphviz installed (type "sudo apt-get install graphviz" or similar)


