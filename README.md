This is the working directory of FDPS on Crystal project.

The current code is still more like a proof of concept but this sample
at least can be compiled and run without relying on fortran sample code.

If you want to compile the current code, you first need to download
FDPS. Better to test if the fortran sample works or not.  Then you
place the source files in some directory, edit Makefile so that
FDPS_LOC points to reasonable place, and then try

   make fdpscr

# TODO

* enum types shuld also be generated

* calling interfaces to xxxx_0000x functions shuld also be generated

* get_psys_cptr  shuld also be generated

* Makefile needs some cleaning up
