This is working directory of FDPS on Crystal project

Current code is still proof-of-concept and is not fully working yet

If you want to compile the current code, you first need to download
FDPS and compile Fortran version of nbody test program. Then you place
the source files (Makefile as well) to that directory
(presumably $FDPS/sample/fortran/nbody), edit the Makefile
to correct FDPS_LOC (uncomment  the first commented one and coment out
the current active one), and then try

   make fdpscr

