# crystalfdps

A working  sample code using FDPS, writen entirely in 
Crystal language.

To cmpile, after you download and place the files (either from
zip file or by git clone), 

   shards install
   make fdpscr

should create an executable. Try

  ./fdpscr -h

to get the list of command-line options.

Edit Makefile to use OpenMP (enabled by default) and MPI

