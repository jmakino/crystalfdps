/* Standard headers */
#include <iostream>
#include <fstream>
/* FDPS headers */
#include <particle_simulator.hpp> 
/* User-defined headers */
#include "FDPS_Manipulators.h"
extern "C"{
void crystal_init(void);
void crmain(void);

}
int main(int argc, char *argv[])
{
   
   //* Initialize fdps_manip
   FDPS_Manipulators::Initialize(argc,argv);
   //* Call Fortran main subroutine
   //   f_main_();
   crystal_init();
   crmain();
   

   return 0;

}
