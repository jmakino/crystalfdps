/* Standard headers */
#include <iostream>
#include <fstream>
/* FDPS headers */
#include <particle_simulator.hpp> 
/* User-defined headers */
#include "FDPS_Manipulators.h"
extern "C"{
    void crystal_init(void);
    void crmain(int argc, char *argv[]);

}
int main(int argc, char *argv[])
{
   //* Initialize fdps_manip
   FDPS_Manipulators::Initialize(argc,argv);
   crystal_init();
   crmain(argc, argv);
   return 0;
}
