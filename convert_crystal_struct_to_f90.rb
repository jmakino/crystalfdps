#
# convert_crystal_struct_to_f90.rb
#
def print_header
  print <<-EOF
!===============================
!   MODULE: User defined types
!===============================
module user_defined_types
   use, intrinsic :: iso_c_binding
   use fdps_vector
   use fdps_super_particle
   implicit none
EOF
end

def print_trailer
  print <<-EOF
end module user_defined_types
EOF
end

$type_conversion_table =<<-EOF
Int64         integer(kind=c_long_long)
Float64       real(kind=c_double)  
Cvec_Float64  type(fdps_f64vec)
EOF


$type_conversion_hash = $type_conversion_table.split("\n").map{|s| s.split}.to_h

open("crmain.cpp","w"){|f|
       f.print <<EOF
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
EOF
     }
def print_fortran_member(s)
  a = s.split
  ss = "     "
  if a[0][0] != "#"
    # member variable line
    varname = a[0]
    typename = $type_conversion_hash[a[2]]
    raise "unparsable line: #{s}"    if a[1] != ":" || typename == nil
    ss += " "+typename + " :: "+ varname
    a = a[3..(a.length)]
  end
  if a.length >0
    if a[0] == "#!fdps" || a[0] == "#$fdps"
      a[0] = "!$fdps"
    else
      a[0][0]= "!"
    end
    ss += " "+a.join(" ")
  end
  print ss,"\n"
end
        
      
      
    
print_header

while s=gets
  if s =~ /struct (.*) #!fdps (.*)/
#    print s
    struct_name = $1.downcase
    struct_types=$2
    lines=[]
    instruct=true
    while instruct
      s=gets
      if s=~/^(\s*)end(\s*)$/
        instruct =false
      else
        lines.push s.chomp
      end
    end
    print "   type, public, bind(c) :: #{struct_name}  !$fdps #{struct_types}\n"
    lines.each{|s| print_fortran_member(s)}
    print "   end type #{struct_name}\n"
  end
end

print_trailer

