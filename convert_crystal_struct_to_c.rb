#
# convert_crystal_struct_to_c.rb
#
def print_header
  print <<-EOF
#pragma once
/* Standard headers */
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"
EOF
end

$type_conversion_table =<<-EOF
Int64         long long
Float64       double
Cvec_Float64  fdps_f64vec
Cvec_Float32  fdps_f32vec
EOF


$type_conversion_hash = $type_conversion_table.split("\n").map{|s|
  a= s.split
  [a[0], a[1..(a.size-1)].join(" ")]}.to_h

open("cppmain.cpp","w"){|f|
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
    ss += " "+typename + " "+ varname+";"
    a = a[3..(a.length)]
  end
  if a.length >0
    if a[0] == "#!fdps" || a[0] == "#$fdps"
      a[0] = "//$fdps"
    else
      a[0][0]= "//"
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
    print "typedef struct #{struct_name}{  //$fdps #{struct_types}\n"
    lines.each{|s| print_fortran_member(s)}
    print "} #{struct_name.capitalize};\n"
  end
end
