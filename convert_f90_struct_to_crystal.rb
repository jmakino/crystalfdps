#
# convert_crystal_struct_to_f90.rb
#
def print_header
  print <<EOF
module FDPS_vector
lib FDPS
EOF
end
def print_trailer
  print <<EOF
end
end
EOF
end


$type_conversion_table =<<-EOF
integer(kind=c_long_long) Int64
real(kind=c_double)   Float64      
real(kind=c_float)   Float32
type(fdps_f64vec) Cvec_Float64
type(fdps_f32vec) Cvec_F32
type(fdps_f64mat) Cmat_F64
type(fdps_f32mat) Cmat_F32
EOF


$type_conversion_hash = $type_conversion_table.split("\n").map{|s| s.split}.to_h


def print_crystal_member(s)
  a = s.split
  return if a.length==0
  ss = "     "
  if a[0][0] != "!"
    # member variable line
    varname = a[2]
    typename = $type_conversion_hash[a[0]]
    raise "unparsable line: #{s}"    if a[1] != "::" || typename == nil
    ss += " "+varname + " : "+ typename
    a = a[3..(a.length)]
  end
  if a.length >0
    a[0][0]= "#"
    ss += " "+a.join(" ")
  end
  print ss,"\n"
end
        
      
      
    
print_header

while s=gets
  if s =~ /!\*\*\*\* PS::(\S+)/
    print "# "+ s
    struct_name = $1
    s = gets
    lines=[]
    instruct=true
    while instruct
      s=gets
      print "# "+ s
      if s=~/^(\s*)end type/
        instruct =false
      else
        lines.push s.chomp
      end
    end
    print "  struct #{struct_name}\n"
    lines.each{|s| print_crystal_member(s)}
    print "   end\n"
  end
end

print_trailer

