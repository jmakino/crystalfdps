#
# convert_f90_if_to_crystal.rb
#

$type_conversion_table =<<-EOF
void          Void
int           Int32
bool          Bool
PS::S64       Int64
PS::S32       Int32
PS::U32       UInt32
PS::U64       UInt64
PS::F32       Float32
PS::F64       Float64
PS::F64vec    Cvec_Float64
PS::F32vec    Cvec_Float32
double        Float64
PS::TimeProfile TimeProfile
PS::F32vec Cvec_F32
char          UInt8
size_t        LibC::SizeT
full_particle  Full_particle
float         Float32
EOF


$type_conversion_hash = $type_conversion_table.split("\n").map{|s| s.split}.to_h

$functions_to_skip=<<-EOF
get_psys_cptr
fdps_calc_force_all
fdps_calc_force_making_tree
fdps_calc_force_and_write_back
fdps_sort_particle
fdps_calc_force_all_and_write_back
EOF
def convert_functype(s)
  ss = $type_conversion_hash[s]
  if ss == nil
    raise "unknown functype #{s}"
  end
  ss
end

def preprocess_single_line_def(s)
  a=[]
  if s =~/^([a-z,A-Z]\w*\:*\w*) fdps_(\S*)\((.*)\)/
    a.push $1
    a.push $2
    a += $3.split(",")
  end
  a
end

def preprocess_funarg(a)
  b=[]
#  print "Funargs=",a,"\n"
  s= a.shift
  while s
#    print "Funarg=",s,"\n"
    if s =~/([a-z,A-Z]\w*\:*\w*) \(\*/
#      print "Funarg complex type"
      ss = s
      s= a.shift
      if s
        while s !~ /\)/ 
          ss += s
          s= a.shift
        end
        ss += s
      end
      s=ss
    end
    b.push s.chomp
    s = a.shift
  end
  b
end
def preprocess_multi_line_def(a)
  b=[]
  s = a.shift
  if s =~/^([a-z]\w*) fdps_(\S*)\((.*)/
    b.push $1
    b.push $2
    a.unshift $3
  else
    raise "unparsable line: #{s}"
  end
  a=preprocess_funarg(a)
  b+a
  
end

class String
  def tocap
    s=self
    s=self.capitalize if self[0] =~/[a-z]/
    s
  end
end      
    
def convert_function_pointers(s)
  if s=~/bool \(\*pfunc_comp\)/
    s = " pfunc_comp : Pointer(Full_particle), Pointer(Full_particle) -> Bool"
  elsif s =~/(pfunc_ep\S+)\)\(void/
    funname = $1
    s = <<EOF
         #{funname}:
                         Pointer(Full_particle), Int32,
                         Pointer(Full_particle), Int32,
                         Pointer(Full_particle) -> Void
EOF
  end
  s.chomp
end

def convert_to_crystal_arg_def(s)
  result=""
  s.gsub!(/\*(\s+)/," *")
  if s =~ /pfunc/
    result= convert_function_pointers(s)
  else
    a=s.split(")")[0].split(",")[0].split
    a.shift if a[0] == "const"
    varname = a.pop
    a.pop if a[a.length-1]== "const"
    typename = $type_conversion_hash[a.join(" ")]
    if typename  == nil
      a.shift if a[0] == "enum"
      typename = a.join(" ")
    end
    if varname[0]=="*"
      if varname[1]=="*"
        varname = varname[1,varname.length-1]
      end
      varname = varname[1,varname.length-1]
      typename = "Pointer(#{typename})"
    end
    result= "#{varname} : #{typename}"
  end
  result
end

def convert_one_func(a)
  funtype = convert_functype(a[0])
  funname = a[1]
  if $functions_to_skip.index(funname)
    print "# function #{funname} skipped --- should be written manually\n"
    return
  end
  print "#" + a.join("|"), "\n"
  funargs = a[2..(a.length-1)].map{|s|  convert_to_crystal_arg_def(s)}
  print "fun #{funname}=fdps_#{funname}("
  print funargs[0] if funargs.length >0
  funargs.shift
  funargs.each{|s|
    print ",\n           ", s
  }
  print ") : #{funtype}\n"
  
end

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
print_header
while s=gets
  if s =~  /^([a-z,A-Z]\w*\:*\w*) fdps_(\S*)\(/
#    p s
    a=[]
    a.push s
    while ! ( a[a.length-1] =~ /\{/ )
      a.push gets
    end
#    print "func body:\n"
    #    print  a.join
    print "\n"
    a.each{|s|print "# "+s}
    print "\n"
      
    if a.length==1
      a =preprocess_single_line_def(a[0])
    else
      a =preprocess_multi_line_def(a)
    end
    convert_one_func(a)
  end
end
print_trailer
