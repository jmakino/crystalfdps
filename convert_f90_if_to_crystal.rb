#
# convert_f90_if_to_crystal.rb
#

$type_conversion_table =<<-EOF
void          Void
int           Int32
bool          Bool
PS::S64       Int64
PS::S32       Int32
PS::F32       Float32
PS::F64       Float64
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
  if s =~/^([a-z]\w*) fdps_(\S*)\((.*)\)/
    a.push $1
    a.push $2
    a += $3.split(",")
  end
  a
end

def preprocess_funarg(a)
  b=[]
  s= a.shift
#  print "Funarg[0]=",s,"\n"
  while s
    if s =~/([a-z]\w*) \(\*/
#      print "Funarg complex type"
      ss = s
      s= a.shift
      while s !~ /\)/ 
        ss += s
        s= a.shift
      end
      ss += s
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
  if s.index("(*pfunc_comp)")
    s =~ /pfunc_comp\)\(const struct (\S+) (.*)\n(\s*) const struct (\S+)/
    args = [$1, $4].map{|s| "FDPS::"+ s.tocap}
    s = " pfunc_comp : Pointer(#{args[0]}), Pointer(#{args[1]}) -> Bool"
  elsif s.index("(*pfunc_ep")
    if  s =~/(pfunc_ep\S+)\)\(struct (\S+) .*\n\s* int ,\n\s* struct (\S+) .*\n\s*int.*\n\s* struct (\S+)/
    args=[$1,$2,$3,$4]
    elsif s =~/(pfunc_ep\S+)\)\(struct (\S+) .*\n\s* int ,\n\s* PS::(\S+) .*\n\s*int.*\n\s* struct (\S+)/
      args=[$1,$2,$3,$4]
    end
    funname = args.shift
    args = args.map{|s| "FDPS::"+ s.tocap}
    s = <<EOF
         #{funname}:
                         Pointer(#{args[0]}), Int32,
                         Pointer(#{args[1]}), Int32,
                         Pointer(#{args[2]}) -> Void
EOF
  end
  s.chomp
end

def convert_to_crystal_arg_def(s)
  result=""
  s.gsub!(/\*(\s+)/," *")
  if s =~ /\n/
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
  if s =~  /^([a-z]\w*) fdps_(\S*)\(/
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
