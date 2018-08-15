#==============================
#   MODULE: FDPS vector
#==============================

macro use_fdps_vector(ndim,etype)
module FDPS_vector
lib FDPS
  struct Cvec_{{etype}}
    x : {{etype}}
    y : {{etype}}
    {% if ndim == 3 %} z : {{etype}} {% end %}
  end
  struct Testp
         pos : FDPS::Cvec_{{etype}}
  end
end
lib FDPS
  struct Testp2
         pos : FDPS::Cvec_{{etype}}
  end
end 
  class Vec_{{etype}}
    property x, y {% if ndim == 3 %},  z {% end %}
    def to_unsafe
      cvec = FDPS::Cvec_{{etype}}.new
      cvec.x = @x
      cvec.y = @y
       {% if ndim == 3 %} cvec.z = @z {% end %}
      cvec
    end
    def initialize(@x : {{etype}}, @y : {{etype}}     {% if ndim == 3 %} , @z : {{etype}}{% end %})
    end
    def initialize()
       @x ={{etype}}.new(0)
       @y ={{etype}}.new(0)
       {% if ndim == 3 %}   @z ={{etype}}.new(0) {% end %}
    end
    def initialize(a : {{etype}} )
       @x =a
       @y =a
      {% if ndim == 3 %}   @z =a {% end %}
    end
    def initialize(a : Array({{etype}}) )
       @x =a[0]
       @y =a[1]
      {% if ndim == 3 %}   @z =a[2] {% end %}
    end
    def initialize(a : FDPS::Cvec_{{etype}})
       @x =a.x
       @y =a.y
      {% if ndim == 3 %}   @z =a.z {% end %}
    end
    def +(other : Vec_{{etype}})
       Vec_{{etype}}.new(@x+ other.x,
                       @y+other.y {% if ndim == 3 %} , @z+other.z{% end %})
    end
    def +
       self
    end
    def -(other : Vec_{{etype}})
       Vec_{{etype}}.new(@x- other.x,
                       @y-other.y {% if ndim == 3 %} , @z-other.z{% end %})
    end
    def -
       Vec_{{etype}}.new(-@x, -@y {% if ndim == 3 %} , -@z{% end %})
    end
    def *(other : Vec_{{etype}}) 
       @x*other.x + @y*other.y {% if ndim == 3 %} + @z*other.z{% end %}
    end
    def *(other : {{etype}}) 
       Vec_{{etype}}.new(@x* other,
                       @y*other {% if ndim == 3 %} , @z*other{% end %})
    end
    def /(other : {{etype}})
       scale = {{etype}}.new(1)/other
       self*scale
    end

  end #class
end #module
end #macro


# use_fdps_vector(3,Float64)
# include FDPS_vector

# module  FDPS_vector
# lib FDPS
#   struct Testp3
#          pos : FDPS::Cvec_Float64
#   end
# end 
# end

#a = Vec_Float64.new(1,2,4)
#particle = FDPS::Testp2.new
#p particle
#c : FDPS::Cvec_Float64 = a.to_unsafe
#p a
#p a.to_unsafe
#c =a.to_unsafe
#p a+Vec_Float64.new(c)
#p a+c
# b = FDPS_vec(Float64).new
# p b
# c = FDPS_vec(Float64).new(32)
# p c+a
# p (a+b+c)/10

