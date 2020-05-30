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
  struct Cvec_F32
    x : Float32
    y : Float32
    {% if ndim == 3 %} z : Float32 {% end %}
  end
  struct Cmat_F32
    xx : Float32
    yy : Float32
    {% if ndim == 3 %} zz : Float32 {% end %}
    xy : Float32
    {% if ndim == 3 %} xz : Float32
    yz : Float32{% end %}
  end
  struct Cmat_F64
    xx : Float64
    yy : Float64
    {% if ndim == 3 %} zz : Float64 {% end %}
    xy : Float64
    {% if ndim == 3 %} xz : Float64
    yz : Float64{% end %}
  end

end

def add(a : FDPS_vector::FDPS::Cvec_{{etype}}, b : FDPS_vector::FDPS::Cvec_{{etype}})
  cvec = FDPS::Cvec_{{etype}}.new
  cvec.x =  a.x +  b.x
  cvec.y = a.y + b.y
  {% if ndim == 3 %} cvec.z = a.z + b.z {% end %}
  cvec
end    

  def mul(a : FDPS_vector::FDPS::Cvec_{{etype}}, b : {{etype}})
      cvec = FDPS::Cvec_{{etype}}.new
      cvec.x = a.x * b
      cvec.y = a.y * b
      {% if ndim == 3 %} cvec.z = a.z*b {% end %}
      cvec
  end    

# Here, stuct must be used instead of class,
# since struct is allocated on stack (no malloc) and thus far faster.
  struct Vec_{{etype}}
    property x, y {% if ndim == 3 %},  z {% end %}
    def to_a
       {% if ndim == 2 %}  [@x,@y]{% end %}
       {% if ndim == 3 %}  [@x,@y,z]{% end %}
    end
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

