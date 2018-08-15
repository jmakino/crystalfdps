#===============================
#   MODULE: User defined types
#===============================
module FDPS_vector
lib FDPS
   struct Full_particle #!fdps FP,EPI,EPJ,Force
      #!fdps copyFromForce full_particle (pot,pot) (acc,acc)
      #!fdps copyFromFP full_particle (id,id) (mass,mass) (eps,eps) (pos,pos) 
      #!fdps clear id=keep, mass=keep, eps=keep, pos=keep, vel=keep
      id : Int64 
      mass : Float64
      eps : Float64
      pos : FDPS::Cvec_Float64 #!fdps position
      vel : Cvec_Float64 #!fdps velocity
      pot : Float64
      acc : Cvec_Float64 
    end
end
end

#   !**** Interaction function (particle-particle)
include Math
def calc_gravity(ep_i,n_ip,ep_j,n_jp,f)
  STDERR.print "in calc gravity ni=#{n_ip}  nj=#{n_jp} \n"
  n_ip.times{|i|
    pi = (ep_i + i).value
    eps2 = pi.eps*pi.eps
    xi = Vec_Float64.new(pi.pos)
    ai =  Vec_Float64.new(0)
    poti = 0_f64
    n_jp.times{|j|
      pj = (ep_j + j).value
      xj = Vec_Float64.new(pj.pos)
      rij = xi - xj
      r2 = rij*rij+eps2
      rinv = 1_f64/sqrt(r2)
      mrinv = pj.mass*rinv
      mr3inv = mrinv*rinv*rinv
      ai -= rij*mr3inv
      poti = poti - mrinv
    }
    pfi = (f+i)
    pfi.value.pot =  pfi.value.pot + poti
    pfi.value.acc =   Vec_Float64.new(pfi.value.acc)+ ai
  }
end
