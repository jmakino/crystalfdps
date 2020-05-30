require "clop"
require "nacsio"
require "./FDPS_vector"
use_fdps_vector(3,Float64)
include FDPS_vector
require "./FDPS_types"
require "./user_defined"

module FLOCAL
Optionstr= <<-END
  Description: Sample program for the use of FDPS from Crystal
  Long description:
    Sample program for the use of FDPS from Crystal
    (c) 2020-, Jun Makino

  Short name:           -n
  Long name:  		--numner-of-particles
  Value type:		int
  Default value: 	0
  Variable name: 	n
  Description:		Number of particles
  Long description:     Number of particles

  Short name: 		-T
  Long name:		--opening_tolerance
  Value type:		float
  Default value: 	0.5
  Variable name: 	tol
  Description:		Opening tolerance
  Long description:
    This option sets the tolerance value that governs the maximum size
    of a tree cell that can remain closed; cells (nodes) with a size
    large than the product of tolerance and distance to that cell will
    be opened, and acceleration to its children will be computed.

  Short name: 		-s
  Long name:		--softening_length
  Value type:		float
  Default value: 	0.05
  Variable name: 	eps
  Description:		Softening length
  Long description:
    This option sets the softening length used to calculate the force
    between two particles.  The calculation scheme comforms to standard
    Plummer softening, where rs2=r**2+eps**2 is used in place of r**2.

  Short name: 		-c
  Long name:		--step_size
  Value type:		float
  Default value:	0.0078125
  Variable name:	dt
  Description:		Time step size
  Long description:
    This option sets the size of the time step, which is constant and
    shared by all particles.  It is wise to use option -s to specify a
    softening length that is significantly larger than the time step size.


  Short name: 		-d
  Long name:		--diagnostics_interval
  Value type:		float
  Default value:	0.25
  Variable name:	dt_dia
  Description:		Interval between diagnostics output
  Long description:
    The time interval between successive diagnostics output.
    The diagnostics include the kinetic and potential energy,
    and the absolute and relative drift of total energy, since
    the beginning of the integration.
        These diagnostics appear on the standard error stream.
    For more diagnostics, try option "-x" or "--extra_diagnostics".

  Short name: 		-o
  Long name:		--output_interval
  Value type:		float
  Default value:	2
  Variable name:	dt_out
  Description:		Time interval between snapshot output
  Long description:
    The time interval between output of a complete snapshot
    A snapshot of an N-body system contains the values of the
    mass, position, and velocity for each of the N particles.

  Short name: 		-O
  Long name:		--output_file_name
  Value type:		string
  Default value:	
  Variable name:	ofname
  Description:	Name for the output snappshot (nacsio) file
  Long description:
    Name for the output snappshot (nacsio) file. Should handle MPI process
    in some way but not yet.

  Short name: 		-t
  Long name:		--duration
  Value type:		float
  Default value:	1
  Variable name:	dt_end
  Description:		Duration of the integration
  Long description:
    This option sets the duration t of the integration, the time period
    after which the integration will halt.  If the initial snapshot is
    marked to be at time t_init, the integration will halt at time
    t_final = t_init + t.

  Short name:		-i
  Long name:  		--init_out
  Value type:  		bool
  Variable name: 	init_out
  Description:		Output the initial snapshot
  Long description:
    If this flag is set to true, the initial snapshot will be output

  Short name:		-x
  Long name:  		--extra_diagnostics
  Value type:  		bool
  Variable name:	x_flag
  Description:		Extra diagnostics
  Long description:
    If this flag is set to true, the following extra diagnostics
    will be printed;
      acceleration (for all integrators)
END
end
clop_init(__LINE__, __FILE__, __DIR__, "FLOCAL::Optionstr")

module Nacsio
  class Particle
    YAML.mapping(
      id: {type: Int64, default: 0i64,},
      time: {type: Float64, key: "t", default: 0.0,},
      mass: {type: Float64, key: "m",default: 0.0,},
      pos: {type: Vector3, key: "r",default: [0.0,0.0,0.0].to_v,},
      vel: {type: Vector3, key: "v",default: [0.0,0.0,0.0].to_v,},
      pot: {type: Float64, key: "pot",default: 0.0,},
    )  
  end
end

include Math
def calc_gravity(ep_i,n_ip,ep_j,n_jp,f)
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
      if r2 > 1e-20
        rinv = 1_f64/sqrt(r2)
        mrinv = pj.mass*rinv
        mr3inv = mrinv*rinv*rinv
        ai -= rij*mr3inv
        poti = poti - mrinv
      end
    }
    pfi = (f+i)
    pfi.value.pot =  pfi.value.pot + poti
    pfi.value.acc =   Vec_Float64.new(pfi.value.acc)+ ai
  }
end

fun init = crystal_init : Void
   GC.init
   LibCrystalMain.__crystal_main(0, Pointer(Pointer(UInt8)).null)
end

module FDPS_vector
  lib FDPS

   enum PS_BOUNDARY_CONDITION 
      BC_OPEN
      BC_PERIODIC_X
      BC_PERIODIC_Y
      BC_PERIODIC_Z
      BC_PERIODIC_XY
      BC_PERIODIC_XZ
      BC_PERIODIC_YZ
      BC_PERIODIC_XYZ
      BC_SHEARING_BOX
      BC_USER_DEFINED
   end

enum PS_INTERACTION_LIST_MODE
  MAKE_LIST
  MAKE_LIST_FOR_REUSE
  REUSE_LIST
end
end
end
require "./FDPS_cr_if.cr"

include Nacsio
fun crmain = crmain(argc : Int32, argv : UInt8**) :  Void
   av = Array(String).new
   argc.times{|i|
     p=argv[i]
     l=0
     while (p+l).value != 0
       l+=1
     end
     av.push String.new(p,l)
   }
   pname = av.shift
   options=CLOP.new(FLOCAL::Optionstr,av)
   crbody(pname, av,options)
   exit 0
end

def return_psys_cptr(psys_num)
  FDPS_calc.get_psys_cptr(psys_num)
end

def  nacs_write(psys_num,time,of)
  # work on non-mpi only yet
  # Local parameters
  ptcl = return_psys_cptr(psys_num)
  p = Particle.new
  FDPS.get_nptcl_loc(psys_num).times{|i|
    q = (ptcl+i)
    p.id = q.value.id
    p.pos =   Vec_Float64.new(q.value.pos).to_a.to_v
    p.vel =   Vec_Float64.new(q.value.vel).to_a.to_v
    p.mass =   q.value.mass
    p.pot =   q.value.pot
    p.time =  time
    of.print p.to_nacs
  }
end

def  setup_IC(psys_num,nptcl_glb, body, options)
  # Local parameters
  nprocs = FDPS.get_num_procs()
  myrank = FDPS.get_rank()
  time=0.0
  if myrank == 0 
    FDPS.set_nptcl_loc(psys_num,nptcl_glb)
    ptcl = return_psys_cptr(psys_num)
    nptcl_glb.times{|i|
      q = (ptcl+i)
      q.value.id = body[i].id
      q.value.pos = Vec_Float64.new(body[i].pos.to_a)
      q.value.vel = Vec_Float64.new(body[i].vel.to_a)
      q.value.eps = options.eps
      q.value.mass = body[i].mass
    }
    time = body[0].time
  else
    FDPS.set_nptcl_loc(psys_num,0)
  end
  FDPS.broadcast_f64(pointerof(time),1,0)
  time
end

def  setup_IC(psys_num,nptcl_glb)
  # Local parameters
    m_tot=1.0_f64
    rmax=3.0_f64
    r2max=rmax*rmax
    #   Get # of MPI processes and rank number
   nprocs = FDPS.get_num_procs()
   myrank = FDPS.get_rank()
   # Make an initial condition at RANK 0
   if myrank == 0 
      # Set # of local particles
      FDPS.set_nptcl_loc(psys_num,nptcl_glb)
      #!* Create an uniform sphere of particles
      #!** get the pointer to full particle data
      ptcl = return_psys_cptr(psys_num)
      #!** initialize Mersenne twister
      FDPS.mt_init_genrand(0)
      zerov= Vec_Float64.new(0)
      nptcl_glb.times{|i|
        q = (ptcl+i)
        q.value.id = i
        q.value.mass = m_tot/nptcl_glb
        r2 = r2max*2
        while r2 >= r2max
          rv = Vec_Float64.new(Array.new(3){
                                 (2_f64*FDPS.mt_genrand_res53()-1.0) * rmax})
          r2 = rv*rv
          q.value.pos = rv
        end
        q.value.vel = zerov
        q.value.eps = 0.1
      }
      
      cm_pos = zerov
      cm_vel = zerov
      cm_mass = 0_f64
      nptcl_glb.times{|i|
        pi = (ptcl+i).value
        cm_pos +=  Vec_Float64.new(pi.pos)* pi.mass
        cm_vel +=  Vec_Float64.new(pi.vel)* pi.mass
        cm_mass += pi.mass
      }
      cm_pos = cm_pos/cm_mass
      cm_vel = cm_vel/cm_mass
      nptcl_glb.times{|i|
        q = ptcl+i
        q.value.pos = Vec_Float64.new(q.value.pos)-cm_pos
        q.value.vel = Vec_Float64.new(q.value.vel)-cm_vel
      }
   else
     FDPS.set_nptcl_loc(psys_num,0)
   end
end
def  setup_binary(psys_num)
  # Local parameters
    m_tot=1.0_f64
   nprocs = FDPS.get_num_procs()
   myrank = FDPS.get_rank()
   if myrank == 0 
      FDPS.set_nptcl_loc(psys_num,2)
      ptcl = return_psys_cptr(psys_num)
      q0 = (ptcl)
      q1 = (ptcl+1)
      q0.value.id = 0
      q1.value.id = 1
      q0.value.mass = 0.5
      q1.value.mass = 0.5
      q0.value.pos = Vec_Float64.new([0.5,0.0,0.0])
      q1.value.pos = Vec_Float64.new([-0.5,0.0,0.0])
      q0.value.vel = Vec_Float64.new([0.0,0.5,0.0])
      q1.value.vel = Vec_Float64.new([0.0,-0.5,0.0])
      q0.value.eps = 0.0
      q1.value.eps = 0.0
   else
     FDPS.set_nptcl_loc(psys_num,0)
   end
end

def  dump_particles(psys_num)
  ptcl = return_psys_cptr(psys_num)
  FDPS.get_nptcl_loc(psys_num).times{|i|
    pp! i
    pi = (ptcl+i).value
    p pi
  }
end
def calc_force_all_and_write_back(tree_num,
                                      calc_gravity_pp,
                                      calc_gravity_psp,
                                      psys_num,  
                                      dinfo_num)

  FDPS_calc.calc_force_all_and_write_back_local(tree_num,
                                             calc_gravity_pp,
                                             calc_gravity_psp,
                                             psys_num,  
                                             dinfo_num,
                                             true,
                                             FDPS::PS_INTERACTION_LIST_MODE::MAKE_LIST)
end

def calc_energy(psys_num)
  etot = 0_f64
  ekin = 0_f64
  epot = 0_f64
   # Get # of local particles
  nptcl_loc = FDPS.get_nptcl_loc(psys_num)
  dummy = FDPS::Full_particle.new
  ptcl = return_psys_cptr(psys_num)
  
  ekin_loc = 0_f64
  epot_loc = 0_f64
  nptcl_loc.times{|i|
    pi = (ptcl+i).value
    v = Vec_Float64.new(pi.vel)
    ekin_loc += pi.mass * (v*v)
    epot_loc += pi.mass * pi.pot 
    epot_loc -= pi.mass * pi.mass/pi.eps if pi.eps != 0.0
  }
   ekin_loc *= 0.5_f64
   epot_loc *= 0.5_f64
   etot_loc = ekin_loc + epot_loc
   ekin= FDPS.get_sum_f64(ekin_loc)
   epot= FDPS.get_sum_f64(epot_loc)
   etot=FDPS.get_sum_f64(etot_loc)
   [etot,ekin,epot]
end

def kick(psys_num,dt)
  ptcl = return_psys_cptr(psys_num)
  FDPS.get_nptcl_loc(psys_num).times{|i|
    q = ptcl+i
    q.value.vel= add(q.value.vel, mul(q.value.acc, dt))

  }
end
def drift(psys_num,dt)
  ptcl = return_psys_cptr(psys_num)
  FDPS.get_nptcl_loc(psys_num).times{|i|
    q = ptcl+i
    q.value.pos= add(q.value.pos, mul(q.value.vel, dt))
  }
end

lib FDPS_calc
fun calc_force_all_and_write_back_local=fdps_calc_force_all_and_write_back(tree_num : Int32,
                    pfunc_ep_ep:
                         Pointer(FDPS::Full_particle), Int32,
                         Pointer(FDPS::Full_particle), Int32,
                         Pointer(FDPS::Full_particle) -> Void,
                    pfunc_ep_sp:
                         Pointer(FDPS::Full_particle), Int32,
                         Pointer(FDPS::SPJMonopole), Int32,
                         Pointer(FDPS::Full_particle) -> Void,
           psys_num : Int32,
           dinfo_num : Int32,
           clear : Bool,
           list_mode : FDPS::PS_INTERACTION_LIST_MODE) : Void
fun get_psys_cptr = fdps_get_psys_cptr(psys_num : Int32): Pointer(FDPS::Full_particle)

end

def crbody(pname, av,options)
  STDERR.print "FDPS on Crystal test code\n"
  FDPS.initialize()
  pp! options
  dinfo_num : Int32 =1
  coef_ema : Float32 =0.3
  FDPS.create_dinfo(pointerof(dinfo_num))
  FDPS.init_dinfo(dinfo_num,coef_ema)
  psys_num : Int32 = 0
  FDPS.create_psys(pointerof(psys_num),"full_particle")
  FDPS.init_psys(psys_num)
  #  Create tree object
  tree_num = 0
  FDPS.create_tree(pointerof(tree_num), 
                   "Long,full_particle,full_particle,full_particle,Monopole")

  of=STDOUT
  if options.ofname.size > 0
    nprocs = FDPS.get_num_procs()
    if nprocs > 1
      maxdigits = nprocs.to_s.size
      myrank = FDPS.get_rank()
      ofname = options.ofname+"0"* (maxdigits- myrank.to_s.size)+myrank.to_s
    else
      ofname = options.ofname
    end
    of=File.open(ofname, "w")
  end
  time_sys  = 0_f64
  if options.n > 0
    if FDPS.get_rank()==0
      of.print CommandLog.new("treecode with internal IC",pname,av).to_nacs
    end
    ntot=options.n
    if ntot == 2
      setup_binary(psys_num)
    else
      setup_IC(psys_num,ntot)
    end
  else
    body = Array(Particle).new
    ybody = Array(YAML::Any).new
    ntot=0
    if FDPS.get_rank()==0
      update_commandlog(pname, av,of)
      while (sp= CP(Particle).read_particle).y != nil
        body.push sp.p
        ybody.push sp.y
      end
      ntot = body.size
    end
    time_sys = setup_IC(psys_num,ntot,body,options)
  end
    
  n_leaf_limit = 8
  n_group_limit = 16
  FDPS.init_tree(tree_num,ntot, options.tol, n_leaf_limit, n_group_limit)

  # !* Domain decomposition and exchange particle
  FDPS.decompose_domain_all(dinfo_num,psys_num,1)
  FDPS.exchange_particle(psys_num,dinfo_num)
  # !* Compute force at the initial time
  ppf = ->(ep_i : Pointer(FDPS::Full_particle),
           n_ip : Int32,
           ep_j : Pointer(FDPS::Full_particle),
           n_jp : Int32,
           f : Pointer(FDPS::Full_particle)){
    #     C.calc_gravity_c_epep(ep_i,n_ip,ep_j,n_jp,f)
    calc_gravity(ep_i,n_ip,ep_j,n_jp,f)
  }
   psf = ->(ep_i : Pointer(FDPS::Full_particle),
                    n_ip : Int32,
                    ep_j : Pointer(FDPS::SPJMonopole),
                    n_jp : Int32,
            f : Pointer(FDPS::Full_particle)){
#     C.calc_gravity_c_epsp(ep_i,n_ip,ep_j,n_jp,f)
     calc_gravity(ep_i,n_ip,ep_j,n_jp,f)
   }
   tprof = FDPS::TimeProfile.new
   FDPS.clear_tree_time_prof(tree_num)
   calc_force_all_and_write_back(tree_num, ppf,psf,
                                      psys_num,  
                                      dinfo_num)
   # !* Compute energies at the initial time
   etot0,ekin0,epot0 = calc_energy(psys_num)
   # !* Time integration
   time_diag = time_sys
   time_snap = time_sys
   time_end = time_sys + options.dt_end
   dt = options.dt
   dt_diag = options.dt_dia
   dt_snap = options.dt_out
   num_loop = 0
   while time_sys <= time_end
     if time_sys + dt/2 >= time_snap
       nacs_write(psys_num,time_sys,of)
       time_snap += dt_snap
     end
     etot1,ekin1,epot1 = calc_energy(psys_num)
     if FDPS.get_rank() == 0
       if time_sys + dt/2 >= time_diag
         STDERR.print "time: #{time_sys}, energy error: #{(etot1-etot0)/etot0}\n"
         time_diag = time_diag + dt_diag
#         FDPS.get_tree_time_prof(tree_num,pointerof(tprof))
#         p tprof
       end
     end
     kick(psys_num,0.5_f64*dt)
     time_sys +=  dt
     drift(psys_num,dt)
     #    !* Domain decomposition & exchange particle
     if num_loop%4 == 0
       FDPS.decompose_domain_all(dinfo_num,psys_num,1)
     end
     FDPS.exchange_particle(psys_num,dinfo_num)
     calc_force_all_and_write_back(tree_num, ppf,psf,
                                   psys_num,  
                                   dinfo_num)
     kick(psys_num,0.5_f64*dt)
     num_loop += 1
   end
   of.close if options.ofname.size > 0
   FDPS.finalize()
end
