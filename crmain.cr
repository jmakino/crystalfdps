require "./FDPS_vector"
use_fdps_vector(3,Float64)
include FDPS_vector
require "./FDPS_super_particle"
require "./user_defined"


fun init = crystal_init : Void
#  GC.init
  STDERR.print "Crystal initialization end\n"
end

lib C
  fun cos(value : Float64) : Float64
end
module FDPS_vector
lib FDPS
fun PS_Initialize  = fdps_initialize : Void
fun PS_Finalize  = fdps_finalize : Void
fun create_dinfo = fdps_create_dinfo(id : Pointer(Int32)) : Void
fun init_dinfo = fdps_init_dinfo(dinfo_num : Int32,
                                 coef_ema : Float32) : Void
fun create_psys = fdps_create_psys(psys_num : Int32*,
                      psys_info: UInt8*) : Void
fun init_psys = fdps_init_psys(psys_num : Int32) : Void
fun create_tree = fdps_create_tree(tree_num : Int32*,
                      tree_info : UInt8*) : Void
fun init_tree =  fdps_init_tree(tree_num : Int32,
                    nptcl : Int32,
                    theta : Float32,
                    n_leaf_limit : Int32,
                                n_group_limit : Int32) : Void
fun get_num_procs = fdps_get_num_procs : Int32
fun get_rank = fdps_get_rank : Int32
fun set_nptcl_loc = fdps_set_nptcl_loc(psys_num : Int32, nptcl: Int32) : Void
fun get_psys_cptr = fdps_get_psys_cptr(psys_num : Int32,
                                       cptr:  (Full_particle)**) : Void
fun MT_init_genrand = fdps_mt_init_genrand(s : Int32) : Void
fun MT_genrand_res53= fdps_mt_genrand_res53 : Float64
fun decompose_domain_all=fdps_decompose_domain_all(dinfo_num : Int32,
                                                   psys_num : Int32,
                                       weight : Float32) : Void
fun exchange_particle = fdps_exchange_particle(psys_num : Int32, 
                                               dinfo_num : Int32) : Void

fun fdps_calc_force_all_and_write_back_l00000(tree_num : Int32,
                                            pfunc_ep_ep :
                                              Pointer(FDPS::Full_particle),
                                              Int32,
                                              Pointer(FDPS::Full_particle),
                                              Int32,
                                           Pointer(FDPS::Full_particle) -> Void,
                                            pfunc_ep_sp :
                                              Pointer(FDPS::Full_particle),
                                              Int32,
                                              Pointer(FDPS::Spj_monopole),
                                              Int32,
                                           Pointer(FDPS::Full_particle) -> Void,
                                           psys_num : Int32,
                                            dinfo_num : Int32,
                                            clear : Bool,
                                            list_mode : Int32)
fun  get_nptcl_loc=fdps_get_nptcl_loc(psys_num : Int32)  : Int32
fun  get_sum_f64 = fdps_get_sum_r64(lv : Float64, gv : Float64*) : Void 
end
end


fun crmain :  Void
  crbody
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
      dummy = FDPS::Full_particle.new
      ptcl = pointerof(dummy)
      FDPS.get_psys_cptr(psys_num,pointerof(ptcl))
      #!** initialize Mersenne twister
      FDPS.MT_init_genrand(0)
      zerov= Vec_Float64.new(0)
      nptcl_glb.times{|i|
        q = (ptcl+i)
        q.value.id = i
        q.value.mass = m_tot/nptcl_glb
        r2 = r2max*2
        while r2 >= r2max
          rv = Vec_Float64.new(Array.new(3){
                                 (2_f64*FDPS.MT_genrand_res53()-1.0) * rmax})
          r2 = rv*rv
          q.value.pos = rv
        end
        q.value.vel = zerov
        q.value.eps = 1_f64/32_f64
      }
      
      # nptcl_glb.times{|i|
      #   pi = (ptcl+i).value
      #   p pi
      # }
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

def  dump_particles(psys_num)
  dummy = FDPS::Full_particle.new
  ptcl = pointerof(dummy)
  FDPS.get_psys_cptr(psys_num,pointerof(ptcl))
  FDPS.get_nptcl_loc(psys_num).times{|i|
    pi = (ptcl+i).value
    p pi
  }
end

def calc_force_all_and_write_back(tree_num,
                                      calc_gravity_pp,
                                      calc_gravity_psp,
                                      psys_num,  
                                      dinfo_num)

#  p "call   FDPS.fdps_calc_force_and_write_back_l00000\n"
  FDPS.fdps_calc_force_all_and_write_back_l00000(tree_num,
                                             calc_gravity_pp,
                                             calc_gravity_psp,
                                             psys_num,  
                                             dinfo_num,
                                             true,
                                             0)

#  p "return from  FDPS.fdps_calc_force_and_write_back_l00000\n"
end


def calc_energy(psys_num)
  etot = 0_f64
  ekin = 0_f64
  epot = 0_f64
   # Get # of local particles
  nptcl_loc = FDPS.get_nptcl_loc(psys_num)
  dummy = FDPS::Full_particle.new
  ptcl = pointerof(dummy)
  FDPS.get_psys_cptr(psys_num,pointerof(ptcl))
  
  ekin_loc = 0_f64
  epot_loc = 0_f64
  nptcl_loc.times{|i|
    pi = (ptcl+i).value
#    p pi
    v = Vec_Float64.new(pi.vel)
    ekin_loc += pi.mass * (v*v)
    epot_loc += pi.mass * (pi.pot + pi.mass/pi.eps)
#    p pi.pot
  }
   ekin_loc *= 0.5_f64
   epot_loc *= 0.5_f64
#   p epot_loc
#   p ekin_loc
   etot_loc = ekin_loc + epot_loc
   FDPS.get_sum_f64(ekin_loc,pointerof(ekin))
   FDPS.get_sum_f64(epot_loc,pointerof(epot))
   FDPS.get_sum_f64(etot_loc,pointerof(etot))
   [etot,ekin,epot]
end

def kick(psys_num,dt)
  dummy = FDPS::Full_particle.new
  ptcl = pointerof(dummy)
  FDPS.get_psys_cptr(psys_num,pointerof(ptcl))
  FDPS.get_nptcl_loc(psys_num).times{|i|
    q = ptcl+i
    q.value.vel = Vec_Float64.new(q.value.vel) +
                  Vec_Float64.new(q.value.acc)* dt
  }
end
def drift(psys_num,dt)
  dummy = FDPS::Full_particle.new
  ptcl = pointerof(dummy)
  FDPS.get_psys_cptr(psys_num,pointerof(ptcl))
  FDPS.get_nptcl_loc(psys_num).times{|i|
    q = ptcl+i
    q.value.pos = Vec_Float64.new(q.value.pos) +
                  Vec_Float64.new(q.value.vel)* dt
  }
end

def crbody
  STDERR.print "FDPS on Crystal test code\n"
  FDPS.PS_Initialize()
  dinfo_num : Int32 =1
  coef_ema : Float32 =0.3

  FDPS.create_dinfo(pointerof(dinfo_num))
  p dinfo_num
  FDPS.init_dinfo(dinfo_num,coef_ema)

  psys_num : Int32 = 0
  FDPS.create_psys(pointerof(psys_num),"full_particle")
  FDPS.init_psys(psys_num)

  #  Create tree object
  tree_num = 0
  FDPS.create_tree(pointerof(tree_num), 
                   "Long,full_particle,full_particle,full_particle,Monopole")
  ntot=1024
  theta = 0.5
   n_leaf_limit = 8
   n_group_limit = 64
  FDPS.init_tree(tree_num,ntot, theta, n_leaf_limit, n_group_limit)

   # !* Make an initial condition
   setup_IC(psys_num,ntot)
   # !* Domain decomposition and exchange particle
   FDPS.decompose_domain_all(dinfo_num,psys_num,1)
   FDPS.exchange_particle(psys_num,dinfo_num)

   # !* Compute force at the initial time
   ppf = ->(ep_i : Pointer(FDPS::Full_particle),
                    n_ip : Int32,
                    ep_j : Pointer(FDPS::Full_particle),
                    n_jp : Int32,
            f : Pointer(FDPS::Full_particle)){
     calc_gravity(ep_i,n_ip,ep_j,n_jp,f)
   }
   psf = ->(ep_i : Pointer(FDPS::Full_particle),
                    n_ip : Int32,
                    ep_j : Pointer(FDPS::Spj_monopole),
                    n_jp : Int32,
            f : Pointer(FDPS::Full_particle)){
     calc_gravity(ep_i,n_ip,ep_j,n_jp,f)
   }
   calc_force_all_and_write_back(tree_num, ppf,psf,
                                      psys_num,  
                                      dinfo_num)
#   p "after force calculation\n"
#   dump_particles(psys_num)
   # !* Compute energies at the initial time
   etot0,ekin0,epot0 = calc_energy(psys_num)
   print "Energies = ", [etot0,ekin0,epot0].map{|x| x.to_s}.join(" "), "\n"
   # !* Time integration
   time_diag = 0_f64
   time_snap = 0_f64
   time_sys  = 0_f64
   time_end = 10.0_f64
   dt = 1.0_f64/128.0_f64
   dt_diag = 1.0_f64
   dt_snap = 1.0_f64
   num_loop = 0
   while time_sys <= time_end
     if time_sys + dt/2 >= time_snap
       #       output(psys_num)
       time_snap += dt_snap
     end
     
     etot1,ekin1,epot1 = calc_energy(psys_num)
     if FDPS.get_rank() == 0
       if time_sys + dt/2 >= time_diag
         print "time: #{time_sys}, energy error: #{(etot1-etot0)/etot0}\n"
         time_diag = time_diag + dt_diag
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
     #    !* Force calculation
     
     calc_force_all_and_write_back(tree_num, ppf,psf,
                                   psys_num,  
                                   dinfo_num)
     kick(psys_num,0.5_f64*dt)
     num_loop += 1
   end
   FDPS.PS_Finalize()
end

