require "./FDPS_vector"
use_fdps_vector(3,Float64)
include FDPS_vector
require "./FDPS_super_particle"
require "./user_defined"


fun init = crystal_init : Void
  GC.init
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


def calc_energy(psys_num,etot,ekin,epot,clear)
  if clear 
      etot = 0_f64
      ekin = 0_f64
      epot = 0_f64
   end

   # Get # of local particles
  nptcl_loc = FDPS.get_nptcl_loc(psys_num)
  dummy = FDPS::Full_particle.new
  ptcl = pointerof(dummy)
  FDPS.get_psys_cptr(psys_num,ptcl)
  
  ekin_loc = 0_f64
  epot_loc = 0_f64
  nptcl_loc.times{|i|
    pi = (ptcl+i).value
    ekin_loc += pi.mass * (pi.vel * pi.vel)
    epot_loc += pi.mass * (pi.pot + pi.mass/pi.eps)
  }
   ekin_loc *= 0.5_f64
   epot_loc *= 0.5_f64
   etot_loc = ekin_loc + epot_loc
   call fdps_ctrl%get_sum(ekin_loc,ekin)
   call fdps_ctrl%get_sum(epot_loc,epot)
   call fdps_ctrl%get_sum(etot_loc,etot)

   !* Release the pointer
   nullify(ptcl)

end subroutine calc_energy

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
  ntot=16
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
   p "after force calculation\n"
   dump_particles(psys_num)
   # !* Compute energies at the initial time
   clear = true
   calc_energy(fdps_ctrl,psys_num,etot0,ekin0,epot0,clear)

   # !* Time integration
   # time_diag = 0.0d0
   # time_snap = 0.0d0
   # time_sys  = 0.0d0
   # num_loop = 0
   # do 
   #    !* Output
   #   !if (FDPS.get_rank() == 0) then
   #   !   write(*,50)num_loop,time_sys
   #   !   50 format('(num_loop, time_sys) = ',i5,1x,1es25.16e3)
   #   !end if
   #    if ( (time_sys >= time_snap) .or. &
   #         (((time_sys + dt) - time_snap) > (time_snap - time_sys)) ) then
   #       call output(fdps_ctrl,psys_num)
   #       time_snap = time_snap + dt_snap
   #    end if
      
   #    !* Compute energies and output the results
   #    clear = .true.
   #    call calc_energy(fdps_ctrl,psys_num,etot1,ekin1,epot1,clear)
   #    if (FDPS.get_rank() == 0) then
   #       if ( (time_sys >= time_diag) .or. &
   #            (((time_sys + dt) - time_diag) > (time_diag - time_sys)) ) then
   #          write(*,100)time_sys,(etot1-etot0)/etot0
   #          100 format("time: ",1es20.10e3,", energy error: ",1es20.10e3)
   #          time_diag = time_diag + dt_diag
   #       end if
   #    end if

   #    !* Leapfrog: Kick-Drift
   #    call kick(fdps_ctrl,psys_num,0.5d0*dt)
   #    time_sys = time_sys + dt
   #    call drift(fdps_ctrl,psys_num,dt)

   #    !* Domain decomposition & exchange particle
   #    if (mod(num_loop,4) == 0) then
   #       FDPS.decompose_domain_all(dinfo_num,psys_num)
   #    end if
   #    FDPS.exchange_particle(psys_num,dinfo_num)

   #    !* Force calculation
   #    pfunc_ep_ep = c_funloc(calc_gravity_pp)
   #    pfunc_ep_sp = c_funloc(calc_gravity_psp)
   #    FDPS.calc_force_all_and_write_back(tree_num,    &
   #                                                 pfunc_ep_ep, &
   #                                                 pfunc_ep_sp, &
   #                                                 psys_num,    &
   #                                                 dinfo_num)
   #    !* Leapfrog: Kick
   #    call kick(fdps_ctrl,psys_num,0.5d0*dt)

   #    !* Update num_loop
   #    num_loop = num_loop + 1

   #    !* Termination
   #    if (time_sys >= time_end) then
   #       exit
   #    end if
   # end do


  FDPS.PS_Finalize()
end

