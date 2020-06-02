#pragma once
/* Standard headers */
#include <math.h>
/* FDPS headers */
#include "FDPS_c_if.h"
typedef struct full_particle{  //$fdps FP,EPI,EPJ,Force
      //$fdps copyFromForce full_particle (pot,pot) (acc,acc)
      //$fdps copyFromFP full_particle (id,id) (mass,mass) (eps,eps) (pos,pos)
      //$fdps clear id=keep, mass=keep, eps=keep, pos=keep, vel=keep
      long long id;
      double mass; //$fdps charge
      double eps;
      fdps_f64vec pos; //$fdps position
      fdps_f64vec vel; //$fdps velocity
      double pot;
      fdps_f64vec acc;
} Full_particle;
