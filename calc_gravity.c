#include <math.h>
typedef struct cvec{
    double x,y,z;
}CVECF64, *PCVECF64;

typedef  struct spjmonopole{
    double mass;
    CVECF64  pos;
}SPJMONOPOLE, *PSPJMONOPOLE;

typedef struct cFull_particle{
    long long id;
    double   mass;
    double eps;
    CVECF64  pos;
    CVECF64  vel;
    double pot;
    CVECF64 acc;
}CFULL_PARTICLE,*PCFULL_PARTICLE;

void  calc_gravity_c_epep(PCFULL_PARTICLE ep_i,
			   int n_ip,
			   PCFULL_PARTICLE ep_j,
			   int n_jp,
			   PCFULL_PARTICLE f)
{
    int i, j;
    for(i=0;i<n_ip;i++){

	PCFULL_PARTICLE pi = ep_i + i;
	double eps2 = pi->eps*pi->eps;
	double xi = pi->pos.x;
	double yi = pi->pos.y;
	double zi = pi->pos.z;
	double ax, ay, az, pot;
	ax = ay = az = pot = 0;
	for(j=0;j<n_jp;j++){
	    PCFULL_PARTICLE pj = ep_j + j;
	    double dx = xi - pj->pos.x;
	    double dy = yi - pj->pos.y;
	    double dz = zi - pj->pos.z;
            double r2 = dx*dx+dy*dy+dz*dz+eps2;
	    double rinv = 1.0/sqrt(r2);
	    double mrinv = pj->mass* rinv;
	    double mr3inv = mrinv*rinv*rinv;
	    ax -= dx*mr3inv;
	    ay -= dy*mr3inv;
	    az -= dz*mr3inv;
	    pot = pot - mrinv;
	}
	PCFULL_PARTICLE  pfi = f+i;
	pfi->pot += pot;
	pfi->acc.x += ax;
	pfi->acc.y += ay;
	pfi->acc.z += az;
    }
}

void  calc_gravity_c_epsp(PCFULL_PARTICLE ep_i,
			   int n_ip,
			   PSPJMONOPOLE ep_j,
			   int n_jp,
			   PCFULL_PARTICLE f)
{
    int i, j;
    for(i=0;i<n_ip;i++){

	PCFULL_PARTICLE pi = ep_i + i;
	double eps2 = pi->eps*pi->eps;
	double xi = pi->pos.x;
	double yi = pi->pos.y;
	double zi = pi->pos.z;
	double ax, ay, az, pot;
	ax = ay = az = pot = 0;
	for(j=0;j<n_jp;j++){
	    PSPJMONOPOLE pj = ep_j + j;
	    double dx = xi - pj->pos.x;
	    double dy = yi - pj->pos.y;
	    double dz = zi - pj->pos.z;
            double r2 = dx*dx+dy*dy+dz*dz+eps2;
	    double rinv = 1.0/sqrt(r2);
	    double mrinv = pj->mass* rinv;
	    double mr3inv = mrinv*rinv*rinv;
	    ax -= dx*mr3inv;
	    ay -= dy*mr3inv;
	    az -= dz*mr3inv;
	    pot = pot - mrinv;
	}
	PCFULL_PARTICLE  pfi = f+i;
	pfi->pot += pot;
	pfi->acc.x += ax;
	pfi->acc.y += ay;
	pfi->acc.z += az;
    }
}
