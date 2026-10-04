/* formation_run.c  [CL draft, 4 October 2026]
 * Navigator's Log R&D, Tunable Rocky Planet Water Model, Stage 2.
 * Late-stage rocky-planet accretion with fragmentation, for benchmark B0 and comparison P-D12.
 * Toolchain (pinned): REBOUND 4.4.10 (C) + Childs & Steffen fragmentation module, commit d67ab2b.
 * Build: scripts/setup.sh  (copies fragmentation.c next to this file as frag.c)
 * Usage: ./formation_run <seed> <expansion_factor> <tmax_yr> <walltime_s>
 *   Resumes automatically from run.bin + state.txt if they exist.
 * Outputs (current directory): run.bin (simulation archive), state.txt (fragment counter),
 *   progress.csv (t, N, loose-material mass, collision-resolve calls, walltime), collision_report.txt (module).
 * Disk: Chambers (2013) bimodal setup as used in Tajer et al. (2025, arXiv:2511.01842):
 *   14 embryos and 140 planetesimals, 0.3 to 2.0 AU, Jupiter and Saturn on current orbits.
 *   Embryo and planetesimal masses 0.093 and 0.0093 Earth masses: [CL inference] Tajer et al. round
 *   these to 0.1 and 0.01, but their stated total (2.604) matches 14 x 0.093 + 140 x 0.0093 exactly.
 *   Surface density: rises linearly from zero at 0.3 AU to 0.7 AU, then a^-3/2 to 2.0 AU
 *   [CL assumption: the starting value at 0.3 AU is not stated in the source read; zero is assumed].
 *   Semi-major axes are drawn from that profile for each population separately.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "frag.c"   /* defines min_frag_mass, tot_no_frags, reb_collision_resolve_fragment */

#define MEARTH 3.003e-6           /* Earth mass in solar masses */
#define LOOSE_CUT (0.1*MEARTH)    /* loose material: bodies under 0.1 Earth mass ([102] yardstick) */

static double f_exp = 1.0;
static double wall_limit = 0;
static long ncoll = 0;
static double mloose0 = 0;
static FILE* prog = NULL;
static int stop_flag = 0;

int resolve(struct reb_simulation* const r, struct reb_collision c){ ncoll++; return reb_collision_resolve_fragment(r,c); }

static int is_big(struct reb_particle* p){ return p->hash==0 || p->hash==reb_hash("JUPITER") || p->hash==reb_hash("SATURN"); }

static double loose_mass(struct reb_simulation* r){
    double m=0; for(int i=0;i<r->N;i++){ struct reb_particle* p=&r->particles[i]; if(!is_big(p) && p->m<LOOSE_CUT) m+=p->m; } return m; }

static void save_state(struct reb_simulation* r){
    reb_simulation_save_to_file(r,"run.bin");
    FILE* s=fopen("state.txt","w"); fprintf(s,"%d %ld %.17g\n",tot_no_frags,ncoll,mloose0); fclose(s); }

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r,1.e4)){
        fprintf(prog,"%.6e,%d,%.10e,%ld,%.1f\n",r->t,r->N,loose_mass(r)/MEARTH,ncoll,r->walltime); fflush(prog); }
    if(reb_simulation_output_check(r,1.e5)) save_state(r);
    if(wall_limit>0 && r->walltime>wall_limit){ stop_flag=1; reb_simulation_stop(r); }
}

/* inverse-CDF draw of semi-major axis from the assumed surface density profile */
static double draw_a(struct reb_simulation* r){
    static double grid[4001], cdf[4001]; static int init=0;
    if(!init){ double s0=pow(0.7,-1.5); cdf[0]=0;
        for(int k=0;k<=4000;k++){ double a=0.3+1.7*k/4000.; grid[k]=a;
            double sig = (a<0.7)? s0*(a-0.3)/0.4 : pow(a,-1.5);
            if(k>0){ double ap=grid[k-1]; double sp=(ap<0.7)? s0*(ap-0.3)/0.4 : pow(ap,-1.5);
                cdf[k]=cdf[k-1]+0.5*(sig*a+sp*ap)*(a-ap); } }
        for(int k=0;k<=4000;k++) cdf[k]/=cdf[4000]; init=1; }
    double u=reb_random_uniform(r,0,1); int k=1; while(k<4000 && cdf[k]<u) k++;
    double w=(u-cdf[k-1])/(cdf[k]-cdf[k-1]+1e-300); return grid[k-1]+w*(grid[k]-grid[k-1]);
}

int main(int argc, char* argv[]){
    if(argc<5){ fprintf(stderr,"usage: %s seed factor tmax_yr walltime_s\n",argv[0]); return 1; }
    int seed=atoi(argv[1]); f_exp=atof(argv[2]); double tmax=atof(argv[3]); wall_limit=atof(argv[4]);
    min_frag_mass = 0.5*0.0093*MEARTH;      /* half a planetesimal mass, as in [104] */
    const double rho = 5.05e6;              /* 3 g/cm^3 in Msun/AU^3, as in the module's example */
    struct reb_simulation* r;
    int resumed = (access("run.bin",F_OK)==0 && access("state.txt",F_OK)==0);
    if(resumed){
        r = reb_simulation_create_from_file("run.bin",-1);
        FILE* s=fopen("state.txt","r"); if(fscanf(s,"%d %ld %lg",&tot_no_frags,&ncoll,&mloose0)!=3){fprintf(stderr,"bad state\n");return 2;} fclose(s);
        r->walltime = 0;   /* per-session wall clock */
    } else {
        r = reb_simulation_create();
        r->G = 39.476926421373; r->rand_seed = seed;
        reb_simulation_add_fmt(r,"m r hash",1.0,0.00465,0u);   /* solar radius in AU */
        struct reb_particle star=r->particles[0];
        for(int i=0;i<154;i++){
            double m = (i<14)? 0.093*MEARTH : 0.0093*MEARTH;
            double a = draw_a(r);
            double e = reb_random_uniform(r,0,0.01), inc = reb_random_uniform(r,0,0.0175);
            double om = reb_random_uniform(r,0,2*M_PI), Om = reb_random_uniform(r,0,2*M_PI), fa = reb_random_uniform(r,0,2*M_PI);
            struct reb_particle p = reb_particle_from_orbit(r->G,star,m,a,e,inc,Om,om,fa);
            p.r = get_radii(m,rho)*f_exp; p.hash = i+1;
            reb_simulation_add(r,p);
        }
        struct reb_particle J = reb_particle_from_orbit(r->G,star,9.543e-4,5.20349,0.048381,0.365*M_PI/180,0.0,68.3155*M_PI/180,227.0537*M_PI/180);
        J.r = get_radii(J.m,rho); J.hash = reb_hash("JUPITER"); reb_simulation_add(r,J);
        struct reb_particle S = reb_particle_from_orbit(r->G,star,2.857e-4,9.54309,0.052519,0.8892*M_PI/180,M_PI,324.5263*M_PI/180,256.9188*M_PI/180);
        S.r = get_radii(S.m,rho); S.hash = reb_hash("SATURN"); reb_simulation_add(r,S);
        reb_simulation_move_to_com(r);
        mloose0 = loose_mass(r);
    }
    r->dt = 6./365.25;
    r->integrator = REB_INTEGRATOR_TRACE;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = resolve;
    r->heartbeat = heartbeat;
    prog = fopen("progress.csv", resumed? "a":"w");
    if(!resumed){ fprintf(prog,"# seed=%d factor=%.3f mloose0_Mearth=%.6f\nt_yr,N,loose_Mearth,resolve_calls,walltime_s\n",seed,f_exp,mloose0/MEARTH); }
    reb_simulation_integrate(r,tmax);
    save_state(r);
    fprintf(prog,"# session end t=%.6e N=%d stopped_by_walltime=%d\n",r->t,r->N,stop_flag); fclose(prog);
    reb_simulation_free(r);
    return 0;
}
