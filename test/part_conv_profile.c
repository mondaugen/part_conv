/* Firstly, is the algorithm correct?
   For which value of D is the algorithm fastest? */

#include "part_conv.h"

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h> 
#include <time.h>

#define ERR_TOL 1.e-6

#ifndef  FABS
#define FABS(x) ((x)<0.)?((-1)*(x)):(x) 
#endif  

#define TESTCHK(x,s) \
    if (x) { puts(s " passed."); } else { puts(s " failed."); }

#define PC_DEF(M,N) \
    static part_conv_t *pc_new_ ## M ## _ ## N (float *ir) { \
      part_conv_t *pc = part_conv_new(M,N); part_conv_set_ir_td(pc,ir); return pc;\
    }; \
    static void pc_call_ ## M ## _ ## N (part_conv_t *pc, float *x, size_t L) { \
      while ((L -= M) >= M) { part_conv_proc(pc,x); x += M; } \
    }

#define PC_NEW_CALL(M,N,ir)  pc_new_ ## M ## _ ## N (ir)
#define PC_CALL_CALL(M,N,pc,x,L)  pc_call_ ## M ## _ ## N (pc, x, L)

int check_close(float *x, float *y, size_t n)
{
    int ret = 1;
    while (n--) { ret = ret && (FABS(x[n]-y[n]) < ERR_TOL); }
    return ret;
}

int part_conv_correct_test (void)
{
    float in[] = {1.,0,.5,0,0,0,0,0},
          ir[] = {0.5, 0.25, -0.5, -0.25},
          outc[] = {0.5, 0.25, -0.5+0.25, -0.25+0.125,-0.25,-0.125,0.,0.};
    part_conv_t *pc = part_conv_new(2,4); if (!pc) { return -1; }
    part_conv_set_ir_td(pc,ir);
    part_conv_proc(pc,in);
    part_conv_proc(pc,in+2);
    part_conv_proc(pc,in+4);
    part_conv_proc(pc,in+6);
    return check_close(in,outc,8);
}

int part_conv_simple_test(void)
{
    size_t M = 64;
    size_t N_x = M*1000;
    size_t N_ir = 8192;
    size_t N_xp = 11;
    float *x = (float*)calloc(N_x+N_ir-1,sizeof(float));
    float *ir = (float*)calloc(N_ir,sizeof(float));
    float **xs = (float**)calloc(N_xp+1,sizeof(float*));
    size_t WTF = N_x+N_ir-1;
    if (!x) { return 0; }
    srandom(time(NULL));
    {size_t n; for (n=0;n<N_x;n++) {x[n] = 2.*random()/(float)RAND_MAX - 1.; }}
    {size_t n; for (n=0;n<N_ir;n++) {ir[n] = 2.*random()/(float)RAND_MAX - 1.; }}
    {size_t n; for (n=0;n<N_xp;n++) { xs[n] = (void*)malloc((WTF)*sizeof(float));
                                       memcpy(xs[n],x,sizeof(float)*(WTF)); }}
    part_conv_t *pc = part_conv_new(M,N_ir);
    if (!pc) { return 0; }
    part_conv_set_ir_td(pc,ir);
    size_t L = WTF;
    float *x_ = x;
    while ((L -= M) >= M) { part_conv_proc(pc,x_); x_+=M; }
    return 1;
}

int part_conv_correct_test_2(void)
{
    size_t M = 8;
    size_t N_x = 64;
    size_t N_ir = 16;
    float x[] = { 0.08377305,  0.84906741,  0.08804433,  0.18773621,  0.0555924 ,
         0.68319716,  0.54277709,  0.8686232 ,  0.11818004,  0.35001514,
         0.17432048,  0.62206198,  0.9518995 ,  0.39974543,  0.24092933,
         0.97722383,  0.30898167,  0.74782438,  0.52189031,  0.56526228,
         0.05248395,  0.16949527,  0.46590737,  0.62066344,  0.9610273 ,
         0.0877524 ,  0.57863575,  0.96633556,  0.27341699,  0.18526249,
         0.42843129,  0.22942503,  0.94893393,  0.39566267,  0.70052593,
         0.53125384,  0.30041993,  0.82365971,  0.9858987 ,  0.75307173,
         0.74418494,  0.40327832,  0.98333168,  0.55933785,  0.62888525,
         0.30888951,  0.83069315,  0.1712736 ,  0.09555093,  0.96634021,
         0.41527469,  0.1380426 ,  0.72644052,  0.60837661,  0.10090554,
         0.12326398,  0.96052138,  0.08219205,  0.22688496,  0.38074259,
         0.01297886,  0.32387302,  0.42417337,  0.73018909,
         0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0.,
         0.,  0.,  0.,  0.,  0., 0.};
    float ir[] = { 0.99851366,  0.48396177,  0.26439248,  0.29038465,  0.39348695,
         0.3234635 ,  0.76323747,  0.30602051,  0.35906104,  0.42804507,
         0.62423667,  0.09478329,  0.46949767,  0.76674858,  0.0530853 ,
         0.55425212};
    float z[] = {0.08364854,  0.88834837,  0.5209786 ,  0.4788807 ,  0.44916492,
        1.14548353,  1.31505093,  2.10281808,  1.31999207,  1.59172872,
        1.60872908,  2.53123044,  2.62516707,  2.93745329,  2.64080227,
        3.30207012,  3.24031435,  3.4099925 ,  3.64368553,  3.58203695,
        3.30403619,  3.2671478 ,  3.4561074 ,  3.26540457,  4.10770226,
        3.60773753,  3.11368385,  4.20237017,  3.69169733,  2.96158022,
        3.75248192,  2.90284984,  3.93216292,  3.42890529,  3.51024235,
        3.13083032,  3.61558128,  3.99142341,  4.04422647,  4.40410901,
        4.225771  ,  3.69006258,  4.56701045,  4.11806354,  4.37361484,
        4.30745333,  4.37771531,  4.47956987,  3.87786397,  4.41085414,
        4.50170343,  3.65011336,  4.55706284,  4.08102819,  3.33251889,
        3.61909431,  3.82055851,  3.0910765 ,  2.99813151,  3.48032296,
        2.02432559,  2.6156058 ,  3.38007839,  2.62722002,  2.40901092,
        2.4151884 ,  1.96696851,  1.62203218,  1.86412751,  1.90567553,
        0.84743063,  1.53660705,  0.96366293,  0.80398781,  0.72840369,
        0.69244286,  0.76189612,  0.2738613 ,  0.40470885, 0.};
    part_conv_t *pc = part_conv_new(M,N_ir);
    if (!pc) { return 0; }
    part_conv_set_ir_td(pc,ir);
    size_t WTF, L = N_ir + N_x - 1;
    WTF = L;
    float *x_ = x;
    size_t l = 0;
    while (l < L) { part_conv_proc(pc,x_); x_+=M; l+=M; }
    return check_close(x,z,L);
}


int main (void)
{
    TESTCHK(part_conv_correct_test(),"Test correct:");
    TESTCHK(part_conv_simple_test(),"Test simple:");
    TESTCHK(part_conv_correct_test_2(),"Test correct 2:");
    return 0;
}
