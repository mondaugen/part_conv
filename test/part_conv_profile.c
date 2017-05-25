/* Firstly, is the algorithm correct?
   For which value of D is the algorithm fastest? */

#include "part_conv.h"

#include <stdio.h> 

#define ERR_TOL 1.e-6 

//#ifndef  FABS
#define FABS(x) ((x)<0.)?((-1)*(x)):(x) 
//#endif  

int part_conv_correct_test (void)
{
    float in[] = {1.,0,.5,0,0,0,0,0},
          ir[] = {0.5, 0.25, -0.5, -0.25},
          outc[] = {0.5, 0.25, -0.5+0.25, -0.25+0.125,-0.25,-0.125,0.,0.};
    part_conv_t *pc = part_conv_new(2,4,2); if (!pc) { return -1; }
    part_conv_set_ir_td(pc,ir);
    part_conv_proc(pc,in);
    part_conv_proc(pc,in+2);
    part_conv_proc(pc,in+4);
    part_conv_proc(pc,in+6);
    int ret = 1;
    size_t n;
    for (n=0;n<8;n++) { ret = ret && (FABS(in[n]-outc[n]) < ERR_TOL); }
    return ret;
}

int main (void)
{
    if (part_conv_correct_test()) { puts("Passed."); }
    else { puts("Failed."); }
    return 0;
}
