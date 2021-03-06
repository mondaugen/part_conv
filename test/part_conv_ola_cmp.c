/* Compare conv_ola and part_conv */

#include <stdlib.h>
#include <time.h> 

#include "part_conv.h"
#include "conv_ola.h"

#if SAMP_FLOAT_NBITS!=32
#error "Only works with 32 bit floats."
#endif  

#define ERR_TOL 1.e-4

#ifndef  FABS
#define FABS(x) (((x)<0.)?((-1)*(x)):(x))
#endif  

#define TESTCHK(x,s) \
    if (x) { puts(s " passed."); } else { puts(s " failed."); }

int check_close(float *x, float *y, size_t n, float *maxe, size_t *maxen)
{
    int ret = 1;
    while (n--) { ret = ret && (FABS(x[n]-y[n]) < ERR_TOL);
        *maxe = FABS(x[n]-y[n]) > *maxe ? FABS(x[n]-y[n]) : *maxe;
        *maxen = FABS(x[n]-y[n]) == *maxe ? n : *maxen; }
    return ret;
}

int ConvOLA_test(void)
{
    /* FFT size */
//    size_t N = 128, H = 8, n, N_x = 100000, N_x_, N_=256;
//    size_t N = 16, H = 8, n, N_x = 100000, N_x_, N_=32;
    size_t N = 8192, H = 64, n, N_x = 100000, N_x_, N_=16384;
    N_x_ = (N_x/H + 1)*H;

    srandom(time(NULL));

    float *ir_td = (float*)malloc(sizeof(float)*N);
    for(n=0;n<N;n++){ir_td[n]=2.*random()/(float)RAND_MAX - 1.;}
    float *x_td = (float*)calloc(N_x_,sizeof(float));
    for(n=0;n<N_x;n++){x_td[n]=2.*random()/(float)RAND_MAX - 1.;}

   /* split frequency-domain representation of IR */
    DSPSplitComplex ir = {
        (samp_float_t*)malloc(sizeof(samp_float_t)*(N_/2+1)),
        (samp_float_t*)malloc(sizeof(samp_float_t)*(N_/2+1))
    };
    vDSP_vclr(ir.realp,1,N_/2+1);
    vDSP_vclr(ir.imagp,1,N_/2+1);
    vDSP_DFT_Setup ir_dft_plan = vDSP_DFT_zrop_CreateSetup(NULL,N_,vDSP_DFT_FORWARD);
    vDSP_ctoz((DSPComplex*)ir_td,2,&ir,1,N/2);
    /* Scale down because IFFT multiplies by N_ and FFT multiplies by 2 (for each
     * transformed input, both of which are multiplied together, therefore scale
     * by 1/4)*/
    samp_float_t _scalar = 1./(4.*N_);
    vDSP_vsmul(ir.realp,1,&_scalar,ir.realp,1,N_/2);
    vDSP_vsmul(ir.imagp,1,&_scalar,ir.imagp,1,N_/2);
    vDSP_DFT_Execute(ir_dft_plan,ir.realp,ir.imagp,ir.realp,ir.imagp);
    /* extract Nyquist */
    ir.realp[N_/2] = ir.imagp[0];
    ir.imagp[0] = 0.;

    /* Set up convolution system */
    ConvOLADesc colad;

    /* Set up for real now */
    ConvOLADesc_setup(&colad, N_, H);

    colad.ir = ir;

    /* Allocate space for output */
    float *out_c = (float*)calloc(N_x_,sizeof(float));

    /* Copy input to output */
    memcpy(out_c, x_td, sizeof(float)*N_x);

    /* Do processing */
    size_t h = 0;
    while (h < N_x_) {
        ConvOLA_tick(out_c + h, &colad);
        h += H;
    }

    ConvOLADesc_destroy(&colad);
    
    /* Now compare with partition convolver */

    part_conv_t *pc = part_conv_new(H,N);
    if (!pc) { fprintf(stderr,"Error allocating part_conv.\n"); return 0; }
    part_conv_set_ir_td(pc,ir_td);

    /* Allocate space for output */
    float *out_p = (float*)calloc(N_x_,sizeof(float));

    /* Copy input to output */
    memcpy(out_p, x_td, sizeof(float)*N_x);

    h = 0;
    while (h < N_x_) { part_conv_proc(pc,out_p+h); h += H; }
    float maxe = 0; size_t maxen = N_x_; int ret;
    ret = check_close(out_c,out_p,N_x_,&maxe,&maxen);
    printf("N_x_: %zu\n",N_x_);
    printf("max error: %g, index %zu\n",maxe,maxen);
    return ret;
}

int main (void)
{
    TESTCHK(ConvOLA_test(),"conv_ola and part_conv compare:");
}
