#include "part_conv.h"
#include <math.h> 
#ifdef FFT_LIB_VDSP
#include <Accelerate/Accelerate.h> 
#endif  

#define VALIGN_BYTES 64 

typedef struct f32_part_t {
    struct f32_part_t *next;
    float *f;
} f32_part_t;

typedef struct z32_part_t {
    struct z32_part_t *next;
#ifdef FFT_LIB_VDSP
    DSPSplitComplex z;
#endif  
} z32_part_t;

struct part_conv_t {
    /* Size of incoming processed vectors */
    size_t M; 
    /* Size of IR to convolve incoming vectors with. */
    size_t N_ir; 
    /* Total size of convolution buffer rounded to next greater
       power of 2 */
    size_t N_c; 
    /* Number of times M fits into N_c */
    size_t K;
    /* Number of parts that IR has been partitioned into */
    size_t D;
    /* Previous convolutions (time domain) */
    f32_part_t *convs;
    /* Parts of IR (frequency domain, updated when IR changes) */
    z32_part_t *ir_parts;
    /* Previous input blocks (frequency domain) */
    z32_part_t *past_ins;
#ifdef FFT_LIB_VDSP
    /* Temporary processing buffer */
    DSPSplitComplex convt;
    /* DFT plans */
    vDSP_DFT_Setup dft_setup_fw;
    vDSP_DFT_Setup dft_setup_bw;
#endif  
};

static int is_pow_2(size_t x)
{
    int bitidx = (int)(sizeof(size_t) * 8 - 1);
    while (((x & (1 << bitidx)) == 0) && (bitidx >= 0)) {
        bitidx--;
    }
    if (bitidx == -1) {
        /* All bits were 0. 0 is not a power of 2. */
        return 0;
    }
    bitidx--;
    while (((x & (1 << bitidx)) == 0) && (bitidx >= 0)) {
        bitidx--;
    }
    if (bitidx == -1) {
        /* All bits after the initial set bit were 0, is a power of 2. */
        return 1;
    }
    /* Some bit was set after the initial set bit. Not power of 2. */
    return 0;
}

/* Return the next power of 2 greater than or equal to x */
static int next_pow_2(float x)
{
    return (int)pow(2.,ceil(log(x)/log(2.)));
}

size_t part_conv_opt_D(size_t M, size_t N_ir) {
    /* TODO */
    return 1;
}

/* Create new partition convolver that can process vectors of length M and
   convolves with IRs of length N_ir partitioned into parts of size N_ir/D. M
   must divide next_pow_2(N_ir/D+M-1). N_ir must be divisible by D. The best
   value of D can be estimated using part_conv_opt_D. */
part_conv_t *part_conv_new(size_t M,
                           size_t N_ir,
                           size_t D)
{
    if (N_ir % D) { return NULL; }
    size_t N_c = (int)next_pow_2((float)N_ir/D + M - 1);
    if (N_c % M) { return NULL; }
    size_t K = N_c / M;
    /* Need to store previous convolutions = K*N_c floats,
       Parts of IR = D*(N_c/2+1) complex,
       Previous transformed inputs = D*(N_c/2+1) complex,
       and a temporary buffer for computations = N_c/2+1 complex */
    /* Pad floats to multiple of VALIGN_BYTES */
    size_t conv_sz = (sizeof(float)*N_c/VALIGN_BYTES
            + (((sizeof(float)*N_c)%VALIGN_BYTES) != 0))*VALIGN_BYTES;
#ifdef FFT_LIB_VDSP
    /* For split complex type, we ensure real and imaginary have proper alignment */
    size_t irs_sz = 2*((sizeof(float)*N_c/2+1)/VALIGN_BYTES
            + (((sizeof(float)*N_c/2+1)%VALIGN_BYTES) != 0)*VALIGN_BYTES);
    /* Size of data structures to get to start of data segment */
    size_t head_sz = sizeof(part_conv_t) + K*sizeof(f32_part_t)
        + 2*D*sizeof(z32_part_t);
    /* Allocate space for part_conv_t, its parts and their data and pad so
       proper alignment can be enforced. */
    size_t ms = sizeof(part_conv_t) + K*sizeof(f32_part_t)
        + 2*D*sizeof(z32_part_t) + K*conv_sz
        + (2*D+1)*irs_sz + VALIGN_BYTES - 1;
#endif  
    void *mem = PART_CONV_CALLOC(ms,sizeof(char)); if (!mem) { return NULL; }
    void *dat = mem + head_sz + VALIGN_BYTES 
        - ((size_t)mem + head_sz) % VALIGN_BYTES;
    /* Assign data structures to parts of allocated memory */
    part_conv_t *pc = (part_conv_t*)mem;
    *pc = (part_conv_t) {
        .M = M, 
        .N_ir = N_ir, 
        .N_c = N_c, 
        .K = K,
        .D = D,
    };
    mem += sizeof(part_conv_t);
    f32_part_t *first_conv = NULL;
    while (K--) {
        f32_part_t *tmp = pc->convs;
        pc->convs = mem; pc->convs->next = tmp;
        first_conv = first_conv ? first_conv : pc->convs;
        mem += sizeof(f32_part_t);
        pc->convs->f = dat; dat += conv_sz;
    }
    /* Make into ring */
    first_conv->next = pc->convs;
    /* Only past inputs (not ir_parts) need be ring */
    z32_part_t *first_past_in;
    /* Paranoia to keep list iteration from advancing into unknown memory */
    pc->ir_parts = NULL; 
    while (D--) {
        z32_part_t *tmp = pc->ir_parts;
        pc->ir_parts = mem; pc->ir_parts->next = tmp;
        mem += sizeof(z32_part_t);
        tmp = pc->past_ins; 
        pc->past_ins = mem; pc->past_ins->next = tmp;
        first_past_in = first_past_in ? first_past_in : pc->past_ins;
        mem += sizeof(z32_part_t);
#ifdef  FFT_LIB_VDSP 
        pc->ir_parts->z.realp = dat; dat += irs_sz/2;
        pc->ir_parts->z.imagp = dat; dat += irs_sz/2;
        pc->past_ins->z.realp = dat; dat += irs_sz/2;
        pc->past_ins->z.imagp = dat; dat += irs_sz/2;
#endif  
    }
    /* make into ring */
    first_past_in->next = pc->past_ins;
#ifdef  FFT_LIB_VDSP 
    /* Assign memory for temporary DFT processing buffer */
    pc->convt.realp = dat; dat += irs_sz/2;
    pc->convt.imagp = dat; dat += irs_sz/2;
    /* Allocate DFT setup for forward and inverse */
    pc->dft_setup_fw = vDSP_DFT_zrop_CreateSetup(0,N_c,vDSP_DFT_FORWARD);
    pc->dft_setup_bw = vDSP_DFT_zrop_CreateSetup(0,N_c,vDSP_DFT_INVERSE);
#endif  
    return pc;
}

void part_conv_free(part_conv_t *pc)
{
    vDSP_DFT_DestroySetup(pc->dft_setup_fw);
    vDSP_DFT_DestroySetup(pc->dft_setup_bw);
    free(pc);
}

/* Set the IR from a time domain representation.
   ir must have length pc->N_ir.
   TODO: Eventually there will also be a way to set using the frequency-domain paritions.  */
void part_conv_set_ir_td (part_conv_t *pc,
                          const float *ir)
{
    size_t D = pc->D;
    z32_part_t *tmp = pc->ir_parts;
    while (D--) {
#ifdef FFT_LIB_VDSP
        vDSP_vclr(tmp->z.realp,1,pc->N_c/2 + 1); 
        vDSP_vclr(tmp->z.imagp,1,pc->N_c/2 + 1); 
        /* Stride of 2 because of vDSP legacy style */
        vDSP_ctoz((DSPComplex*)ir,2,&tmp->z,1,pc->N_ir/pc->D/2);
        /* Do forward DFT */
        vDSP_DFT_Execute(pc->dft_setup_fw,tmp->z.realp,tmp->z.imagp,tmp->z.realp,tmp->z.imagp);
        /* Scale by 1/(4*N_c) because it will get multiplied by other DFT'd input
           also scaled by 2 and will be scaled by N_c on inverse transform. */
        float scale = 1/(4.*pc->N_c);
        vDSP_vsmul(tmp->z.realp,1,&scale,tmp->z.realp,1,pc->N_c/2);
        vDSP_vsmul(tmp->z.imagp,1,&scale,tmp->z.imagp,1,pc->N_c/2);
        /* Extract nyquist from DC bin (will have to be transferred back for inverse transform */
        tmp->z.realp[pc->N_c/2] = tmp->z.imagp[0]; tmp->z.imagp[0] = 0;
#endif 
        tmp = tmp->next;
        ir += pc->N_ir/pc->D;
    }
}

/* Process a vector of input with the partition convolve algorithm. x must have
   length pc->M. */
void part_conv_proc(part_conv_t *pc, float *x)
{
#ifdef FFT_LIB_VDSP
    /* Clear next past_ins */
    vDSP_vclr(pc->past_ins->z.realp,1,pc->N_c/2+1); 
    vDSP_vclr(pc->past_ins->z.imagp,1,pc->N_c/2+1); 
    /* Read in input as split complex */
    vDSP_ctoz((DSPComplex*)x,2,&pc->past_ins->z,1,pc->M/2);
    /* Do forward DFT in place */
    vDSP_DFT_Execute(pc->dft_setup_fw,
            pc->past_ins->z.realp,
            pc->past_ins->z.imagp,
            pc->past_ins->z.realp,
            pc->past_ins->z.imagp);
    /* Extract nyquist from DC bin */
    pc->past_ins->z.realp[pc->N_c/2] = pc->past_ins->z.imagp[0]; pc->past_ins->z.imagp[0] = 0;
    /* Sum in this and past inputs convolved (multiplied) with repective IR
       parts. First zero block of temporary convolution buffer.*/
    vDSP_vclr(pc->convt.imagp,1,pc->N_c/2+1);
    vDSP_vclr(pc->convt.realp,1,pc->N_c/2+1);
    size_t D = pc->D;
    z32_part_t *irs = pc->ir_parts, *past = pc->past_ins;
    while (D--) {
        /* Do multiply and add using IRs into pc->convt. Must multiply
           pc->N_c/2+1 elements because Nyquist has been extracted. */
        vDSP_zvma(&irs->z,1,&past->z,1,&pc->convt,1,&pc->convt,1,pc->N_c/2+1);
        irs = irs->next;
        /* Keep track of past to move this pointer to the next available buffer
           for writing. */
        pc->past_ins = past;
        past = past->next;
    }
    /* Put nyquist back to DC */
    pc->convt.imagp[0] = pc->convt.realp[pc->N_c/2]; pc->convt.realp[pc->N_c/2] = 0;
    /* Backward transform in place */
    vDSP_DFT_Execute(pc->dft_setup_bw,
            pc->convt.realp,
            pc->convt.imagp,
            pc->convt.realp,
            pc->convt.imagp);
    /* Deinterleave split format into current convolution buffer */
    vDSP_ztoc(&pc->convt,1,(DSPComplex*)pc->convs->f,2,pc->N_c/2);
    /* Accumulate past convolutions into output */
    size_t K = pc->K, x_n = 0;
    f32_part_t *conv = pc->convs;
    vDSP_vclr(x,1,pc->M);
    while (K--) {
        /* Add convs into x */
        vDSP_vadd(x,1,conv->f+x_n,1,x,1,pc->M);
        pc->convs = conv; x_n += pc->M;
        conv = conv->next;
    }
#endif
    /* Pointers already advanced */
}



