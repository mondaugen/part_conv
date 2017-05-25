#ifndef PART_CONV_H
#define PART_CONV_H 

#include <stddef.h> 
#include <stdlib.h> 

struct part_conv_t;
typedef struct part_conv_t part_conv_t;

/* Redefine to allocate otherwise. Should behave like calloc (mustn't
   necessarily give allocation aligned to particular address). TODO: These could
   be implemented using the "weak" attribute instead. */
#ifndef PART_CONV_CALLOC
#define PART_CONV_CALLOC(c,s) calloc(c,s) 
#endif  

/* Must behave like free */
#ifndef PART_CONV_FREE
#define PART_CONV_FREE(x) free(x) 
#endif  

part_conv_t *part_conv_new(size_t M,
                           size_t N_ir,
                           size_t D);
void part_conv_free(part_conv_t *pc);
void part_conv_set_ir_td (part_conv_t *pc,
                          const float *ir);
void part_conv_proc(part_conv_t *pc, float *x);

#endif /* PART_CONV_H */
