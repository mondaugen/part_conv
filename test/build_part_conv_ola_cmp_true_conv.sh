#/bin/bash
rm -f test/part_conv_ola_cmp_true_conv.bin
gcc -DSAMP_FLOAT_NBITS=32 part_conv.c test/part_conv_ola_cmp_true_conv.c \
    ../soundalivelib/src/conv_ola.c \
    ../soundalivelib/src/emalloc.c \
    ../soundalivelib/src/sa_utils.c \
    -DFFT_LIB_VDSP -DVMATH_VDSP -g -framework Accelerate \
    -I. -I../soundalivelib/inc -o test/part_conv_ola_cmp_true_conv.bin
