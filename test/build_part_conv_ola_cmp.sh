#/bin/bash
rm -f test/part_conv_ola_cmp.bin
gcc -DSAMP_FLOAT_NBITS=32 part_conv.c test/part_conv_ola_cmp.c \
    ../soundalivelib/src/conv_ola.c \
    ../soundalivelib/src/emalloc.c \
    ../soundalivelib/src/sa_utils.c \
    -DFFT_LIB_VDSP -g -framework Accelerate \
    -I. -I../soundalivelib/inc -o test/part_conv_ola_cmp.bin
