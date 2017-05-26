gcc test/part_conv_profile.c part_conv.c -o test/part_conv_profile.bin -g \
    -DFFT_LIB_VDSP -I. -framework Accelerate
