#include <math.h>
#include <stdint.h>
#include <string.h>
#include "fft.h"


void radix2_1024(fpr f[1024]);
void ir2_1024(fpr f[1024]);

void radix2_512(fpr f[512]);
void ir2_512(fpr f[512]);
void print_tf2();
