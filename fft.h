#ifndef FFT_H
#define FFT_H

#include <math.h>
#include <stdint.h>
#include <string.h>



typedef struct { double v; } fpr;

static inline fpr
FPR(double v)
{
	fpr x;

	x.v = v;
	return x;
}

static inline fpr
fpr_add(fpr x, fpr y)
{
	return FPR(x.v + y.v);
}

static inline fpr
fpr_sub(fpr x, fpr y)
{
	return FPR(x.v - y.v);
}

static inline fpr
fpr_mul(fpr x, fpr y)
{
	return FPR(x.v * y.v);
}

static inline fpr
fpr_neg(fpr x)
{
	return FPR(-x.v);
}

#define FPC_ADD(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_re, fpct_im; \
		fpct_re = fpr_add(a_re, b_re); \
		fpct_im = fpr_add(a_im, b_im); \
		(d_re) = fpct_re; \
		(d_im) = fpct_im; \
	} while (0)

#define FPC_SUB(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_re, fpct_im; \
		fpct_re = fpr_sub(a_re, b_re); \
		fpct_im = fpr_sub(a_im, b_im); \
		(d_re) = fpct_re; \
		(d_im) = fpct_im; \
	} while (0)

#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_b_re, fpct_b_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_b_re = (b_re); \
		fpct_b_im = (b_im); \
		fpct_d_re = fpr_sub( \
			fpr_mul(fpct_a_re, fpct_b_re), \
			fpr_mul(fpct_a_im, fpct_b_im)); \
		fpct_d_im = fpr_add( \
			fpr_mul(fpct_a_re, fpct_b_im), \
			fpr_mul(fpct_a_im, fpct_b_re)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

void fft1024(fpr f[1024]);
void fft8(fpr f[8]);
void ifft1024(fpr f[1024]);
void multconvolution(fpr h[8], fpr f[8], fpr g[8]);

void r2fft(fpr f[512]);
void r2ifft(fpr f[512]);
void ifft8(fpr f[8]);

#endif /* FFT_H */