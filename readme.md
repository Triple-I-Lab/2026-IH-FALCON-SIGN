Fast Fourier Transform (FFT) and Fast Fourier Sampling (ffSampling) is most challenging in Falcon for implementation, due to the hard in floating-point arithmetic and recursive function. This project presents the source code for a FFT/IFFT computations in Falcon and the proposed efficient ffSampling algorithm, which base on finite-state-machine (FSM) to control, more friendly with user and hardware implementations than the recursive algorithm proposed in Falcon specification.  . The implementation is developed and verified based on the Falcon team’s open-source code, available at: https://falcon-sign.info/

fft.cpp, fft.h: Adapted from the original fft.c provided by the Falcon team and modified to operate as a standalone implementation. These files include two FFT and two IFFT functions corresponding to the two parameter sets of $n$ in Falcon.

radix2.cpp, radix2.h: Our optimized radix-2 FFT/IFFT implementation for Falcon

radix2.cpp, radix2.h: Our optimized radix-2 FFT/IFFT implementation for Falcon.
ffSampling.cpp, ffSampling.h: The code design for FSM-based ffSampling algorithm in Falcon level I parameter sets.


The Fast Fourier Transform for a polynomial computes $f$($\zeta$) for all roots of $\phi$ (over C). $\phi$ has distinct roots in C, leading to $n$ values $f(\zeta)$, each being a complex number, with a real and an imaginary part. Storage space requirements are then 2$n$ floating-point numbers. However, if $f$ is real, then, for every root $\zeta$ of $\phi$, $\bar{\zeta}$ is also a root of $\phi$ and $\bar{f(\zeta)}$ = $f(\bar{\zeta})$. Thus, the FFT representation is redundant, and half of the values can be omitted, reducing storage space requirements to $n$/2 complex numbers, hence $n$ floating-point values.

The Hermitian adjoint of $f$ is obtained in FFT representation by simply computing the conjugate
of each $f(\zeta)$, i.e. negating the imaginary part. This means that when a polynomial is equal to
its Hermitian adjoint (e.g. $ff^* + gg^*$), then its FFT representation contains only real values.

In Falcon, since the input polynomials are real, this Hermitian property is directly exploited. The second half of the spectrum is simply the complex conjugate of the first half, so only $n$/2 complex values must be computed and stored. This reduces both computation and storage requirements by approximately half compared to a conventional complex FFT. This property is clearly described in the Falcon specification and reflected in the original fft.c implementation provided by the Falcon team. 

Our ffSampling is divide into 4 algorithmic states: splitfft, SamplerZ, Calculate corresponding $t_0$ from $t_1$, mergefft as in order present in ffSampling.cpp. The recursive computation flow is tracked by 4 variables: $state$ controls current FSM operation, $n$ represents recursion level of $\mathsf{ffSampling}_n$, $idx$ tracks position of the current subproblem in the recursion tree. 

Besides the software code, the implementation of configurable Falcon signature generation has been developed by our team detail both in C code and RTL Verilog code, which will be available in here together with the time of paper is publish. We warmly welcome your contribution and please feel free to contact us for discussion (email: dangquang@inha.edu).
