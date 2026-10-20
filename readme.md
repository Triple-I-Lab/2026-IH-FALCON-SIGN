Overview

The Fast Fourier Transform (FFT) and fast Fourier sampling (ffSampling) are among the most challenging components of \textsmaller{FALCON} to implement efficiently because of their intensive floating-point arithmetic and recursive computational structure. This project provides source code for FFT/IFFT computation in \textsmaller{FALCON} and our proposed efficient ffSampling algorithm. The proposed algorithm uses a finite-state machine (FSM) to control the sampling process, providing a more implementation-friendly alternative to the recursive algorithm specified in \textsmaller{FALCON}, particularly for hardware implementation.

The implementation was developed and verified based on the \textsmaller{FALCON} team's open-source code, available at \url{https://falcon-sign.info/}.

Source Files
    fft.cpp, fft.h: Adapted from the original \texttt{fft.c} provided by the \textsmaller{FALCON} team and modified to support standalone execution. These files provide two FFT and two IFFT functions corresponding to the two polynomial degrees, $n=512$ and $n=1024$, used by \textsmaller{FALCON}.
    \texttt{radix2.cpp}, \texttt{radix2.h}: Our optimized radix-2 FFT/IFFT implementation for \textsmaller{FALCON}.
    \texttt{ffSampling.cpp}, \texttt{ffSampling.h}: Our FSM-based implementation of the ffSampling algorithm for the \textsmaller{FALCON} security level I parameter set.

FFT Representation and Hermitian Symmetry

The FFT of a polynomial $f$ evaluates $f(\zeta)$ at all roots $\zeta$ of the polynomial $\phi$ over $\mathbb{C}$. Since $\phi$ has distinct complex roots, the FFT produces $n$ complex values, each consisting of a real part and an imaginary part. Storing all these values requires $2n$ floating-point numbers.

However, when $f$ has real coefficients, its FFT representation exhibits Hermitian symmetry. For every root $\zeta$ of $\phi$, its complex conjugate $\bar{\zeta}$ is also a root, and
    $\bar{f(\zeta)}$ = $f(\bar{\zeta})$.
Consequently, the FFT representation contains redundant information, allowing half of the complex values to be omitted. The representation can therefore be stored using $n/2$ complex values, equivalent to $n$ floating-point numbers.

The Hermitian adjoint of $f$ is obtained in the FFT representation by taking the complex conjugate of each value, which negates its imaginary part. Therefore, when a polynomial is equal to its Hermitian adjoint, such as $ff^* + gg^*$, its FFT representation contains only real values.

\textsmaller{FALCON} exploits this property because its input polynomials have real coefficients. The second half of the spectrum is determined by the complex conjugates of the first half, so only $n/2$ complex values need to be computed and stored. This reduces the storage requirements of the FFT representation by half and can also reduce the associated computational cost. This property is described in the \textsmaller{FALCON} specification and reflected in the original \texttt{fft.c} implementation provided by the \textsmaller{FALCON} team.

\subsection{FSM-Based ffSampling Algorithm}

Our ffSampling implementation is organized into four algorithmic states, executed in the following order: \texttt{splitfft}, \texttt{SamplerZ}, calculation of the corresponding $t_0$ from $t_1$, and \texttt{mergefft}, as implemented in \texttt{ffSampling.cpp}.

The recursive computation flow is managed using three variables:
    \texttt{state} controls the current FSM operation.
    $n$ represents the recursion level of $\mathsf{ffSampling}_n$.
    \texttt{idx} tracks the position of the current subproblem in the recursion tree.

This FSM-based organization replaces recursive function calls with explicit state and index management, making the algorithm easier to trace and more suitable for hardware implementation.

Hardware Implementation

In addition to the software implementation, our team has developed a configurable \textsmaller{FALCON} signature-generation architecture in both C and RTL Verilog. The corresponding source code will be made available in this repository upon publication of the associated paper.

Contributions, feedback, and technical discussions are welcome. Please feel free to contact us at \texttt{dangquang@inha.edu}.
