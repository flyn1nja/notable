#pragma once
#ifndef FFTransform
#define FFTransform

#define MaxiFft     true
#define DjFft       false
#define FFTW        false

#if MaxiFft

#include "libs/maxim.h"
#include "maximilian.h"
#include "maximilianRecording.h"

#elif DjFft

#ifndef FFT_FUNC
#define FFT_FUNC
#endif

#include "dj_fft.h"

#elif FFTW

//#include ...

#endif

#include <complex>
#include <vector>

namespace FastFourierTransform
{

#if MaxiFft
    maxiFFT FFTres;
    maxiIFFT inverseFFTres;
    std::vector<float> mags = std::vector<float>(512);
    std::vector<float> mags2 = std::vector<float>(512);
    std::vector<float> phases = std::vector<float>(512);
    std::vector<float> phases2 = std::vector<float>(512);
    maxiEnvGen shiftEnvr;
#endif

    template <typename T>
    std::vector<std::complex<T>> CombineToComplexVect(T* reals, T* imags, int len);

    template <typename T>
    std::vector<std::complex<float>> DoFFT(std::vector<T>& data);

    template <typename T>
    std::vector<std::complex<float>> DoFFT(const std::vector<std::complex<T>>& data, bool);

    template <typename T>
    std::vector<std::complex<float>> DoFFTRealTime(std::vector<T>& data);
}
#endif