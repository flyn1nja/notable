#define MaxiFft     true
#define DjFft       false
#define FFTW        false


#include <complex>
#include <math.h>
#include <vector>

#if MaxiFft

#include "maximilian.h"
#include "maximilianRecording.h"
#include "libs/maxim.h"

#elif DjFft

#ifndef FFT_FUNC
#define FFT_FUNC
#endif

#include "dj_fft.h"

#elif FFTW

//#include ...

#endif

namespace FastFourierTransform
{

    template <typename T>
    std::vector<std::complex<T>> CombineToComplexVect(T* reals, T* imags, int len)
    {
        std::vector<std::complex<T>> outVec(len);

        for (int i = 0; i < len; ++i)
            outVec[i] = std::complex<T>(reals[i], imags[i]);

        return outVec;        
    }

#if MaxiFft

maxiFFT FFTres;
maxiIFFT inverseFFTres;
std::vector<float> mags = std::vector<float>(512);
std::vector<float> mags2 = std::vector<float>(512);
std::vector<float> phases = std::vector<float>(512);
std::vector<float> phases2 = std::vector<float>(512);
maxiEnvGen shiftEnvr;
    
    template <typename T>
    std::vector<std::complex<float>> DoFFT(std::vector<T>& data)
    {
        maxiFFT FFTres = maxiFFT();
        for (T& dataVal : data)
        {
            FFTres.process(dataVal);
            /*
            if (FFTres.process(dataVal)) 
            {
                // std::cout << "SC: " << FFTres.spectralCentroid() << endl;
                
                //shift some bins and phases around
                mags = FFTres.getMagnitudes();
                phases = FFTres.getMagnitudes();
                for (int b=0; b < mags.size(); ++b) {
                    int binIdx = b-static_cast<int>(shiftEnvr.play(1));
                    if (binIdx > mags.size()) binIdx = mags.size();
                    if (binIdx < 0 || binIdx >= mags.size())
                    {
                        mags2[b] = 0;
                        phases2[b]=0;
                    }
                    else
                    {
                        mags2[b] = mags[binIdx];
                        phases2[b] = phases[binIdx];
                    }
                }        
            }
            */
        }
        return CombineToComplexVect(FFTres.getReal(), FFTres.getImag(), FFTres.getFFTSize());
    }

    
    // REALLY NOT IDEAL
    template <typename T>
    std::vector<std::complex<float>> DoFFT(const std::vector<std::complex<T>>& data, bool)
    {

        maxiFFT FFTResReal;
        maxiFFT FFTResImag;

        std::vector<std::complex<float>> res;
/*
        FFTres.buffer.clear();
        FFTres.buffer.insert(FFTres.buffer.begin(), data.begin(), --data.end());

*/
        for (int i = 0; i < data.size(); ++i)
        {
            FFTResReal.process(data[i].real());
            FFTResImag.process(data[i].imag());
        }

        // ! NOT WORKING BECAUSE MAXIMILIAN CANT PROCESS FFT FROM COMPLEX
        // *https://www.dsprelated.com/showarticle/97.php

        std::vector<std::complex<float>> resR = CombineToComplexVect(FFTResReal.getReal(), FFTResReal.getImag(), FFTResReal.getFFTSize());
        std::vector<std::complex<float>> resI = CombineToComplexVect(FFTResImag.getReal(), FFTResImag.getImag(), FFTResImag.getFFTSize());

        
        std::vector<std::complex<float>>(FFTResReal.getReal() - FFTResImag.getImag());

        // Xc(m) = real[Xr(m)] - imag[Xi(m)] + j{imag[Xr(m)] + real[Xi(m)]} 
        for (int n = 0; n <  FFTResReal.getFFTSize(); ++n)
            res.push_back(std::complex<float>(resR[n].real() - resI[n].imag(),
                                          resR[n].imag() + resI[n].real()));
        


        return res;
    }

    template <typename T>
    std::vector<std::complex<float>> DoFFTRealTime(std::vector<T>& data)
    {
        if (FFTres.process(data)) 
        {
            // std::cout << "SC: " << FFTres.spectralCentroid() << endl;
            
            //shift some bins and phases around
            mags = FFTres.getMagnitudes();
            phases = FFTres.getMagnitudes();
            for (int b=0; b < mags.size(); ++b) {
                int binIdx = b-static_cast<int>(shiftEnvr.play(1));
                if (binIdx > mags.size()) binIdx = mags.size();
                if (binIdx < 0 || binIdx >= mags.size())
                {
                    mags2[b] = 0;
                    phases2[b]=0;
                }
                else
                {
                    mags2[b] = mags[binIdx];
                    phases2[b] = phases[binIdx];
                }
            }        
        }

        return CombineToComplexVect(FFTres.getReal(), FFTres.getImag(), FFTres.getFFTSize());
    }


#elif DjFft

    template <typename T>
    void DoFFT(std::vector<T>& data)
    {
        //...
    }

#elif FFTW

    template <typename T>
    void DoFFT(std::vector<T>& data)
    {
        //...
    }

#else
// NO FFT DEFINED!!
#endif


} // namespace FastFourierTransform
