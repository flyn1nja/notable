#pragma once
#ifndef QCONST_TRANSFORM
#define QCONST_TRANSFORM

#include <complex>
#include <math.h>
#include <vector>
#include <tuple>

#ifndef MaxiFft
#define MaxiFft true
#endif


#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif


#if MaxiFft

// #include "maximilian.h"
// #include "maximilianRecording.h"
// #include "libs/maxim.h"

class maxiFFT;
class maxiRecorder;

#endif

typedef float number;

typedef std::vector<number> Vector;
typedef std::vector<std::complex<number>> ComplexVector;
typedef std::vector<std::vector<number>> Matrix;
typedef std::vector<std::vector<std::complex<number>>> ComplexMatrix;

struct QConstTransKernel
{
    std::complex<number>** coeffCQT; // The Q-constant spectral coefficients in the form of a sparse matrix. Rasterized
    number freqKernel; // Spectral Kernel
    number fMin; // Minimum frequency (freq of the lowest bin)
    number fMax; // Maximum frequency (freq of the highest bin)
    int octaveCount; // Count of the number of octaves
    int binPerOct; // Count of the bin per octave
    int params; // Params for the inverse transform
};

template <typename T>
struct sparseMat {
    // using matValType = std::conditional_t<std::is_arithmetic<T>::value>;
    std::vector<std::tuple<T, int, int>> values;
    int width = 0;
    int height = 0;
    bool transposed = false;

    // void appendValues(std::vector<std::tuple<T, int, int>> newValues)
    // {
    //     this->values.push_back(newValues);
    // }
};

struct CQTKernel
{
    sparseMat<std::complex<number>> freqKernel;

    int fftLen;
    int fftHop;

    number nkMax;
    number nkMin;

    int atomHop;
    int atomNr;
    number firstcenter;
};

struct CQTResult
{
    // // Matrix spCQT;
    sparseMat<std::complex<number>> spCQT;
    // // ComplexMatrix fKernel;
    CQTKernel fKernel;
    
    number fmin;
    number fmax;
    int octaveNr;
    int bins;
    // Add additional fields as needed
};



namespace QConstTrans
{

    const int COEFFS_SIZE = 7;
    // Computed from matlab:  butter(6, 0.5, 'low') 
    const number BUTTER_B_COEFFS[COEFFS_SIZE] = { 0.029588, 0.177529, 0.443823, 0.59174, 0.443823, 0.177529, 0.29588 };
    const number BUTTER_A_COEFFS[COEFFS_SIZE] = { 1.000000, -6.6613e-16, 7.7770e-1, -2.8192e-16, 1.1420e-1, -1.1472e-17, 1.7509e-3 };


    inline int nextpow2(int x) { return round(pow(2.0f, ceil(log2(x)))); }
    inline number sumVect(const std::vector<number>& v)
    {
        number sum = 0;
        for(int i=0; i<v.size(); ++i) sum += v[i];
        return sum;
    }

    ComplexVector conjugate(std::complex<number>* v, int length)
    {
        ComplexVector conjVect(length,0.0);

        for (int i = 0; i < length; ++i)
            conjVect[i] = std::conj(v[i]);
        
        return conjVect;
    }

    inline ComplexVector conjugate(ComplexVector v, int length)
    {
        ComplexVector conjVect(length);

        for (int i = 0; i < length; ++i)
            conjVect[i] = std::conj(v[i]);
        
        return conjVect;
    }

    template <class T>
    inline T** transposeMat(T** mat, int lenx, int leny)
    {
        T copyMat[lenx][leny];
        for (int i = 0; i < leny; ++i)
        for (int j = 0; j < lenx; ++j)
            copyMat[j][i]= mat[i][j];
        
        return copyMat;
    }
    template <class T>
    inline T** transposeMat(sparseMat<T>& sparseM)
    {
        sparseM.transposed = true;
    }

    template <class T>
    inline number meanOfAbs(std::vector<T> vect)
    {
        number sum = static_cast<number>(0); 
        
        for (int i = 0; i < vect.size(); ++i)
            sum += static_cast<number>(abs(vect[i]));
        
        return sum / vect.size();
    }

    template <class T>
    inline number meanOfAbs(std::vector<std::vector<T>> mat)
    {
        number sum = static_cast<number>(0); 
        
        for (int i = 0; i < mat.size(); ++i)
        for (int j = 0; j < mat[i].size(); ++j)
            sum += static_cast<number>(abs(mat[i][j]));
        
        return sum / (mat.size() * mat[0].size());
    }

    inline std::vector<number> modhann(number length)
    {

        number left = 0.5 - floor(length/2.0f)*1/length;
        // number right = 0.5 + floor(length/2.0f)*1/length;
        std::vector<number> output(length+1);
        
        output[0] = left;

        for (int i = 1; i < length+1; ++i)
            output[i] = output[i-1] + 1.0f/length;

        for (int i = 0; i < output.size(); ++i)
            output[i] = 0.5 * (1- cos(2*PI*output[i]));

        return output;
    }

    // Only store positions of non zero values in a matrix.
    // Other values are assumed to be 0 when reconstructing it
    template <class T>
    inline sparseMat<T> sparse(std::vector<std::vector<T>> mat, float thresh=0.0005f)
    {
        sparseMat<T> sparseMat;
        sparseMat.values = std::vector<std::tuple<T, int, int>>();
        sparseMat.height = mat.size();
        sparseMat.width = mat[0].size();

        // if (std::is_arithmetic<T>)
        // {
            for (int i = 0; i < mat.size(); ++i)
            for (int j = 0; j < mat[i].size(); ++j)
            if (abs(mat[i][j]) >= thresh)
                sparseMat.values.push_back(std::make_tuple(mat[i][j], i, j));
            
        // }
        return sparseMat;
    }

    template <class T>
    inline std::vector<std::vector<T>> unsparse(const sparseMat<T>& sparseM)
    {
        // Creates the matrix with 0 as def. value
        std::vector<std::vector<T>> matr = 
                std::vector<std::vector<T>>(sparseM.height,
                std::vector<T>(sparseM.width, static_cast<T>(0)));

            for (int i = 0; i < sparseM.values.size(); ++i)
                matr[std::get<1>(sparseM.values[i])][std::get<2>(sparseM.values[i])] = std::get<0>(sparseM.values[i]);
            
        // }
        return matr;
    }

    template <class T>
    inline void sparseAppend(sparseMat<T>& sparseM, std::vector<T> vec, float thresh=0.0005f)
    {
        sparseM.height = vec.size();
        ++sparseM.width;

        // if (std::is_arithmetic<T>)
        // {
            for (int j = 0; j < vec.size(); ++j)
            if (abs(vec[j]) >= thresh)
                sparseM.values.push_back(std::make_tuple(vec[j], sparseM.width-1, j));
                // sparseM.values.push_back(std::tuple<T, int, int>(vec[j], sparseM.width-1, j));
            
        // }
        // return sparseMat;
    }

    // Function to multiply a transposed matrix by a vector
    template <class T>
    inline std::vector<T> multiplyTransposedMatrixByVec(const sparseMat<T>& spMat, const std::vector<T>& vec)
    {
        std::vector<T> result(vec.size(), 0.0);

        // !! > REVOIR < !!
        for (int i = 0; i < spMat.height; ++i)
        {
            for (int j = 0; j < spMat.values.size(); ++j)
            {
                if (std::get<1>(spMat.values[j]) == i)
                    result[i] += spMat.values[j] * vec[i];
            }
        }

        return result;
    }

    template <class T>
    inline sparseMat<T> bufferVec(const std::vector<T>& vec, int winSize, int hop, number thresh)
    {
        sparseMat<T> result;

        for (int i = 0; i < vec.size(); i += hop)
        for (int j = i; j < i + winSize; ++j)
        {
            if (abs(vec[j]) >= thresh)
                result.values.push_back(std::make_tuple(vec[j], i, j));
        }
        

        return result;
    }


    // Function to multiply two matrices
    template <class T>
    inline sparseMat<T> multiplyMatrices(const sparseMat<T>& spMat1, const sparseMat<T>& spMat2, bool transposedLast=false)
    {
        sparseMat<T> result;

        const int rows1 = spMat1.height;
        const int cols1 = spMat1.width;
        const int rows2 = spMat2.height;
        const int cols2 = spMat2.width;

        T sum = T{ };

        std::tuple<T, int, int> val1;
        std::tuple<T, int, int> val2;

        // Multiply matrix by its transposed form
        for (int i = 0; i != spMat1.values.size(); ++i) {
            val1 = spMat1.values[i];

            for (int j = 0; j != spMat2.values.size(); ++j) {
                val2 = spMat2.values[j];

                if (transposedLast)
                {
                    // if both rows are the same
                    if (std::get<1>(val1) == std::get<1>(val2))
                        sum += std::get<2>(val1) * std::get<2>(val2);  // Multiply each element

                }
                else
                {
                    // if rows and columns are the same
                    if (std::get<1>(val1) == std::get<2>(val2))
                        sum += std::get<2>(val1) * std::get<2>(val2);  // Multiply each element
                }

            }
        }

        return result;
    }

    // Function to multiply matrix by its transposed form and returns diagonal elems
    template <class T>
    inline std::vector<T> diagOfMultiplyMatrixByTranspose(sparseMat<T>& spMat)
    {
        const int rows = spMat.height;
        const int cols = spMat.width;

        std::vector<T> result;
        result.reserve(rows);

        T sum = T{ };
        std::tuple<T, int, int> val1;
        std::tuple<T, int, int> val2;


        // Multiply matrix by its transposed form
        for (int i = 0; i != spMat.values.size(); ++i) {
            val1 = spMat.values[i];

            for (int j = 0; j != spMat.values.size(); ++j) {
                val2 = spMat.values[j];

                // if both columns are the same
                if (std::get<2>(val1) == std::get<2>(val2))
                {
                    sum += std::get<2>(val1) * std::get<2>(val2);  // Square each element
                    break;
                }
            }
        }


        return result;
    }

    inline std::vector<number> antialias(const std::vector<number>& x) {
        // std::vector<number> y1 = forwardFilter(b, x);

        const int sizex = x.size();

        //* Maybe rendre static? Ou declarer a l'exterieur? *
        std::vector<number> y1(sizex, 0.0);
        std::vector<number> y2(sizex, 0.0);
        
        // Do forward filtering
        for (int i = COEFFS_SIZE; i < sizex; ++i)
        for (int j = 0; j < COEFFS_SIZE; ++j)
            y1[i] += BUTTER_B_COEFFS[j] * x[i-j];

        // Do backward filtering
        for (int i = sizex-1; i >= COEFFS_SIZE; --i)
        for (int j = COEFFS_SIZE-1; j >= 0; --j)
            y2[i] += BUTTER_B_COEFFS[j] * y1[i-j];

    /*
        // Apply normalization
        number scale = 0.0;
        for (int i = 0; i < a.size(); ++i) {
            scale += a[i];
        }
        for (int i = 0; i < y.size(); ++i) {
            y[i] /= scale;
        }
    */
        return y2;
    }



    template <class T>
    inline ComplexVector getWK(sparseMat<T>& sparseM, const number q)
    {
        /*
        std::tuple<T, int, int> wx1;
        std::tuple<T, int, int> wx2;
        // Normalize the magnitudes of the atoms
        for (std::tuple<T, int, int> val : sparseM.values)
        {
            if (std::get<2>(val) == 0 &&
                std::get<0>(val) > std::get<0>(wx1))
                wx1 = val;
            
            if (std::get<0>(val) > std::get<0>(wx2))
                wx2 = val;
            
            // [ignore,wx1] = max(sparseKernel(:,1));
            // [ignore,wx2] = max(sparseKernel(:,end));        
        }
        
        sparseMat<T> wk = sparseMat<T>();
        wk.width = abs(std::get<1>(wx2) - std::get<1>(wx1));
        wk.height = sparseM.height;

        // wK = sparseKernel(wx1:wx2,:);
        for (std::tuple<T, int, int> val : sparseM.values)
        {
            int posXVal = std::get<1>(val);
            if (clamp(posXVal,wx1,wx2) == posXVal)
                wk.values.push_back(val);
        }
        */
                // Optimisation -- not done according to max values
        
        // std::vector<T> wKVect = diagOfMultiplyMatrixByTranspose(wk);

        ComplexVector wKVect = diagOfMultiplyMatrixByTranspose(sparseM);

        // wK = wK(round(1/q)+1:(end-round(1/q)-2));
        return ComplexVector((wKVect.begin() + static_cast<std::ptrdiff_t>(1.0f/q + 0.5f)),
                            (wKVect.end() - static_cast<std::ptrdiff_t>(1.0f/q + 0.5f - 3)));
    }

    inline void normalizeByWeight(CQTKernel& spKernel, ComplexVector& wk)
    {
        number weight = sqrt((spKernel.fftHop/static_cast<number>(spKernel.fftLen)) / meanOfAbs(wk));

        for(std::tuple<std::complex<number>, int, int>& elem : spKernel.freqKernel.values)
            elem = std::make_tuple(std::get<0>(elem) * weight, std::get<1>(elem), std::get<2>(elem));
    }

    CQTKernel generateCQTkernel(const number fmax, const int bins, const number sampleFreq,
                                const number atomHopFactor=0.25, const number q=1.0,
                                const number thresh=0.0005, const bool oversampleTwo=false,
                                const bool allowSevAtoms=true, const bool perfRast=true);

    CQTResult cqt(std::vector<number>& signal, const int initSignalLen, number fmin, const number fmax,
            const int bins, const number sampleFreq, const number atomHopFactor=0.25,
            const number q=1.0, const number thresh=0.0005, const bool fromKernel=false,
            const bool oversampleTwo=false, const bool allowSevAtoms=true,
            const number coeffB=-1, const number coeffA=-1);

    CQTResult DoCQT(std::vector<number> samples, float minFreq, float maxFreq, int binCount);


#if MaxiFft

    // Compute the constant-Q transform using maximilian stuff
    CQTResult cqtMaxi(maxiRecorder recorder, maxiFFT fft, number minFreq, number maxFreq, int binCount);

    // Compute the constant-Q transform using maximilian stuff -- NOT REALTIME
    CQTResult cqtMaxi(std::vector<number>& signal, maxiFFT fft, number minFreq, number maxFreq, int binCount);

#endif

}

#endif