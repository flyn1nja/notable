// Q Constant Transform
// By Alexis Giguère-Lebel 05-31-2023
// Based on the lib taken from https://users.iem.at/schoerkhuber/cqt2010/
//// Based on the lib taken from https://web.medis.mit.edu/~brown/cqtrans.htm

#include "QConstTrans.h"
#include "FastFourierTransform.h"

CQTKernel QConstTrans::generateCQTkernel(const number fmax, const int bins, const number sampleFreq,
                            const number atomHopFactor, const number q,
                            const number thresh, const bool oversampleTwo,
                            const bool allowSevAtoms, const bool perfRast)
{
    // define
    number fmin = (fmax/2.0)*pow(2.0,(1.0/static_cast<number>(bins)));
    number Q = (1.0/(pow(2.0,(1.0/static_cast<number>(bins))-1)))*q;
    number NkMax = Q * sampleFreq / fmin;
    NkMax = 2 * round((NkMax+1)/2.0) - 1; //modhann windows are always of odd length

    ComplexVector specKernel; // TODO -- Revoir pour le type de datastructure

    // Compute FFT size, FFT hop, atom hop, ...
    number NkMin = Q * sampleFreq / (fmin*pow(2.0,((bins-1)/bins))); //length of the shortest atom [samples]
    int atomHOP = (2*round(NkMin*atomHopFactor/2.0f)); //force atomHOP to be even so that for the shifted kernel (oversampling 2->atomHop/2.0f) no fractional delay occurs

    int firstCenter = ceil(NkMax/2.0f); //first possible center position within the frame
    firstCenter = atomHOP * ceil(firstCenter/atomHOP); //lock the first center to an integer multiple of the atom hop size
    number FFTLen = pow(2.0,nextpow2(firstCenter+ceil(NkMax/2.0f))); //use smallest possible FFT size (increase sparsity)

    number winNr;
    if (perfRast)
    {
        winNr = floor((FFTLen-ceil(NkMax/2.0f)-firstCenter)/atomHOP); //number of temporal atoms per FFT Frame
        if (winNr == 0)
        {
            FFTLen *= 2;
            winNr = floor((FFTLen-ceil(NkMax/2.0f)-firstCenter)/atomHOP);
        }
    }
    else if (oversampleTwo)
    {
        winNr = floor((FFTLen-firstCenter+1-atomHOP/2.0f-(NkMax-1)/2.0f)/atomHOP)+1; //number of temporal atoms per FFT Frame
        if (winNr == 0)
        {
            FFTLen = FFTLen * 2;
            winNr = floor((FFTLen-firstCenter+1-atomHOP/2.0f-(NkMax-1)/2.0f)/atomHOP)+1;
        }
    }
    else
        winNr = floor((FFTLen-ceil(NkMax/2.0f)-firstCenter)/atomHOP)+1; //number of temporal atoms per FFT Frame

    if (allowSevAtoms == 0)
        winNr = 1;

    number last_center = firstCenter + (winNr-1)*atomHOP;
    number fftHOP = (last_center + atomHOP) - firstCenter; // hop size of FFT frames
    number fftOLP = (FFTLen-fftHOP/FFTLen)*100; //overlap of FFT frames in percent ***AK:needed?

    // init variables
    std::complex<number> tempKernel[static_cast<int>(FFTLen)] = { }; // Vertical, zero-filled
    // ComplexVector sparseKernel = ComplexVector(); 
    sparseMat<std::complex<number>> sparseKernel = sparseMat<std::complex<number>>(); 

    // Compute kernel
    int atomInd = 0;
    number Nk;
    for (int k=1; k < bins; ++k)
    {
        number fk = fmin*pow(2.0,((k-1)/static_cast<number>(bins)));
        Nk = Q*sampleFreq/fk;

        std::vector<number> winFct = modhann(Nk);
        std::vector<number> n(winFct.size(), 0.0);

        for (int i = 0; i < winFct.size(); ++i)
            n[i] = (winFct.size()-1) / 2.0 + i;
        
        for (int j = 0; j < winFct.size(); ++j)
            winFct[j] /= sumVect(winFct);


        using namespace std::complex_literals;

        ComplexVector tempKernelBin;
        for (int l = 0; l < winFct.size(); ++l)
            tempKernelBin.push_back(exp(static_cast<float>(2.0 * PI) * n[l] * fk/sampleFreq * -1if) * winFct[l]);
        

        // std::complex<number> tempKernelBin = (/*winFct/sumVect(winFct)*/winFct).* exp(-std::complex_literals::i*2*pi*n*fk/sampleFreq);
        const number atomOffset = firstCenter-ceil(Nk/2.0f);

        number shift;

        for (int m = 0; m < winNr; ++m)
        {
            shift = atomOffset + (m * atomHOP);
            // tempKernel(1+shift:length(tempKernelBin)+shift) = tempKernelBin;

            for (int n = shift; n < tempKernelBin.size() + shift - 1 /*TODO: revoir pour -1*/; ++n)
                tempKernel[n] = tempKernelBin[n];

            ++atomInd;

            // * specKernel = fft(conjugate(tempKernel, static_cast<int>(FFTLen))); // Do fft
            specKernel = FastFourierTransform::DoFFT(conjugate(tempKernel, static_cast<int>(FFTLen)),true); // Do fft

                for (int ker = 0; ker < specKernel.size(); ++ker)
                {
                    if (abs(specKernel[ker]) <= thresh)
                        specKernel[ker] = 0.0;
                }
            
            // specKernel(abs(specKernel)<=thresh) = 0;

            // Sparse removes 0 and stores only the data to save space.
            // sparseKernel = sparse([sparseKernel; specKernel]);
            sparseAppend(sparseKernel, specKernel, thresh);

            // tempKernel = zeros(1,FFTLen); //reset window
            std::complex<number> tempKernel[static_cast<int>(FFTLen)] = { };
        }

    }

    // sparseKernel = (sparseKernel.')/FFTLen;
    
    // // TRANSPOSER LE sparseKernel!!
    
    transposeMat(sparseKernel); // ! Revoir pour les get<1>() et get<2>() si transposée

    // Normalize the magnitudes of the atoms
    // [ignore,wx1] = max(sparseKernel(:,1));
    // [ignore,wx2] = max(sparseKernel(:,end));
    // wK=sparseKernel(wx1:wx2,:);
    // wK = diag(wK * wK');
    // wK = wK(round(1/q)+1:(end-round(1/q)-2));

    ComplexVector wk = getWK(sparseKernel, q);
    // weight = 1./mean(abs(wk));
    // weight = weight.*(fftHOP/FFTLen); 
    // weight = sqrt(weight); //sqrt because the same weight is applied in icqt again

    // return
    // cqtKernel = struct('fKernel',sparseKernel,'fftLEN',FFTLen,'fftHOP',fftHOP,'fftOverlap',fftOLP,'perfRast',perfRast,...
    //     'bins',bins,'firstcenter',firstCenter,'atomHOP',atomHOP,'atomNr',winNr,'NkMax',NkMax,'Q',Q,'fmin',fmin);
    
    CQTKernel qcKernel;

    qcKernel.freqKernel = sparseKernel;

    qcKernel.fftLen = FFTLen;
    qcKernel.fftHop = fftHOP;

    qcKernel.atomHop = atomHOP;
    qcKernel.atomNr = winNr;
    qcKernel.firstcenter = firstCenter;

    qcKernel.nkMax = NkMax;
    qcKernel.nkMin = NkMin;

    normalizeByWeight(qcKernel, wk);
    
    return qcKernel;
}

// void cellToSparse() 
// {   
// #if true
// // this version has big memory consumption but is very fast
//     emptyHops = firstcenter/atomHOP;
//     drop = emptyHops*2^(octaves-1)-emptyHops; //distance between first value in highest octave and first value in lowest octave
//     spCQT = zeros(bins*octaves,size(Xcq{1},2)*atomNr-drop);

//     for i=1:octaves
//         drop = emptyHops*2^(octaves-i)-emptyHops; %first coefficients of all octaves have to be in synchrony
//         X = Xcq{i}; 
//         if  atomNr > 1 %more than one atom per bin --> reshape
//             Xoct = zeros(bins,atomNr*size(X,2)-drop);
//             for u=1:bins %reshape to continous windows for each bin (for the case of several wins per frame)
//                octX_bin = X((u-1)*atomNr+1:u*atomNr,:);
//                Xcont = reshape(octX_bin,1,size(octX_bin,1)*size(octX_bin,2));
//                Xoct(u,:) = Xcont(1+drop:end);
//             end
//             X = Xoct;
//         else
//             X = X(:,1+drop:end);
//         end
//         binVec = bins*octaves-bins*i+1:bins*octaves-bins*(i-1);
//         spCQT(binVec,1:2^(i-1):size(X,2)*2^(i-1)) = X;

//     end
//     spCQT = sparse(spCQT); %storing as sparse matrix at the end is the fastest way. Big memory consumption though!

// #else
// // this version uses less memory but is noticable slower
//     emptyHops = firstcenter/atomHOP;
//     drops = emptyHops*2.^(octaves-(1:octaves))-emptyHops;
//     len = max(((atomNr*cellfun('size',Xcq,2)-drops).*2.^(0:octaves-1))); %number of columns of output matrix
//     spCQT = [];

//     for i=octaves:-1:1
//         drop = emptyHops*2^(octaves-i)-emptyHops; %first coefficients of all octaves have to be in synchrony
//         X = Xcq{i}; 
//         if  atomNr > 1 %more than one atom per bin --> reshape
//             Xoct = zeros(bins,atomNr*size(X,2)-drop);
//             for u=1:bins %reshape to continous windows for each bin (for the case of several wins per frame)
//                octX_bin = X((u-1)*atomNr+1:u*atomNr,:);
//                Xcont = reshape(octX_bin,1,size(octX_bin,1)*size(octX_bin,2));
//                Xoct(u,:) = Xcont(1+drop:end);
//             end
//             X = Xoct;
//         else
//             X = X(:,1+drop:end);
//         end
//         X = upsample(X.',2^(i-1)).';
//         X = [X zeros(bins,len-size(X,2))];
//         spCQT = sparse([spCQT; X]);  
//     end

// #endif
// }

CQTResult QConstTrans::cqt(std::vector<number>& signal, const int initSignalLen, number fmin, const number fmax,
         const int bins, const number sampleFreq, const number atomHopFactor,
         const number q, const number thresh, const bool fromKernel,
         const bool oversampleTwo, const bool allowSevAtoms,
         const number coeffB, const number coeffA)
{
    const number OCTAVE_NR = ceil(log2(fmax/fmin)); // int?
    fmin = (fmax/(pow(2.0,OCTAVE_NR))) * pow(2.0, 1.0/static_cast<number>(bins));

    // const number q = 1.0;
    // const number atomHopFactor = 0.25;
    // const number thresh = 0.0005;
    // const bool allowSevAtoms = true;
    // const bool oversampTwo = false;
    const bool perfRast = true;


    // Do once
    static bool kernelInited = false;
    static CQTKernel cqtKernel;

    // Design the kernel for one octave
    if (!fromKernel && !kernelInited)
        cqtKernel = generateCQTkernel(fmax, bins, sampleFreq,
                         atomHopFactor, q, thresh, oversampleTwo,
                         allowSevAtoms, perfRast);

    kernelInited = true;
    
    // Calculate the CQT
    const int maxBlock = cqtKernel.fftLen * std::pow(2, OCTAVE_NR - 1);

    std::vector<number> paddedSignal(maxBlock + initSignalLen + maxBlock, 0.0);
    std::copy(signal.begin(), signal.end(), paddedSignal.begin() + maxBlock);

    const int OVRLP = cqtKernel.fftLen - cqtKernel.fftHop;
    const number totBinsInOctaves = bins * OCTAVE_NR; // int?

    int emptyHops = cqtKernel.firstcenter / cqtKernel.atomHop;
    int fftBlockNr = std::ceil((paddedSignal.size() - cqtKernel.fftLen) / (cqtKernel.fftLen - OVRLP));
    // Todo: rendre sparse
    ComplexMatrix spCQT(totBinsInOctaves, ComplexVector(fftBlockNr * cqtKernel.atomNr, 0.0));
   
    int drops = 0;
    std::pair<int,int> binVecBorders(0,0);
    // std::vector<number> bufferedSignal(OVRLP,0.0);
    ComplexVector signalFft(cqtKernel.fftLen,0.0);
    std::copy(signal.begin(), signal.end(), paddedSignal.begin() + maxBlock);
    
    // ComplexMatrix signalOctTrans(cqtKernel.fftLen,ComplexVector(std::complex<number>(0.0)));

    std::vector<int> tVec(cqtKernel.fftLen,0.0);

    for (int i = 0; i < OCTAVE_NR; ++i)
    {
        binVecBorders.first = bins * (OCTAVE_NR - i);
        binVecBorders.second = bins * (OCTAVE_NR - i+1)-1;
        drops = emptyHops * pow(2, OCTAVE_NR-i-1) - emptyHops; // synchronize first coeffs of all octaves
        
        // ! À Revoir. Utilise la fonction buffer de matlab (pour realtime)...
        signalFft = FastFourierTransform::DoFFT(paddedSignal/* -from-  binVecBorders.first   -to-   binVecBorders.second*/);

        //TODO : Matrix multiplication
        // ComplexMatrix signalOct = cqtKernel.freqKernel*signalFft; //TODO : Matrix multiplication

        sparseMat<std::complex<number>> bufferedFFT = bufferVec(signalFft,cqtKernel.fftHop,1,thresh);

        ComplexMatrix signalOct = unsparse<std::complex<number>>(multiplyMatrices(cqtKernel.freqKernel, bufferedFFT, true));

        //Reshape?

        int k = 0;
        int kInd = 0;
        while (k < signalOct[0].size() * pow(2,i))
        {
            tVec[kInd] = k;
            k += pow(2.0,i);
            ++kInd;
        }
        
        for (int b = binVecBorders.first; b < binVecBorders.second; ++b)
        for (int t = 0; t < tVec.size(); ++t)
            spCQT[b][tVec[t]] = signalOct[b][t]; // Todo: Revoir pour indexation du 2e
            
        
         


        tVec.clear();
        // for (int b = 0; b < drops; ++b)
        // {
        //     for (int j = 0; j < m_p.atomsPerFrame; ++j)
        //     {
        //         int target = base + (b * (totalColumns / blocksThisOctave) +
        //                     (j * ((totalColumns / blocksThisOctave) /
        //                         m_p.atomsPerFrame)));
        //         while (int(out[target].size()) <
        //                 m_p.binsPerOctave * (octave + 1)) {
        //             out[target].push_back(Complex());
        //         }                    
        //         for (int i = 0; i < m_p.binsPerOctave; ++i) {
        //             out[target][m_p.binsPerOctave * octave + i] =
        //                 block[j][m_p.binsPerOctave - i - 1];
        //         }
        //     }
        // }
        
        
        //
        //
        //

/*
        bool enough = true;
        while (enough) {
            // We could have quite different remaining sample counts in
            // different octaves, because (apart from the predictable
            // added counts for decimator output on each block) we also
            // have variable additional latency per octave
            enough = true;
            for (int i = 0; i < m_octaves; ++i) {
                int required = m_p.fftSize * pow(2, m_octaves - i - 1);
                if ((int)m_buffers[i].size() < required) {
                    enough = false;
                }
            }
            if (!enough) break;
            int base = out.size();
            int totalColumns = pow(2, m_octaves - 1) * m_p.atomsPerFrame;
            for (int i = 0; i < totalColumns; ++i) {
                out.push_back(ComplexColumn());
            }
            for (int octave = 0; octave < m_octaves; ++octave) {
                int blocksThisOctave = pow(2, (m_octaves - octave - 1));
                for (int b = 0; b < blocksThisOctave; ++b) {
                    ComplexBlock block = processOctaveBlock(octave);                
                    for (int j = 0; j < m_p.atomsPerFrame; ++j) {
                        int target = base +
                                (b * (totalColumns / blocksThisOctave) +
                                    (j * ((totalColumns / blocksThisOctave) /
                                        m_p.atomsPerFrame)));
                        while (int(out[target].size()) <
                                m_p.binsPerOctave * (octave + 1)) {
                            out[target].push_back(Complex());
                        }                    
                        for (int i = 0; i < m_p.binsPerOctave; ++i) {
                            out[target][m_p.binsPerOctave * octave + i] =
                                block[j][m_p.binsPerOctave - i - 1];
                        }
                    }
                }
            }
        }
*/

        // Anti aliasing filter
        if (i != OCTAVE_NR)
        {
            antialias(signal);

            // Re-sample to interval of 2
            for (int s = 0; s < signal.size(); s+=2)
                signal[s/2] = signal[s];
            
        }

    } 


    // Create CQTResult object and fill in the values
    CQTResult result;
    // result.spCQT = spCQT;
    result.spCQT = sparse(spCQT);
    result.fKernel = cqtKernel;
    result.fmin = fmin;
    result.fmax = fmax;
    result.octaveNr = OCTAVE_NR;
    result.bins = bins;
   
    return result;

}

// Not realtime
CQTResult QConstTrans::DoCQT(std::vector<number> samples, float minFreq, float maxFreq, int binCount)
{
    return cqt(samples, samples.size(), minFreq, maxFreq, binCount, 1.0f);
}


#if MaxiFft

// Compute the constant-Q transform using maximilian stuff
CQTResult QConstTrans::cqtMaxi(maxiRecorder recorder, maxiFFT fft, number minFreq, number maxFreq, int binCount)
{

}

// Compute the constant-Q transform using maximilian stuff -- NOT REALTIME
CQTResult QConstTrans::cqtMaxi(std::vector<number>& signal, maxiFFT fft, number minFreq, number maxFreq, int binCount)
{

    // !! REVOIR pour params
    return cqt(signal, signal.size(), minFreq, maxFreq, binCount,1);
}

#endif

/*
// #include <iostream>
// #include <vector>
// #include <cmath>

CQTResult cqt2(const std::vector<number>& x, number fmin, number fmax, int bins, number fs,
              number q = 1.0, number atomHopFactor = 0.25, number thresh = 0.0005,
              int allowSevAtoms = 1, int oversampTwo = 0) {
   
    // Input checking
    if (x.size() > 1 && x[0].size() > 1) {
        throw std::runtime_error("cqt requires one-dimensional input!");
    }
    if (x.size() > 1) {
        // Convert input signal to column vector
        std::vector<number> x_col(x.size());
        for (int i = 0; i < x.size(); i++) {
            x_col[i] = x[i][0];
        }
        x = x_col;
    }
   
    // Default parameter values
    // Modify these as needed
    number q = 1.0;
    number atomHopFactor = 0.25;
    number thresh = 0.0005;
    int allowSevAtoms = 1;
    int oversampTwo = 0;
    int perfRast = 0;
   
    // Design lowpass filter
    // Modify filter parameters as needed
    int LPorder = 6;
    number cutoff = 0.5;
    std::vector<number> B(LPorder + 1), A(LPorder + 1);
    butter(LPorder, cutoff, B, A);
   
    // Define octave number
    int octaveNr = std::ceil(std::log2(fmax / fmin));
    fmin = (fmax / std::pow(2, octaveNr)) * std::pow(2, 1.0 / bins); // Update fmin to actual value
   
    // Design kernel for one octave
    // Generate the kernel using the genCQTkernel function (not shown here)
    ComplexMatrix cqtKernel = genCQTkernel(fmax, bins, fs,
        "q", q, "atomHopFactor", atomHopFactor, "thresh", thresh,
        "allowSevAtoms", allowSevAtoms, "perfRast", perfRast, "oversampTwo", oversampTwo);
   
    // Calculate CQT
    int maxBlock = cqtKernel[0].size() * std::pow(2, octaveNr - 1);
    int suffixZeros = maxBlock;
    int prefixZeros = maxBlock;
    std::vector<number> x_padded(prefixZeros + x.size() + suffixZeros, 0.0);
    std::copy(x.begin(), x.end(), x_padded.begin() + prefixZeros);
   
    int OVRLP = cqtKernel[0].size() - cqtKernel[1].size();
    int emptyHops = cqtKernel[2][0].number() / cqtKernel[5][0].number();
    int fftBlockNr = std::ceil((x_padded.size() - cqtKernel[0].size()) / (cqtKernel[0].size() - OVRLP));
   
    Matrix spCQT(bins, std::vector<number>(fftBlockNr * (x_padded.size() / cqtKernel[0].size()), 0.0));
   
    for (int i = 0; i < fftBlockNr; i++) {
        int segmentOffset = (cqtKernel[0].size() - OVRLP) * i;
       
        std::vector<number> segment(x_padded.begin() + segmentOffset, x_padded.begin() + segmentOffset + cqtKernel[0].size());
       
        std::vector<number> segment_hann(segment.size(), 0.0);
        for (int j = 0; j < segment.size(); j++) {
            segment_hann[j] = segment[j] * (0.5 - 0.5 * std::cos(2.0 * M_PI * j / (segment.size() - 1)));
        }
       
        std::vector<number> segment_filtered(segment.size(), 0.0);
        filter(segment_hann, B, A, segment_filtered);
       
        Matrix D(cqtKernel.size());
        for (int j = 0; j < cqtKernel.size(); j++) {
            D[j].resize(segment.size(), 0.0);
        }
       
        for (int j = 0; j < segment.size(); j++) {
            for (int k = 0; k < cqtKernel.size(); k++) {
                D[k][j] = segment_filtered[j] * cqtKernel[k][j].number();
            }
        }
       
        ComplexMatrix X(cqtKernel.size());
        for (int j = 0; j < cqtKernel.size(); j++) {
            X[j].resize(segment.size(), std::complex<number>(0.0, 0.0));
        }
       
        fft(D, X);
       
        Matrix Z(X.size());
        for (int j = 0; j < X.size(); j++) {
            Z[j].resize(X[j].size(), 0.0);
        }
       
        for (int j = 0; j < X.size(); j++) {
            for (int k = 0; k < X[j].size(); k++) {
                Z[j][k] = std::pow(std::abs(X[j][k]), 2.0);
            }
        }
       
        for (int j = 0; j < bins; j++) {
            int blockOffset = (cqtKernel[0].size() - cqtKernel[1].size()) * j;
            std::vector<number> Z_bin(Z[0].begin() + blockOffset, Z[0].begin() + blockOffset + Z[1].size());
            std::copy(Z_bin.begin(), Z_bin.end(), spCQT[j].begin() + i * Z_bin.size());
        }
    }
   
    // Create CQTResult object and fill in the values
    CQTResult result;
    result.spCQT = spCQT;
    result.fKernel = cqtKernel;
    result.fmin = fmin;
    result.fmax = fmax;
    result.octaveNr = octaveNr;
    result.bins = bins;
   
    return result;
}
*/

/*
int main() {
    // Example usage
    std::vector<number> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    number fmin = 10.0;
    number fmax = 100.0;
    int bins = 10;
    number fs = 44100.0;
   
    CQTResult result = cqt(x, fmin, fmax, bins, fs);
   
    // Print the spectrogram
    for (int i = 0; i < result.spCQT.size(); i++) {
        for (int j = 0; j < result.spCQT[i].size(); j++) {
            std::cout << result.spCQT[i][j] << " ";
        }
        std::cout << std::endl;
    }
   
    return 0;
}
*/