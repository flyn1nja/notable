#include "maximilian.h"
#include "maximilianRecording.h"
#include "libs/maxim.h"
// #include "libs/maxiFFT.h"
#include "player.h"

//* #include "FastFourierTransform.h"
//* #include "QConstTrans.h"

// #include "instruments.h"

#include <iostream>
#include <fstream>
#include <complex>
#include <tuple>

#include <sndfile.h>

#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <cstdlib>


#include "../constant-q-cpp-master/cq/ConstantQ.h"
#include "../constant-q-cpp-master/cq/CQInverse.h"
#include "../constant-q-cpp-master/cq/Chromagram.h"
#include "../constant-q-cpp-master/cq/CQSpectrogram.h"
#include "../constant-q-cpp-master/src/Pitch.h"

#include "MidiGenerator.h"


// Constant Declaration----------------------------------

// * Set this to false to have it load a file instead
constexpr bool FROM_RECORDING = true;

// * Tuning frequency (440 Hz by default)
constexpr int TUNING_FREQ = 440;

// * How much the fundamental has priority over the harmonics
// * > 1  to prioritize lower octaves (fundamental), < 1 to prioritize higher octaves (harmonics)
constexpr double LOW_OCTAVE_PRIORITY_FACTOR = 1.5; 

// * Minimal amplitude ratio required to update the max note stored
// * 0 to get any higher values, 1 to have at least double the current max
constexpr double MIN_AMP_RATIO_TO_UPDATE = 1.5; 

// * Set this to true to consider all major notes instead of only the main one
constexpr bool POLYPHONIC = false;


// Max frequency of the cqt. Usually samplerate / 3 as rule of thumb when unknown
constexpr int MAX_FREQ = 3520;

// Min frequency of the cqt
constexpr int MIN_FREQ = 110;

// Bins per octave.
constexpr int BPO = 42; // 42

// Inter buffer size
constexpr int IBS = 512;


maxiSample samplePlayback; 
maxiFFT myFFT;
// maxiIFFT myInverseFFT;
std::vector<float> mags = std::vector<float>(512);
std::vector<float> mags2 = std::vector<float>(512);
std::vector<float> phases = std::vector<float>(512);
std::vector<float> phases2 = std::vector<float>(512);
maxiEnvGen shiftEnv;

std::string audioFileName;
double maxFreq = MAX_FREQ;
double minFreq = MIN_FREQ;
int bpo = BPO;

//* CQTResult cqtRes;

// Here we define a double floating value that will contain our
// frame of lovely maximilian generated audio
double out;

// Our oscillator to fill our frame up with my favourite wave
maxiOsc osc;

// Our ramp to modulate the oscillators
maxiOsc ramp;

// We declare our recorder object here, which will call it's
// default constructor. 
maxiRecorder recorder;


// Functions Declaration----------------------------------
void setup();
void play(double *);

int runCQT(const std::string& filename, ConstantQ* constq);
int runCQTSpectrogram(const std::string& filename, CQSpectrogram* cqspect);
void processCQTFrame(std::vector<float> frame, int& frameId, CQBase* cq);
int processCQTFromFile(const std::string& filename, ConstantQ* cq);
// int processCQTSpectrFromFile(const std::string& filename, CQSpectrogram* cqspect);
int processCQTSpectrFromFile(const std::string& filename, CQSpectrogram* cqspect);

// Returns maximum amplitude of file
double findMaxFreqs(const CQSpectrogram& cq, const std::vector<CQBase::RealBlock>& blocks,
                  std::vector<double>& maxs, std::vector<double>& fmaxs, float minFreq=MIN_FREQ);
void findMaxFreqsPolyphonic(const CQSpectrogram& cq, const std::vector<CQBase::RealBlock>& blocks,
                  std::vector<std::vector<double>>& maxsPoly, std::vector<std::vector<double>>& fmaxsPoly);

void outputToFile(const std::vector<double>& maxs, const std::vector<double>& fmaxs, const float duration);
void normalizeOutput(std::vector<double>& maxs, std::vector<double>& fmaxs, double maxMax);
std::string getNoteName(float freq);



//This is main()
int main()
{
    setup();

    // TODO -- Use this one when can read from stream ok
    StartStream(std::function<void(double *)>(play));

    std::cout << "Converting samples" << std::endl;
/*
    // Remove if using double precision
    std::vector<float> samples(samplePlayback.getLength(), 0.0f);

    for (int i = 0; i < samplePlayback.getLength(); i++)
        samples[i] = static_cast<float>(samplePlayback.amplitudes[i]);
*/

    std::cout << "Running CQT" << std::endl;

    //* // For now, q-const transform is done after the recording
    //* cqtRes = QConstTrans::DoCQT(samples, 25.0f, 4500.0f, 10.0f);

    // const std::vector<std::vector<std::complex<float>>> unsparsedMat = QConstTrans::unsparse(cqtRes.spCQT);


    CQParameters params(44100,MIN_FREQ,MAX_FREQ,BPO);



    CQSpectrogram* constq = new CQSpectrogram(params, CQSpectrogram::InterpolateHold);
    processCQTSpectrFromFile(audioFileName, constq); // * Use samples instead?


    // runCQT("../Mixolydian_Mode.wav", constq); // Use samples instead?

    std::cout << "Saving to file..." << std::endl;

    // // Create and open the output file
    // //* ofstream spectFile("output.wav");
    // ofstream spectFile("spectrogram.txt");


    // // Output the spectrogram to file
    // for (int b = 0; b < binsCount; ++b)
    //     spectFile << std::to_string(b) << ": " << constq->getBinFrequency(b) << "\n";

//*
//*    // Output the spectrogram to file
//*    for (int i = 0; i < cqtRes.spCQT.width; ++i)
//*    {
//*        spectFile << i << ":\t";
//*        for (int j = 0; j < cqtRes.spCQT.height; ++j)
//*        {
//*
//*            //TODO --> DO it without unsparsing the matrix
//*            spectFile << "\t" << unsparsedMat[i][j];
//*
//*            /*currVal = cqtRes.spCQT.values;
//*
//*            spectFile << "\t";
//*
//*            if (std::get<1>(currVal) == i 
//*             && std::get<2>(currVal) == j)
//*                spectFile << abs(std::get<0>(currVal));
//*            else
//*                spectFile << "0.0";
//*            */
//*            
//*        }
//*        spectFile << "\n";
//*    }

    // constq->getBinFrequency();

    // Close the file
    // spectFile.close();

    std::cout << "Bye!" << std::endl;

    delete constq;
    return 0;
}


void setup()
{
    // samplePlayback.load("../../../beat2.wav");//load in your samples. Provide the full path to a wav file.
    audioFileName = ReadFromFile(samplePlayback, minFreq, maxFreq);
    
    myFFT.setup(1024, 512, 1024);
    // myInverseFFT.setup(1024, 512, 1024);
    shiftEnv.setup({0,50,0}, {10000,10000}, {1,1}, true);

    //------------------------------------------------

    // Call setup here, make sure you do this so the recorder
    // knows where to write the file. Currently the recorder
    // will write the wav file to the directory that this file
    // is in if you use linux but with mac and windows I 
    // strongly reccomend putting an absolute file path to the
    // directory you want to write to. Also, when in Windows,
    // remember to do double '\' characters because they
    // count as an escape which will nullify any path you write
    recorder.setup("output.wav");

    // This must be called to start the asynchronous thread to
    // manage the recorder's internal memory
    recorder.startRecording();
}

void play(double *output)
{

    // Fill our output buffer
    output[0]=out;
    output[1]=out;

    // After we have filled our output array, send the array
    // and the size of the array (in this case the amount of
    // channels, but in ofx or juce you might need to do 
    // something like channels*bufferSize).
    recorder.passData(output, maxiSettings::channels);

    // float myOut=samplePlayback.play();

    //------------------------------------------------------
    
    // if (myFFT.process(myOut)) {
    //     std::cout << "SC: " << myFFT.spectralCentroid() << endl;
        
    //     //shift some bins and phases around
    //     mags = myFFT.getMagnitudes();
    //     phases = myFFT.getMagnitudes();
    //     for(size_t b=0; b < mags.size(); b++) {
    //         size_t binIdx = b-static_cast<int>(shiftEnv.play(1));
    //         if (binIdx > mags.size()) binIdx = mags.size();
    //         if (binIdx < 0 || binIdx >= mags.size()) {
    //             mags2[b] = 0;
    //             phases2[b]=0;
    //         }else{
    //             mags2[b] = mags[binIdx];
    //             phases2[b] = phases[binIdx];
    //         }
    //     }

    //     // Todo: Do the q-const transform here for realtime!!
    //     // (Might have to use the genCQTKernel)
    // }
    // myOut = myInverseFFT.process(mags2, phases2);
    
    //output[0] is the left output. output[1] is the right output
    // output[0]=myOut;//simple as that!
    // output[1]=output[0];

    // std::cout << output[0] ;

        
    // // A pulse wave!!! Yay
    // out = osc.pulse(90, ramp.phasor(.2));
}

// We don't need to worry about telling the recorder to stop;
// when the stack unwinds and the maximillian program stops,
// the recorder will have its destructor called and the wav
// will be written for you. If you would like to do something
// more dynamic, look at the class definition in maximilian.h -
// the api allows for stricter control of the object.




// Runs the CQT from the lib
int runCQT(const std::string& filename, ConstantQ* constq)
{
    int c;
    SNDFILE *sndfile;
    // SNDFILE *sndfileOut;
    SNDFILE *sndDiffFile = 0;
    SF_INFO sfinfo;
    SF_INFO sfinfoOut;
    SF_INFO sfinfoDiff;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    sndfile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    if (!sndfile)
    {
	    std::cerr << "ERROR: Failed to open input file \"" << filename << "\": "
	         << sf_strerror(sndfile) << std::endl;
	    return 1;
    }

    sfinfoOut.channels = 1;
    sfinfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    sfinfoOut.frames = sfinfo.frames;
    sfinfoOut.samplerate = sfinfo.samplerate;
    sfinfoOut.sections = sfinfo.sections;
    sfinfoOut.seekable = sfinfo.seekable;

    // sndfileOut = sf_open(outFilename.c_str(), SFM_WRITE, &sfinfoOut) ;
    // if (!sndfileOut)
    // {
	//     std::cerr << "ERROR: Failed to open output file \"" << outFilename << "\" for writing: "
	//          << sf_strerror(sndfileOut) << std::endl;
	//     return 1;
    // }

    int ibs = 1024;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    if (maxFreq == 0.0) maxFreq = sfinfo.samplerate / 3;
    if (minFreq == 0.0) minFreq = 100;
    if (bpo == 0) bpo = 60;

    CQParameters params(sfinfo.samplerate, minFreq, maxFreq, bpo);
    ConstantQ cq(params);
    CQInverse cqi(params);

    std::cerr << "max freq = " << cq.getMaxFrequency() << ", min freq = "
	 << cq.getMinFrequency() << ", octaves = " << cq.getOctaves() << std::endl;

    std::cerr << "octave boundaries: ";
    for (int i = 0; i < cq.getOctaves(); ++i)
	    std::cerr << cq.getMaxFrequency() / pow(2, i) << " ";
    std::cerr << std::endl;

    int inframe = 0;
    int outframe = 50;
    int latency = cq.getLatency() + cqi.getLatency();

    std::vector<double> buffer;

    double maxdiff = 0.0;
    int maxdiffidx = 0;

    std::cerr << "forward latency = " << cq.getLatency() << ", inverse latency = " 
	 << cqi.getLatency() << ", total = " << latency << std::endl;

    timeval tv;
    (void)gettimeofday(&tv, 0);

    while (inframe < sfinfo.frames)
    {
        int count = -1;
	
        if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0)
            break;

        std::vector<double> cqin;
        for (int i = 0; i < count; ++i)
        {
            double v = fbuf[i * channels];
            if (channels > 1)
            {
            for (int c = 1; c < channels; ++c)
                v += fbuf[i * channels + c];

            v /= channels;
            }
            cqin.push_back(v);
        }

	// if (doDiff)
	//     buffer.insert(buffer.end(), cqin.begin(), cqin.end());

        std::vector<double> cqout = cqi.process(cq.process(cqin));

        for (int i = 0; i < int(cqout.size()); ++i)
        {
            if (cqout[i] > 1.0) cqout[i] = 1.0;
            if (cqout[i] < -1.0) cqout[i] = -1.0;
        }

        // if (outframe >= latency)
        // {
        //     sf_writef_double(sndfileOut, 
        //             cqout.data(), 
        //             cqout.size());

        // }
        // else if (outframe + (int)cqout.size() >= latency)
        // {

        //     int offset = latency - outframe;
        //     sf_writef_double(sndfileOut, 
        //             cqout.data() + offset,
        //             cqout.size() - offset);
        // }

	// if (doDiff) {
	//     for (int i = 0; i < (int)cqout.size(); ++i) {
	// 	if (outframe + i >= latency) {
	// 	    int dframe = outframe + i - latency;
	// 	    if (dframe >= (int)buffer.size()) cqout[i] = 0;
	// 	    else cqout[i] -= buffer[dframe];
	// 	    if (fabs(cqout[i]) > maxdiff &&
	// 		dframe > sfinfo.samplerate && // ignore first/last sec
	// 		dframe + sfinfo.samplerate < sfinfo.frames) {
	// 		maxdiff = fabs(cqout[i]);
	// 		maxdiffidx = dframe;
	// 	    }
	// 	}
	//     }
	    
	//     if (outframe >= latency) {

	// 	sf_writef_double(sndDiffFile, 
	// 			 cqout.data(), 
	// 			 cqout.size());

	//     } else if (outframe + (int)cqout.size() >= latency) {

	// 	int offset = latency - outframe;
	// 	sf_writef_double(sndDiffFile, 
	// 			 cqout.data() + offset,
	// 			 cqout.size() - offset);
	//     }
	// }

        inframe += count;
        outframe += cqout.size();
    }

    std::vector<double> r = cqi.process(cq.getRemainingOutput());
    std::vector<double> r2 = cqi.getRemainingOutput();

    r.insert(r.end(), r2.begin(), r2.end());

    for (int i = 0; i < int(r.size()); ++i)
    {
        if (r[i] > 1.0) r[i] = 1.0;
        if (r[i] < -1.0) r[i] = -1.0;
    }

    // sf_writef_double(sndfileOut, r.data(), r.size());

    // if (doDiff) {
	// for (int i = 0; i < (int)r.size(); ++i) {
	//     if (outframe + i >= latency) {
	// 	int dframe = outframe + i - latency;
	// 	if (dframe >= (int)buffer.size()) r[i] = 0;
	// 	else r[i] -= buffer[dframe];
	// 	if (fabs(r[i]) > maxdiff &&
	// 	    dframe > sfinfo.samplerate && // ignore first/last sec
	// 	    dframe + sfinfo.samplerate < sfinfo.frames) {
	// 	    maxdiff = fabs(r[i]);
	// 	    maxdiffidx = dframe;
	// 	}
	//     }
	// }
	// sf_writef_double(sndDiffFile, r.data(), r.size());
    // }
    outframe += r.size();

    sf_close(sndfile);
    // sf_close(sndfileOut);

    // if (doDiff) {
	//     sf_close(sndDiffFile);
    // }

    std::cerr << "in: " << inframe << ", out: " << outframe - latency << std::endl;

    // if (doDiff) {
	//     double db = 10 * log10(maxdiff);
    //     std::cerr << "max diff [excluding first and last second of audio] is "
    //         << maxdiff << " (" << db << " dBFS)"
    //         << " at sample index " << maxdiffidx << std::endl;
    // }
    
    timeval etv;
    (void)gettimeofday(&etv, 0);
        
    etv.tv_sec -= tv.tv_sec;
    if (etv.tv_usec < tv.tv_usec) {
        etv.tv_usec += 1000000;
        etv.tv_sec -= 1;
    }
    etv.tv_usec -= tv.tv_usec;
        
    double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
    std::cerr << "elapsed time (not counting init): " << sec << " sec, frames/sec at input: " << inframe/sec << std::endl;


    std::cout << std::endl << std::endl;

    for (int b = 0; b < cq.getTotalBins(); ++b)
        std::cout << std::to_string(b) << ": " << constq->getBinFrequency(b) << "\n";
    
    // // Output the spectrogram to file
    // for (int b = 0; b < cq->getTotalBins(); ++b)
    //     spectFile << std::to_string(b) << ": " << constq->getBinFrequency(b) << "\n";


    constq = std::move(&cq);
    
    return 0;
}

// Runs the CQT with spectrogram from the lib
int runCQTSpectrogram(const std::string& filename, CQSpectrogram* cqspect)
{
    int c;
    SNDFILE *sndfile;
    SNDFILE *sndDiffFile = 0;
    SF_INFO sfinfo;
    SF_INFO sfinfoOut;
    SF_INFO sfinfoDiff;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    sndfile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    if (!sndfile)
    {
	    std::cerr << "ERROR: Failed to open input file \"" << filename << "\": "
	         << sf_strerror(sndfile) << std::endl;
	    return 1;
    }

    sfinfoOut.channels = 1;
    sfinfoOut.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    sfinfoOut.frames = sfinfo.frames;
    sfinfoOut.samplerate = sfinfo.samplerate;
    sfinfoOut.sections = sfinfo.sections;
    sfinfoOut.seekable = sfinfo.seekable;

    int ibs = 1024;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    if (maxFreq == 0.0) maxFreq = sfinfo.samplerate / 3;
    if (minFreq == 0.0) minFreq = 100;
    if (bpo == 0) bpo = 60;

    CQParameters params(sfinfo.samplerate, minFreq, maxFreq, bpo);
    CQSpectrogram cqs(params, CQSpectrogram::InterpolateZeros);

    std::cerr << "max freq = " << cqs.getMaxFrequency() << ", min freq = "
	 << cqs.getMinFrequency() << ", octaves = " << cqs.getOctaves() << std::endl;

    std::cerr << "octave boundaries: ";
    for (int i = 0; i < cqs.getOctaves(); ++i)
	    std::cerr << cqs.getMaxFrequency() / pow(2, i) << " ";
    std::cerr << std::endl;

    int inframe = 0;
    int outframe = 0;
    const int latency = cqs.getLatency();

    std::vector<double> buffer;

    double maxdiff = 0.0;
    int maxdiffidx = 0;

    timeval tv;
    (void)gettimeofday(&tv, 0);

    while (inframe < sfinfo.frames)
    {
        int count = -1;
	
        if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0)
            break;

        std::vector<double> cqin;
        for (int i = 0; i < count; ++i)
        {
            double v = fbuf[i * channels];
            if (channels > 1)
            {
            for (int c = 1; c < channels; ++c)
                v += fbuf[i * channels + c];

            v /= channels;
            }
            cqin.push_back(v);
        }

        cqs.process(cqin);

        inframe += count;
    }

    cqs.getRemainingOutput();

    sf_close(sndfile);

    std::cerr << "in: " << inframe << ", out: " << outframe - latency << std::endl;
    
    timeval etv;
    (void)gettimeofday(&etv, 0);
        
    etv.tv_sec -= tv.tv_sec;
    if (etv.tv_usec < tv.tv_usec) {
        etv.tv_usec += 1000000;
        etv.tv_sec -= 1;
    }
    etv.tv_usec -= tv.tv_usec;
        
    double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
    std::cerr << "elapsed time (not counting init): " << sec << " sec, frames/sec at input: " << inframe/sec << std::endl;

    cqspect = std::move(&cqs);
    
    return 0;
}

int processCQTFromFile(const std::string& filename, ConstantQ* constq)
{
    int c;
    SNDFILE *sndfile;
    SNDFILE *sndDiffFile = 0;
    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    sndfile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    if (!sndfile)
    {
	    std::cerr << "ERROR: Failed to open input file \"" << filename << "\": "
	         << sf_strerror(sndfile) << std::endl;
	    return 1;
    }

    int ibs = 1024;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    if (maxFreq == 0.0) maxFreq = sfinfo.samplerate / 3;
    if (minFreq == 0.0) minFreq = 100;
    if (bpo == 0) bpo = 60;

    CQParameters params(sfinfo.samplerate, minFreq, maxFreq, bpo);
    ConstantQ cq(params);

    std::cerr << "max freq = " << cq.getMaxFrequency() << ", min freq = "
	 << cq.getMinFrequency() << ", octaves = " << cq.getOctaves() << std::endl;

    std::cerr << "octave boundaries: ";
    for (int i = 0; i < cq.getOctaves(); ++i)
	    std::cerr << cq.getMaxFrequency() / pow(2, i) << " ";
    std::cerr << std::endl;

    int inframe = 0;
    const int latency = cq.getLatency();

    double maxdiff = 0.0;
    int maxdiffidx = 0;

    std::cerr << "forward latency = " << cq.getLatency() << std::endl;

    timeval tv;
    (void)gettimeofday(&tv, 0);

    std::vector<CQBase::ComplexBlock> blocks;

    while (inframe < sfinfo.frames)
    {
        int count = -1;
	
        if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0)
            break;

        std::vector<double> cqin;
        for (int i = 0; i < count; ++i)
        {
            double v = fbuf[i * channels];
            if (channels > 1)
            {
            for (int c = 1; c < channels; ++c)
                v += fbuf[i * channels + c];

            v /= channels;
            }
            cqin.push_back(v);
        }

        blocks.push_back(cq.process(cqin));

        inframe += count;

        // std::cout << "Inframe #" << std::to_string(inframe) << std::endl;
    }

    sf_close(sndfile);

    
    timeval etv;
    (void)gettimeofday(&etv, 0);
        
    etv.tv_sec -= tv.tv_sec;
    if (etv.tv_usec < tv.tv_usec) {
        etv.tv_usec += 1000000;
        etv.tv_sec -= 1;
    }
    etv.tv_usec -= tv.tv_usec;
        
    double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
    std::cerr << "elapsed time (not counting init): " << sec << " sec, frames/sec at input: " << inframe/sec << std::endl;

    double max = 0.0;

    for (int i = 0; i < blocks.size(); i-=-1)
    {
        for (int m = 0; m < blocks[i].size(); ++m)
        {
            for (int n = 0; n < blocks[i][m].size(); ++n)
            {
                if (abs(blocks[i][m][n]) > max)
                    max = abs(blocks[i][m][n]);
            }

            if (max != 0.0)
            {
                std::cout << "Block #" << i << std::endl;
                std::cout << max << std::endl;
                max = 0.0;
            }
        }
        // std::cout << max << "\t\t";
        
        // std::cout << std::endl;
    }
    

    constq = std::move(&cq);
    
    return 0;
}

int processCQTSpectrFromFile(const std::string& filename, CQSpectrogram* cqspect)
{    
    SNDFILE *sndfile;
    SNDFILE *sndDiffFile = 0;
    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(SF_INFO));

    sndfile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    if (!sndfile)
    {
	    std::cerr << "ERROR: Failed to open input file \"" << filename << "\": "
	         << sf_strerror(sndfile) << std::endl;
	    return 1;
    }

    int ibs = IBS;
    int channels = sfinfo.channels;
    float *fbuf = new float[channels * ibs];

    // if (maxFreq == 0.0) maxFreq = sfinfo.samplerate / 3;
    if (maxFreq == 0.0) maxFreq = MAX_FREQ;
    if (minFreq == 0.0) minFreq = MIN_FREQ;
    // if (bpo == 0) bpo = BPO;

    CQParameters params(sfinfo.samplerate, minFreq, maxFreq, bpo);
    CQSpectrogram cq(params, CQSpectrogram::InterpolateHold);

    std::cerr << "max freq = " << cq.getMaxFrequency() << ", min freq = "
	 << cq.getMinFrequency() << ", octaves = " << cq.getOctaves() << std::endl;

    std::cerr << "octave boundaries: ";
    for (int i = 0; i < cq.getOctaves(); ++i)
	    std::cerr << cq.getMaxFrequency() / pow(2, i) << " ";
    std::cerr << std::endl;

    int inframe = 0;
    const int latency = cq.getLatency();

    double maxdiff = 0.0;
    int maxdiffidx = 0;

    std::cerr << "forward latency = " << cq.getLatency() << std::endl;

    timeval tv;
    (void)gettimeofday(&tv, 0);

    std::vector<CQBase::RealBlock> blocks;
    std::vector<double> cqin;

    int count = -1;

    while (inframe < sfinfo.frames)
    {
        count = -1;
	
        if ((count = sf_readf_float(sndfile, fbuf, ibs)) < 0)
            break;

        for (int i = 0; i < count; ++i)
        {
            double v = fbuf[i * channels];
            if (channels > 1)
            {
            for (int c = 1; c < channels; ++c)
                v += fbuf[i * channels + c];

            v /= channels;
            }
            cqin.push_back(v);
        }

        blocks.push_back(cq.process(cqin));

        cqin.clear();

        inframe += count;
    }

    timeval etv;
    (void)gettimeofday(&etv, 0);
    
    sf_close(sndfile);    
        
    etv.tv_sec -= tv.tv_sec;
    if (etv.tv_usec < tv.tv_usec) {
        etv.tv_usec += 1000000;
        etv.tv_sec -= 1;
    }
    etv.tv_usec -= tv.tv_usec;
        
    double sec = double(etv.tv_sec) + (double(etv.tv_usec) / 1000000.0);
    std::cerr << "elapsed time (not counting init): " << sec << " sec, frames/sec at input: " << inframe/sec << std::endl;

    const float duration = inframe/static_cast<float>(sfinfo.samplerate);

    if (POLYPHONIC)
    {
        std::vector<std::vector<double>> maxs = std::vector<std::vector<double>>(1, std::vector<double>());
        std::vector<std::vector<double>> fmaxs = std::vector<std::vector<double>>(1, std::vector<double>());

        // Find the maximum frequencies
        findMaxFreqsPolyphonic(cq, blocks, maxs, fmaxs);

        // // Write to file
        // outputToFile(maxs, fmaxs);

        //Generate MIDI
        MidiGenerator::outputToMIDIPolyphonic(maxs, fmaxs);
    }
    else
    {
        std::vector<double> maxs = std::vector<double>();
        std::vector<double> fmaxs = std::vector<double>();

        // Find the maximum frequencies
        double max = findMaxFreqs(cq, blocks, maxs, fmaxs);

        // Normalize the output
        normalizeOutput(maxs, fmaxs, max);

        // Write to file
        outputToFile(maxs, fmaxs, duration);

        //Generate MIDI
        MidiGenerator::outputToMIDI(maxs, fmaxs, duration);
    }

    cqspect = std::move(&cq);
    
    return 0;
}

void processCQTFrame(const std::vector<double>& frame, int& frameId, ConstantQ* cq)
{
    ++frameId;
    CQBase::ComplexBlock block = cq->process(frame);
    std::cout << std::to_string(frameId) << "e frame: " << std::endl;
    for (int i = 0; i < block.size(); ++i)
    {
        for (int j = 0; j < block[0].size(); ++j)
        {
            std::cout << std::to_string(block[i][j].real())
                << " + " << std::to_string(block[i][j].imag())<< "i \t\t";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
}

double findMaxFreqs(const CQSpectrogram& cq, const std::vector<CQBase::RealBlock>& blocks,
                  std::vector<double>& maxs, std::vector<double>& fmaxs, float minFreq)
{
    double max = 0.0;
    double maxMax = 0.0;
    double fmax = -1.0;
    double fprec = -1.0;
    double ampnext = 0.0;
    double fnext = -1.0;
    int bmax = 0;

    for (int i = 0; i < blocks.size(); i-=-1)
    {
        for (int t = 0; t < blocks[i].size(); ++t)
        {
            // Todo: Add time importance

            bmax = blocks[i][t].size()-1;

            max = blocks[i][t][bmax];
            fmax = cq.getBinFrequency(bmax);

            for (int b = blocks[i][t].size()-1; b >= 0; --b)
            // for (int b = 0; b < blocks[i][t].size(); ++b)
            {
                fnext = cq.getBinFrequency(static_cast<double>(b));
                ampnext = blocks[i][t][b];

                // The bigger the gap, the less likely it is we want a new max, especially in higher freqs
                if ((abs(ampnext * log2(b+1))
                   - abs(max * LOW_OCTAVE_PRIORITY_FACTOR * log2(bmax+1)))
                   > MIN_AMP_RATIO_TO_UPDATE * max * log2(bmax-b+1))
                // if (abs(blocks[i][t][b]) - abs(max) > max * MIN_AMP_RATIO_TO_UPDATE * log(b+1))
                {
                    // // Cancel the update if we're skipping too much coeffs and the update is not worth it
                    // if (fnext > 0.25 * (maxFreq - minFreq) + minFreq
                    // && ((fnext / fmax) > 1.75
                    // || (fnext / fprec) > 1.75)
                    // && (ampnext-max) < LOW_OCTAVE_PRIORITY_FACTOR
                    // && max != 0.0)
                    //     continue;

                    max = ampnext;
                    fmax = fnext;
                    bmax = b;

                    maxMax = std::max(max, maxMax);
                }
            }

            fprec = fmax;
            
            maxs.push_back(max);
            fmaxs.push_back(fmax);
        }
    }

    return maxMax;
}

void findMaxFreqsPolyphonic(const CQSpectrogram& cq, const std::vector<CQBase::RealBlock>& blocks,
                  std::vector<std::vector<double>>& maxsPoly, std::vector<std::vector<double>>& fmaxsPoly)
{
    double max = 0.0;
    double fmax = -1.0;
    // std::vector<double> maxRow = std::vector<double>();
    // std::vector<double> fmaxRow = std::vector<double>();

    for (int i = 0; i < blocks.size(); i-=-1)
    {
        for (int t = 0; t < blocks[i].size(); ++t)
        {
            std::vector<double> maxRow = std::vector<double>();
            std::vector<double> fmaxRow = std::vector<double>();

            for (int b = blocks[i][t].size()-1; b >= 0; --b)
            // for (int b = 0; b < blocks[i][t].size(); ++b)
            {
                if (abs(blocks[i][t][b]) - abs(max) > max * MIN_AMP_RATIO_TO_UPDATE * log(b+1))
                {
                    max = blocks[i][t][b];
                    fmax = cq.getBinFrequency(static_cast<double>(b));

                    maxRow.push_back(max);
                    fmaxRow.push_back(fmax);
                    std::cout << fmax << " -> ";
                }
            }

            std::cout << fmax << std::endl;
            max = 0.0;
            fmax = -1.0;
            // std::vector<double> maxRowCopy = maxRow;
            // std::vector<double> fmaxRowCopy = fmaxRow;
            // std::copy(maxRow.begin(), maxRow.end(), std::back_inserter(maxRowCopy));
            // std::copy(fmaxRow.begin(), fmaxRow.end(), std::back_inserter(fmaxRowCopy));
            maxsPoly.push_back(maxRow);
            fmaxsPoly.push_back(fmaxRow);
        }
    }
}


std::string getNoteName(float freq)
{
    static const char* names[] =
    {
        "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"
    };
    float cOffset;

    // float freq = cq->getBinFrequency(cq->getBinsPerOctave() - bin - 1);
    int note = Pitch::getPitchForFrequency(freq, &cOffset, TUNING_FREQ);
    // float nearestFreq = Pitch::getFrequencyForPitch(note, 0, TUNING_FREQ);
    
    // if (fabs(freq - nearestFreq) < 0.01)
    if (fabs(cOffset) < 25)
        return std::string(names[note % 12]);
    else
        if (cOffset > 0)
            return std::string(names[note % 12]) + std::string(" - ") + std::string(names[(note+1) % 12]);
        else
            return std::string(names[(note-1) % 12]) + std::string(" - ") + std::string(names[note % 12]);

    return std::string("?");
}

void outputToFile(const std::vector<double>& maxs, const std::vector<double>& fmaxs, const float duration)
{
        
    ofstream spectFile("notable_spectrogram.txt");
    float centsOffset = 0.0;
    int midiPitch = 0;
    std::string note;
    int lastN = 0;
    float lastMax = maxs[0];

    float msPerBeat = duration / maxs.size();

    for (int n = 0; n < maxs.size()-1; ++n)
        if (fabs(maxs[n] - lastMax) > 0.1)
        {
            if (maxs[n] > 0.02)
            {
                std::cout << std::endl;
                midiPitch = Pitch::getPitchForFrequency(fmaxs[n], &centsOffset);
                note = getNoteName(fmaxs[n]);

                std::printf("%.2f - %.2f |  A: %.3f - f: %.3f -- Note: %s, Midi: %d, C.Offset: %.2f",
                    lastN * msPerBeat, n * msPerBeat, maxs[n], fmaxs[n], note.c_str(), midiPitch, centsOffset);

                // ++n;
                lastN = n;
                lastMax = maxs[n];

                spectFile << "A: " << std::to_string(maxs[n]) << " - f: " << std::to_string(fmaxs[n]) << " -- Note: " << note.c_str()
                    << ",  Midi: " << std::to_string(midiPitch) << ", C.Offset: " << std::to_string(centsOffset) << std::endl;
            }
            else
            {
                std::printf("%.2f - %.2f |  A: %.3f < 0.025 - f: %.3f -- Muted",
                    lastN * msPerBeat, n * msPerBeat, maxs[n], fmaxs[n]);
                spectFile << "A: " << std::to_string(maxs[n]) << " - f: " << std::to_string(fmaxs[n]) << " -- Note: " << note.c_str()
                    << ",  Midi: " << std::to_string(midiPitch) << ", C.Offset: " << std::to_string(centsOffset) << " - MUTED -" << std::endl;
            }
        }

    // Do the last one
    int n = maxs.size()-1;

    std::cout << std::endl;

    std::printf("%.2f - %.2f |  A: %.3f - f: %.3f -- Note: %s, Midi: %d, C.Offset: %.2f",
        lastN * msPerBeat, n * msPerBeat, maxs[n], fmaxs[n], note.c_str(), midiPitch, centsOffset);

    spectFile << "A: " << std::to_string(maxs[n]) << " - f: " << std::to_string(fmaxs[n]) << " -- Note: " << note.c_str()
        << ",  Midi: " << std::to_string(midiPitch) << ", C.Offset: " << std::to_string(centsOffset) << std::endl;


    std::cout << std::endl;

    spectFile.close();
}

void normalizeOutput(std::vector<double>& maxs, std::vector<double>& fmaxs, double maxMax)
{
    // Normalize output
    for (int i = 0; i < maxs.size(); ++i)
        maxs[i] /= maxMax;
    
}