#include "maximilian.h"
#include "maximilianRecording.h"
#include "libs/maxim.h"
#include "player.h"

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
constexpr double LOW_OCTAVE_PRIORITY_FACTOR = 1.0; 

// * Minimal amplitude ratio required to update the max note stored
// * 0 to get any higher values, 1 to have at least double the current max
constexpr double MIN_AMP_RATIO_TO_UPDATE = 0.325;

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

    // Not used, but this is what is used to read samples.
    StartStream(std::function<void(double *)>(play));

    std::cout << "Running CQT" << std::endl;
    CQParameters params(44100,MIN_FREQ,MAX_FREQ,BPO);

    CQSpectrogram* constq = new CQSpectrogram(params, CQSpectrogram::InterpolateHold);
    processCQTSpectrFromFile(audioFileName, constq);

    std::cout << "Bye!" << std::endl;

    delete constq;
    return 0;
}

// Initial setup
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

// Output (not used)
void play(double *output)
{

    // Fill our output buffer
    output[0]=out;
    output[1]=out;
}

// Perform a Const Q transform on the file
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

    if (maxFreq == 0.0) maxFreq = MAX_FREQ;
    if (minFreq == 0.0) minFreq = MIN_FREQ;

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

        std::cout << "Saving to file..." << std::endl;

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