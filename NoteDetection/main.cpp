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



#include "../constant-q-cpp-master/cq/ConstantQ.h"
#include "../constant-q-cpp-master/cq/CQInverse.h"

#include <sndfile.h>

#include <iostream>

#include <getopt.h>
#include <unistd.h>
#include <sys/time.h>
#include <cstdlib>



maxiSample samplePlayback; 
maxiFFT myFFT;
// maxiIFFT myInverseFFT;
std::vector<float> mags = std::vector<float>(512);
std::vector<float> mags2 = std::vector<float>(512);
std::vector<float> phases = std::vector<float>(512);
std::vector<float> phases2 = std::vector<float>(512);
maxiEnvGen shiftEnv;

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


// * Set thisto falase to have it load a file instead
constexpr bool FROM_RECORDING = true;

// bool ReadFromFile(const std::string&);
// bool StartStream(std::function<void(double *)> playCBFn);

//This is main()
int main()
{
    setup();

    // TODO -- Use this one when can read from stream ok
    StartStream(play);

    
    std::cout << "Converting samples" << std::endl;

    // Remove if using double precision
    std::vector<float> samples(samplePlayback.getLength(), 0.0f);

    for (int i = 0; i < samplePlayback.getLength(); i++)
        samples[i] = static_cast<float>(samplePlayback.amplitudes[i]);
    

    std::cout << "Running CQT" << std::endl;

    //* // For now, q-const transform is done after the recording
    //* cqtRes = QConstTrans::DoCQT(samples, 25.0f, 4500.0f, 10.0f);

    std::cout << "Saving to file..." << std::endl;

    // Create and open the output file
    ofstream spectFile("spectrogram.txt");





    // std::tuple<std::complex<double>, int, int> currVal;


    //* const std::vector<std::vector<std::complex<float>>> unsparsedMat = QConstTrans::unsparse(cqtRes.spCQT);

    CQParameters params(0,0,0,0);
    ConstantQ* constq = new ConstantQ(params);


    runCQT("output.wav", constq);
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

    // Close the file
    spectFile.close();

    std::cout << "Bye!" << std::endl;
    return 0;
}


void setup() {
    samplePlayback.load("../../../beat2.wav");//load in your samples. Provide the full path to a wav file.
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

void play(double *output) {

    // Fill our output buffer
    output[0]=out;
    output[1]=out;

    // After we have filled our output array, send the array
    // and the size of the array (in this case the amount of
    // channels, but in ofx or juce you might need to do 
    // something like channels*bufferSize).
    recorder.passData(output, maxiSettings::channels);

    float myOut=samplePlayback.play();

    //------------------------------------------------------
    
    if (myFFT.process(myOut)) {
        std::cout << "SC: " << myFFT.spectralCentroid() << endl;
        
        //shift some bins and phases around
        mags = myFFT.getMagnitudes();
        phases = myFFT.getMagnitudes();
        for(size_t b=0; b < mags.size(); b++) {
            size_t binIdx = b-static_cast<int>(shiftEnv.play(1));
            if (binIdx > mags.size()) binIdx = mags.size();
            if (binIdx < 0 || binIdx >= mags.size()) {
                mags2[b] = 0;
                phases2[b]=0;
            }else{
                mags2[b] = mags[binIdx];
                phases2[b] = phases[binIdx];
            }
        }

        // Todo: Do the q-const transform here for realtime!!
        // (Might have to use the genCQTKernel)
    }
    // myOut = myInverseFFT.process(mags2, phases2);
    
    //output[0] is the left output. output[1] is the right output
    output[0]=myOut;//simple as that!
    output[1]=output[0];

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
    double maxFreq = 0;
    double minFreq = 0;
    int bpo = 0;
    
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
    int outframe = 0;
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

    delete fbuf;

    constq = std::move(&cq);
    
    return 0;
}