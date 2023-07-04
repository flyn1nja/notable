#include "maximilian.h"
#include "maximilianRecording.h"
#include "libs/maxim.h"
#include "player.h"

#include "QConstTrans.h"
#include "instruments.h"

#include <iostream>
#include <fstream>
#include <complex>
#include <tuple>


maxiSample samplePlayback; 
maxiFFT myFFT;
maxiIFFT myInverseFFT;
std::vector<float> mags = std::vector<float>(512);
std::vector<float> mags2 = std::vector<float>(512);
std::vector<float> phases = std::vector<float>(512);
std::vector<float> phases2 = std::vector<float>(512);
maxiEnvGen shiftEnv;

CQTResult cqtRes;

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

void setup();


// Functions Declaration----------------------------------

bool ReadFromFile(const std::string&);


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

    // For now, q-const transform is done after the recording
    cqtRes = DoCQT(samples, 25.0f, 4500.0f, 10.0f);

    std::cout << "Saving to file..." << std::endl;

    // Create and open the output file
    ofstream spectFile("spectrogram.txt");





    // std::tuple<std::complex<double>, int, int> currVal;


    const std::vector<std::vector<std::complex<float>>> unsparsedMat = unsparse(cqtRes.spCQT);

    // Output the spectrogram to file
    for (int i = 0; i < cqtRes.spCQT.width; ++i)
    {
        spectFile << i << ":\t";
        for (int j = 0; j < cqtRes.spCQT.height; ++j)
        {

            //TODO --> DO it without unsparsing the matrix
            spectFile << "\t" << unsparsedMat[i][j];

            /*currVal = cqtRes.spCQT.values;

            spectFile << "\t";

            if (std::get<1>(currVal) == i 
             && std::get<2>(currVal) == j)
                spectFile << abs(std::get<0>(currVal));
            else
                spectFile << "0.0";
            */
            
        }
        spectFile << "\n";
    }

    // Close the file
    spectFile.close();

    std::cout << "Bye!" << std::endl;
    return 0;
}


void setup() {
    samplePlayback.load("../../../beat2.wav");//load in your samples. Provide the full path to a wav file.
    myFFT.setup(1024, 512, 1024);
    myInverseFFT.setup(1024, 512, 1024);
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
    recorder.setup("lovesong.wav");

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
