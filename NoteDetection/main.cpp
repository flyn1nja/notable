
#include "maximilian.h"
#include "maximilianRecording.h"
#include "libs/maxim.h"


/*
maxiSample beats; 
maxiFFT myFFT;
maxiIFFT myInverseFFT;
std::vector<float> mags = std::vector<float>(512);
std::vector<float> mags2 = std::vector<float>(512);
std::vector<float> phases = std::vector<float>(512);
std::vector<float> phases2 = std::vector<float>(512);
maxiEnvGen shiftEnv;


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
*/
void setup() {
/*    
    beats.load("../../../beat2.wav");//load in your samples. Provide the full path to a wav file.
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
    recorder.startRecording();*/
}

void play(double *output) {
/*    
    
    float myOut=beats.play();
    
    if (myFFT.process(myOut)) {
        cout << "SC: " << myFFT.spectralCentroid() << endl;
        
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
    }
    myOut = myInverseFFT.process(mags2, phases2);
    
    //output[0] is the left output. output[1] is the right output
    output[0]=myOut;//simple as that!
    output[1]=output[0];

    //------------------------------------------------------

        
    // A pulse wave!!! Yay
    out = osc.pulse(90, ramp.phasor(.2));
    
    // Fill our output buffer
    output[0]=out;
    output[1]=out;

    // After we have filled our output array, send the array
    // and the size of the array (in this case the amount of
    // channels, but in ofx or juce you might need to do 
    // something like channels*bufferSize).
    recorder.passData(output, maxiSettings::channels);*/
    
}

// We don't need to worry about telling the recorder to stop;
// when the stack unwinds and the maximillian program stops,
// the recorder will have its destructor called and the wav
// will be written for you. If you would like to do something
// more dynamic, look at the class definition in maximilian.h -
// the api allows for stricter control of the object.
