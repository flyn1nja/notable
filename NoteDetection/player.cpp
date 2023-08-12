/*
 *  player.cpp
 *  rtaudiotest
 *
 *  Created by Chris on 23/08/2011.
 *  Copyright 2011 Goldsmiths Creative Computing. All rights reserved.
 *
 */

#include "player.h"
#include "maximilian.h"
#include <iostream>
#include <sys/stat.h>


#if defined( __WIN32__ ) || defined( _WIN32 )
	// #define __WINDOWS_DS__
	#include <dsound.h>
	#include <windows.h>
	#include <fstream>
	#define sleep(t) Sleep(t)
	//#include <dsconf.h>
	//#include <xaudio2.h>
#elif (defined(__APPLE__) && defined(__MACH__))
	// #define __MACOSX_CORE__
	#include "unistd.h"
	#define sleep(t) usleep(t*1000)
#else
	// #define __LINUX_ALSA__
	#include "unistd.h"
	#define sleep(t) usleep(t*1000)
#endif
#include "RtAudio.h"

 // __MACOSX_CORE__,    /*!< Macintosh OS-X Core Audio API. */
 //
 // __LINUX_ALSA__,     /*!< The Advanced Linux Sound Architecture API. */
 // __UNIX_JACK__,      /*!< The Jack Low-Latency Audio Server API. */
 // __LINUX_PULSE__,    /*!< The Linux PulseAudio API. */
 // __LINUX_OSS__,      /*!< The Linux Open Sound System API. */
 //
 // __WINDOWS_ASIO__,   /*!< The Steinberg Audio Stream I/O API. */
 // __WINDOWS_WASAPI__, /*!< The Microsoft WASAPI API. */
 // __WINDOWS_DS__,     /*!< The Microsoft DirectSound API. */

std::function<void(double *)> playCBFunc = run;


inline bool fileExists (const std::string& name);

// maxiSample samplePlayback;

// void setup()
// {	
// 	// std::string filename;
// 	// bool badFile = false;

// 	// const std::string defaultFile("/home/alexis/Bureau/Notable/notable/NoteDetection/snare.wav");
	
// 	// do 
// 	// {
// 	// 	std::cout << "Enter file to use (default to snare.wav): \t";
// 	// 	std::cin >> filename;
		
// 	// 	badFile = !ReadFromFile(filename.length() > 1 ? filename : defaultFile);
		
// 	// 	if (badFile)
// 	// 		std::cout << "File cannot be found or opened. Try another file." << std::endl;

// 	// 	Sleep(1000);
// 	// } while (badFile);
// }

void run(double *output) //run dac! Very very often. Too often in fact. er...
{
// 	output[0] = samplePlayback.playAtSpeed(0.68);
    
// //  output[1] = output[0];
}

int routing	(void *outputBuffer, void *inputBuffer, unsigned int nBufferFrames,
			 double streamTime, RtAudioStreamStatus status, void *userData ) {
	
	
	double *buffer = (double *) outputBuffer;
	double *lastValues = (double *) userData;
	//	double currentTime = (double) streamTime; Might come in handy for control
	if ( status )
		std::cout << "Stream underflow detected!" << std::endl;
	for (size_t i=0; i<nBufferFrames; i++ ) {	
	}
	// Write interleaved audio data.
	for (size_t i=0; i<nBufferFrames; i++ ) {
		playCBFunc(lastValues);			
		for (size_t j=0; j<maxiSettings::channels; j++ ) {
			*buffer++=lastValues[j];
		}
	}
	return 0;
}

void errorCallback( RtAudioErrorType /*type*/, const std::string &errorText )
{
  std::cerr << "\nRtAudio errorCallback: " << errorText << "\n\n";
}

void StartStream(std::function<void(double *)> playCBFn)
{
	RtAudio dac(RtAudio::UNSPECIFIED, &errorCallback);
	if ( dac.getDeviceCount() < 1 ) {
		std::cout << "\nNo audio devices found!\n";
		char input;
		std::cin.get( input );
		exit( 0 );
	}

	RtAudio::StreamParameters parameters;
	parameters.deviceId = dac.getDefaultOutputDevice();
	parameters.nChannels = maxiSettings::channels;
	parameters.firstChannel = 0;
	unsigned int sampleRate = maxiSettings::sampleRate;
	unsigned int bufferFrames = maxiSettings::bufferSize; 
	std::vector<double> data(maxiSettings::channels,0);
	
	if (playCBFn) playCBFunc = playCBFn;

	dac.openStream( &parameters, NULL, RTAUDIO_FLOAT64,
						sampleRate, &bufferFrames, &routing, (void *)&(data[0]));
	
	dac.startStream();
	if ( dac.isStreamOpen()){

		
		char input;
		std::cout << "\nMaximilian is playing ... press <enter> to quit.\n";
		std::cin.get( input );
		
		dac.stopStream();
	}
	if ( dac.isStreamOpen() ) dac.closeStream();
}

// Returns the filename
const std::string ReadFromFile(maxiSample& samplePlayback, double& minFreq, double& maxFreq)
{
	std::string fileNameStr = "";
	std::string minFreqStr = "";
	std::string maxFreqStr = "";

	bool badFile = false;

	const std::string defaultFile("../Mixolydian_Mode.wav");
	
	do 
	{
		std::cout << "Enter file to use (enter \".\" to use default file Mixolydian_Mode.wav): \t";
		std::cin >> fileNameStr;

		if (fileNameStr.length() <= 1)
			fileNameStr = defaultFile;

		std::cout << "File is " << fileNameStr << std::endl;
		
		// badFile = !samplePlayback.load(filename.length() > 1 ? filename : defaultFile, 0);
		badFile = !fileExists(fileNameStr);
		
		if (badFile)
		{
			std::cout << "File cannot be found or opened. Try another file." << std::endl;
			sleep(1000);
		}

	} while (badFile);

	try
	{
		std::cout << "Min freq (enter \".\" to use default 110 Hz): \t";
		std::cin >> minFreqStr;

		if (minFreqStr == ".") minFreq = 0;

		minFreq = std::stod(minFreqStr);
	}
	catch(const std::exception& e)
	{
		minFreq = 0;
	}

	try
	{
		std::cout << "Max freq (enter \".\" to use default 3520 Hz): \t";
		std::cin >> maxFreqStr;

		if (maxFreqStr == ".") maxFreq = 0;

		maxFreq = std::stod(maxFreqStr);
	}
	catch(const std::exception& e)
	{
		maxFreq = 0;
	}

	return fileNameStr;
}


// #if defined(__linux__) || defined(__APPLE__)
inline bool fileExists (const std::string& name)
{
    struct stat buffer;   
    return stat(name.c_str(), &buffer) == 0; 
}
// #elif defined(__MINGW32__) || defined(__WIN32__)
// inline bool fileExists (const std::string& name)
// {
//     if (FILE *file = fopen(name.c_str(), "r"))
// 	{
//         fclose(file);
//         return true;
//     }
// 	else
// 	{
//         return false;
//     }   
// }
// #endif