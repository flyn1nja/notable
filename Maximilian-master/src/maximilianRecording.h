/*
 *  maximilian.h
 *  platform independent synthesis library using portaudio or rtaudio
 *
 *  Created by Mick Grierson on 29/12/2009.
 *  Copyright 2009 Mick Grierson & Strangeloop Limited. All rights reserved.
 *	Thanks to the Goldsmiths Creative Computing Team.
 *	Special thanks to Arturo Castro for the PortAudio implementation.
 * 
 *	Permission is hereby granted, free of charge, to any person
 *	obtaining a copy of this software and associated documentation
 *	files (the "Software"), to deal in the Software without
 *	restriction, including without limitation the rights to use,
 *	copy, modify, merge, publish, distribute, sublicense, and/or sell
 *	copies of the Software, and to permit persons to whom the
 *	Software is furnished to do so, subject to the following
 *	conditions:
 *	
 *	The above copyright notice and this permission notice shall be
 *	included in all copies or substantial portions of the Software.
 *
 *	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,	
 *	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 *	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 *	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 *	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 *	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 *	OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef MAXIMILIAN_H_REC
#define MAXIMILIAN_H_REC

//#define MAXIMILIAN_PORTAUDIO
#ifndef MAXIMILIAN_RT_AUDIO
#define MAXIMILIAN_RT_AUDIO
#endif

#include <iostream>
#include <fstream>
#include <string.h>
#include <cassert>
#include <cstdlib>
#include "math.h"
#include <cerrno>
#include <queue>
#include <vector>

#if !defined(_WIN32) && (defined(unix) || defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__)))
#define OS_IS_UNIX true
#include <pthread.h>
#include <unistd.h>
#endif

#ifdef _WIN32 //|| _WIN64
#define OS_IS_WIN true
#include <algorithm>
#include <Windows.h>
#include <process.h>
#endif

// using namespace std;
#ifndef PI
#define PI  3.1415926535897932384626433832795
#endif
#define TWOPI 6.283185307179586476925286766559

class maxiClock {
public:
    maxiClock();
    void ticker();
    void setTempo(double bpm);
    void setTicksPerBeat(int ticksPerBeat);
    // maxiOsc timer;                       // TODO: A Revoir
    int currentCount;
    int lastCount;
    int playHead;
    double bps;
    double bpm;
    int ticks;
    bool tick;
    
};

class maxiRecorder
{
public:
    maxiRecorder();
    ~maxiRecorder();

    void                setup(std::string _filename);
    void                startRecording();
    void                stopRecording();
    bool                isRecording() const;
    void                passData(double* _in, int _inBufferSize);
    void                passData(float*  _in, int _inBufferSize);
    void                saveToWav();

private:
    template <typename T>
    void                write(std::ofstream& _stream, const T& _t);
    void*               update(void* _context);
    std::vector<double> getProcessedData();
    void                enqueueBuffer();
    void                freeResources();
    bool                threadRunning;
    const int           bufferQueueSize;
    const int           bufferSize;
    long int            bufferIndex;
    long int            recordedAmountFrames;
    std::queue<double*> bufferQueue;
    std::queue<double*> savedBuffers;
    bool                doRecord;
    std::string         filename;
#if defined(OS_IS_UNIX)
	pthread_t           daemon;
	static void*        update_pthread_helper(void* _context)
#elif defined(OS_IS_WIN)
	HANDLE				daemonHandle;
	static unsigned __stdcall
                        update_pthread_helper(void* _context)
#endif
	{
		maxiRecorder* _this = static_cast<maxiRecorder*>(_context);
		_this->update(_this);
		return 0;
	}
};

#endif
