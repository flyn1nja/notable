/*
 *  player.h
 *  rtaudiotest
 *
 *  Created by Chris on 23/08/2011.
 *  Copyright 2011 Goldsmiths Creative Computing. All rights reserved.
 *
 */

//#define MAXIMILIAN_PORTAUDIO
#define MAXIMILIAN_RT_AUDIO

#include <functional>

class maxiSample;

void StartStream(std::function<void(double *)> playCBFn
    = std::function<void(double *)>());
void ReadFromFile(maxiSample& samplePlayback);
void play(double *);