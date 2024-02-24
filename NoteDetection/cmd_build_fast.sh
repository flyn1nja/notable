# Build the NoteDetection & Maximilian lib
cd ../NoteDetection
g++ -g -Wall -D__LINUX_ALSA__ \
-I../Maximilian-master/src -I../Maximilian-master/src/libs -I../constant-q-cpp-master \
-o notabletemp.a \
   RtAudio.cpp \
   ../constant-q-cpp-master/libcq.a \
   ../Maximilian-master/src/maximilian.cpp \
   ../Maximilian-master/src/maximilianRecording.cpp \
   ../Maximilian-master/src/libs/fft.cpp ../Maximilian-master/src/libs/maxiFFT.cpp \
   -lasound -lpthread -lsndfile && \
g++ -g \
-o notable \
   player.cpp main.cpp notabletemp.a && \
   cp notable ./build/ && \
rm notable
   # QConstTrans.h FastFourierTransform.h instruments.h \