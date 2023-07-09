# Build the CQT lib
cd ../constant-q-cpp-master
make -f Makefile.linux

# Build the NoteDetection & Maximilian lib
cd ../NoteDetection
g++ -Wall -D__LINUX_ALSA__ \
-I../Maximilian-master/src -I../Maximilian-master/src/libs -I../constant-q-cpp-master \
-o maximilian \
   player.cpp main.cpp RtAudio.cpp \
   ../constant-q-cpp-master/libcq.a \
   ../Maximilian-master/src/maximilian.cpp \
   ../Maximilian-master/src/maximilianRecording.cpp \
   ../Maximilian-master/src/libs/fft.cpp ../Maximilian-master/src/libs/maxiFFT.cpp \
   -lasound -lpthread -lsndfile && \
cp maximilian ./build/ && \
rm maximilian
   # QConstTrans.h FastFourierTransform.h instruments.h \