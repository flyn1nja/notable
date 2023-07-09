# Build the CQT lib
cd ../constant-q-cpp-master
make -f Makefile.linux

# Build the NoteDetection & Maximilian lib
cd ../NoteDetection
g++ -Wall -D__LINUX_ALSA__ \
-I../Maximilian-master/src -I../Maximilian-master/src/libs \
-o maximilian \
   player.cpp main.cpp RtAudio.cpp \
   QConstTrans.cpp \
   ../Maximilian-master/src/maximilian.cpp \
   ../Maximilian-master/src/libs/fft.cpp ../Maximilian-master/src/libs/maxiFFT.cpp \
   -lasound -lpthread && \
cp maximilian ./build/
   # ../Maximilian-master/src/maximilianRecording.cpp \
   # QConstTrans.h FastFourierTransform.h instruments.h \
rm maximilian