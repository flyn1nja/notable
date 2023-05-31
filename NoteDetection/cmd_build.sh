g++ -Wall -D__LINUX_ALSA__ \
-I../Maximilian-master/src -I../Maximilian-master/src/libs \
-o maximilian player.cpp main.cpp RtAudio.cpp ../Maximilian-master/src/maximilian.cpp -lasound -lpthread && \
cp -i maximilian ./build/
rm maximilian