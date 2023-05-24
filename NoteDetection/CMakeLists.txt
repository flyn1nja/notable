cmake_minimum_required (VERSION 3.15)

project (maximilian)

set(CMAKE_CXX_STANDARD 17)

include_directories(../Maximilian-master/src)
include_directories(../Maximilian-master/src/libs)

# SET(MAXI_SRC_EXTENDED ../Maximilian-master/src/libs/fft.cpp ../Maximilian-master/src/libs/maxiFFT.cpp)

add_executable(maximilian main.cpp player.cpp RtAudio.cpp ../Maximilian-master/src/maximilian.cpp)
target_compile_options(maximilian PUBLIC -Wall)

if (LINUX)
  MESSAGE(STATUS "Linux build")
  target_link_libraries(maximilian PUBLIC -lpthread)
  target_link_libraries(maximilian PUBLIC -lasound)
  add_definitions(-D__LINUX_ALSA__)
endif()

if (OSX)
  MESSAGE(STATUS "OSX build")
  target_link_libraries(maximilian PUBLIC -lpthread)
  target_link_libraries(maximilian PUBLIC -lasound)
  add_definitions(-D__MACOSX_CORE__)
  #incomplete - need to add frameworks -framework CoreAudio -framework CoreFoundation
endif()

if(WIN32)
  MESSAGE(STATUS "Windows build")
  target_link_libraries(maximilian PUBLIC -lole32)
  target_link_libraries(maximilian PUBLIC -lwinmm)
  target_link_libraries(maximilian PUBLIC -ldsound)
  add_definitions(-D__WINDOWS_DS__)
endif(WIN32)

