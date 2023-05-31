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

# __LINUX_ALSA__,     /*!< The Advanced Linux Sound Architecture API. */
# __UNIX_JACK__,      /*!< The Jack Low-Latency Audio Server API. */
# __LINUX_PULSE__,    /*!< The Linux PulseAudio API. */
# __LINUX_OSS__,      /*!< The Linux Open Sound System API. */

add_definitions(-D__LINUX_ALSA__) # Soundcards, self-installed audio, ...
add_compile_definitions(__LINUX_ALSA__)
target_link_libraries(maximilian PUBLIC -lpthread)
target_link_libraries(maximilian PUBLIC -lasound)

# add_definitions(-D__LINUX_PULSE__) # Default for Ubuntu, Mint, openSUSE, BSD, ...
# target_link_libraries(maximilian PUBLIC -lpthread)
# target_link_libraries(maximilian PUBLIC -lpulse-simple)
# target_link_libraries(maximilian PUBLIC -lpulse)

endif()

if (OSX)
  MESSAGE(STATUS "OSX build")

  # __MACOSX_CORE__,    /*!< Macintosh OS-X Core Audio API. */

  add_definitions(-D__MACOSX_CORE__)
  add_compile_definitions(__MACOSX_CORE__)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__MACOSX_CORE__")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__MACOSX_CORE__")

  target_link_libraries(maximilian PUBLIC -lpthread)
  target_link_libraries(maximilian PUBLIC -lasound)

  #incomplete - need to add frameworks -framework CoreAudio -framework CoreFoundation
endif()

if(WIN32)
  MESSAGE(STATUS "Windows build")

  # __WINDOWS_ASIO__,   /*!< The Steinberg Audio Stream I/O API. */
  # __WINDOWS_WASAPI__, /*!< The Microsoft WASAPI API. */
  # __WINDOWS_DS__,     /*!< The Microsoft DirectSound API. */

  add_definitions(-D__WINDOWS_DS__)
  add_compile_definitions(__WINDOWS_DS__)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__WINDOWS_DS__")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__WINDOWS_DS__")

  target_link_libraries(maximilian PUBLIC -lole32)
  target_link_libraries(maximilian PUBLIC -lwinmm)
  target_link_libraries(maximilian PUBLIC -ldsound)
endif(WIN32)

