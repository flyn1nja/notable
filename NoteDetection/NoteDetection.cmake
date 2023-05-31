

# include_directories(.)
include_directories(../Maximilian-master/src)
include_directories(../Maximilian-master/src/libs)

SET (MAXI_SRC ./main.cpp ./player.cpp ./RtAudio.cpp ../Maximilian-master/src/maximilian.cpp)

add_executable(maximilian  ${MAXI_SRC} ${MAXI_SRC_EXTENDED})

target_compile_options(maximilian PUBLIC -Wall)

if (UNIX AND NOT APPLE)
  SET(LINUX TRUE)
endif()

if (LINUX)
  MESSAGE(STATUS "Linux build")

  add_definitions(-D__LINUX_ALSA__) # Soundcards, self-installed audio, ...
  add_compile_definitions(__LINUX_ALSA__)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__LINUX_ALSA__")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__LINUX_ALSA__")

  target_link_libraries(maximilian PUBLIC -lpthread)
  target_link_libraries(maximilian PUBLIC -lasound)
  
  # add_definitions(-D__LINUX_PULSE__) # Native to Ubuntu, Mint, openSUSE, BSD, ...
  # target_link_libraries(maximilian PUBLIC -lpthread)
  # target_link_libraries(maximilian PUBLIC -lpulse-simple)
  # target_link_libraries(maximilian PUBLIC -lpulse)

endif()

if (APPLE)
  MESSAGE(STATUS "OSX BUILD")

  add_definitions(-D__MACOSX_CORE__)
  add_compile_definitions(__MACOSX_CORE__)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__MACOSX_CORE__")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__MACOSX_CORE__")

  target_link_libraries(maximilian PUBLIC -lpthread)
  find_library(CA CoreAudio)
  find_library(CF CoreFoundation)
  target_link_libraries(maximilian PUBLIC ${CA} ${CF})
endif()

if (WIN32)
#TODO
  MESSAGE(STATUS "Windows BUILD")

  add_definitions(-D__WINDOWS_DS__)
  add_compile_definitions(__WINDOWS_DS__)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__WINDOWS_DS__")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__WINDOWS_DS__")

  target_link_libraries(maximilian PUBLIC -lole32)
  target_link_libraries(maximilian PUBLIC -lwinmm)
  target_link_libraries(maximilian PUBLIC -ldsound)
endif()

