

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
  target_link_libraries(maximilian PUBLIC -lpthread)
  target_link_libraries(maximilian PUBLIC -lasound)
  add_definitions(-D__LINUX_ALSA__)
endif()

if (APPLE)
  MESSAGE(STATUS "OSX BUILD")
  target_link_libraries(maximilian PUBLIC -lpthread)
  find_library(CA CoreAudio)
  find_library(CF CoreFoundation)
  target_link_libraries(maximilian PUBLIC ${CA} ${CF})
  add_definitions(-D__MACOSX_CORE__)
endif()

if (WIN32)
#TODO
  MESSAGE(STATUS "Windows BUILD")
  target_link_libraries(maximilian PUBLIC -lole32)
  target_link_libraries(maximilian PUBLIC -lwinmm)
  target_link_libraries(maximilian PUBLIC -ldsound)
  add_definitions(-D__WINDOWS_DS__)
endif()

