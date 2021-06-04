add_library(
  ActsExamplesIoGraphCreation SHARED
  src/ModuleMapCreatorWriter.cpp
  src/GraphCreatorWriter.cpp
  src/GraphCreationOptions.cpp
  )

target_include_directories(
  ActsExamplesIoGraphCreation
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
  
target_link_libraries(
  ActsExamplesIoGraphCreation
  PRIVATE
    ActsExamplesIoCsv
    ActsExamplesFramework ${PYTHON_LIBRARIES}
    Threads::Threads Boost::program_options dfelibs
)
    
target_include_directories(ActsExamplesIoGraphCreation PRIVATE ${PYTHON_INCLUDE_DIRS})

install(
  TARGETS ActsExamplesIoGraphCreation
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
