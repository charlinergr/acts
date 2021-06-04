#GraphCreator
add_executable(
    ActsRecGNNTracks
    RecGNNTracks.cpp)
target_link_libraries(
    ActsRecGNNTracks
    PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDetectorGeneric
    ActsExamplesIoCsv
    ActsExamplesIoGraphCreation
)


#ModuleMap creation
add_executable(
    ActsRecModuleMapCreation
    ModuleMapCreation.cpp)
target_link_libraries(
    ActsRecModuleMapCreation
    PRIVATE
    ActsExamplesFramework
    ActsExamplesCommon
    ActsExamplesDetectorGeneric
    ActsExamplesIoCsv
    ActsExamplesIoGraphCreation
)
install(
  TARGETS
   ActsRecGNNTracks
   ActsRecModuleMapCreation
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

