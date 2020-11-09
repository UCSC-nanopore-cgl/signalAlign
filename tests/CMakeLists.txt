############################################################################################################
# testing
############################################################################################################
add_executable(signalAlignLibTests
        ${PROJECT_SOURCE_DIR}/tests/allTests.c
        ${PROJECT_SOURCE_DIR}/tests/eventAlignerTests.c
        ${PROJECT_SOURCE_DIR}/tests/fastaHandlerTests.c
        ${PROJECT_SOURCE_DIR}/tests/hdpTests.c
        ${PROJECT_SOURCE_DIR}/tests/nanoporeHdpTests.c
        ${PROJECT_SOURCE_DIR}/tests/randomSequences.c
        ${PROJECT_SOURCE_DIR}/tests/signalPairwiseAlignerTest.c
        ${PROJECT_SOURCE_DIR}/tests/stateMachineTests.c
        ${PROJECT_SOURCE_DIR}/tests/variableOrderPairwiseAlignerTests.c)

target_link_libraries(signalAlignLibTests PUBLIC salib)
if (NOT APPLE AND NOT BUILD_SHARED_LIBS)
    #    message(STATUS "not APPLE")
    target_link_libraries(signalAlignLibTests PUBLIC -static -static-libgcc -static-libstdc++)
endif(NOT APPLE AND NOT BUILD_SHARED_LIBS)

target_include_directories(signalAlignLibTests PUBLIC ${sonLib_inc})
add_dependencies(signalAlignLibTests lastz_make)

add_custom_command(
        TARGET signalAlignLibTests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/bin/cPecanLastz_D
        ${CMAKE_CURRENT_BINARY_DIR}/cPecanLastz_D)
add_custom_command(
        TARGET signalAlignLibTests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${PROJECT_SOURCE_DIR}/bin/cPecanLastz
        ${CMAKE_CURRENT_BINARY_DIR}/cPecanLastz)
add_custom_command(
        TARGET signalAlignLibTests POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_BINARY_DIR}/signalAlignLibTests
        ${PROJECT_SOURCE_DIR}/bin/signalAlignLibTests)


############################################################################################################
