file(GLOB TEST_SOURCES ${DEL_TEST_DIR}/*.cpp)
add_executable(unit-tests ${TEST_SOURCES})

file(GLOB TEST_SOURCES test-header-only/*.cpp)
add_executable(header-only-unit-tests ${TEST_SOURCES})
