function(d_add_test _testfile)
    string(REGEX REPLACE "\\.[^.]*$" "" _exe ${_testfile})
    add_executable(${_exe} ${_testfile})
    d_target_compile_settings(${_exe})
    target_include_directories(
        ${_exe}
            PRIVATE
                ${D_INCLUDE_DIR}
                ${D_GTEST_INCLUDE_DIR}
    )
    target_link_libraries(${_exe}
        PRIVATE
            delaunator
            gtest
    )
    add_test(NAME ${_exe}
        COMMAND 
            ${PROJECT_BINARY_DIR}/bin/${_exe}
        WORKING_DIRECTORY
            ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/..
    )
endfunction()

enable_testing()
d_add_test(issue-8.cpp)
d_add_test(delaunator.test.cpp)
