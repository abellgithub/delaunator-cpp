if (MSVC)
    # Since tests link CRT dynamically (/MD[d]), require gtest to
    #link dynamically too (default is /MT[d])
    option(gtest_force_shared_crt "Always use shared Visual C++ run-time DLL" ON)
endif()

set(GOOGLETEST_VERSION 1.10.0)
add_subdirectory(gtest)

