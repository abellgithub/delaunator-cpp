#include "gtest/gtest.h"

// NOTE:  This file contains main() for all tests, and therefore
//   should be included only once for each test executable.  Usually test
//   executables consist of a single source file, but if this isn't the case
//   you may need to take precautions or define your own main().

GTEST_API_ int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
