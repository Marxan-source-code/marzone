#include "CppUTest/CommandLineTestRunner.h"
#include "../logger.hpp"

namespace marzone {
    Logger logger;
}
int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}