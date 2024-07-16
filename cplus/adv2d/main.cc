#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>  // Include for sin and cos functions
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <iomanip> // for std::setprecision


#include "twodproblem.h"

using namespace std;
int main() {

    TwoDProblem twodproblem;
    twodproblem.run();
    return 0;
}