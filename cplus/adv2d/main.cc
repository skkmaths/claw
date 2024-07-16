#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>  // Include for sin and cos functions
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <iomanip> // for std::setprecision

// This code solves the linear advection equation in 2D 
// of the form u_t + u_x + u_y = 0
#include "twodproblem.h"

using namespace std;
int main() {

    TwoDProblem twodproblem;
    twodproblem.run();
    return 0;
}