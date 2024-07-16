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
int main(int argc, char* argv[]) {
    int nx = 50; // Default value for nx
    int ny = 50; // Default value for ny
    double Tf = 2.0; // Default value for Tf
    double cfl = 0.4; // Default value of cfl 
    // Parse command-line arguments if provided
    // Pass this argument to the exe file 
    // ./twodproblem 50 50 2
    if (argc > 1) {
        nx = std::atoi(argv[1]);
    }
    if (argc > 2) {
        ny = std::atoi(argv[2]);
    }
    if (argc > 3) {
        Tf = std::atof(argv[3]);
    }
    if (argc > 4) {
        cfl = std::atof(argv[4]);
    }
    
    TwoDProblem twodproblem(nx, ny, Tf, cfl);
    twodproblem.run();
    return 0;
}