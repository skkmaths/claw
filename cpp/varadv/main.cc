#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>  // Include for sin and cos functions
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <string>
#include <iomanip> // for std::setprecision
#include <cstring> // for strcmp
// This code solves the linear advection equation in 2D 
// u_t + u_x + u_y = 0
#include "twodproblem.h"

using namespace std;
int main(int argc, char* argv[]) {
    int nx = 50; // Default value for nx
    int ny = 50; // Default value for ny
    double Tf = 1.0; // Default value for final time
    double cfl = 0.4; // Default value of cfl 
    unsigned int save_freq = 10; // Default save frequency
    string scheme = "fo"; //Default:first order in time and space scheme
    // Parse command-line arguments if provided
    // Pass this argument to the exe file 
    // ./twodproblem 50 50 2
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-nx") == 0 && i + 1 < argc) {
            nx = std::atoi(argv[i + 1]);
            i++; // Skip the next argument
        } else if (strcmp(argv[i], "-ny") == 0 && i + 1 < argc) {
            ny = std::atoi(argv[i + 1]);
            i++; // Skip the next argument
        } else if (strcmp(argv[i], "-Tf") == 0 && i + 1 < argc) {
            Tf = std::atof(argv[i + 1]);
            i++; // Skip the next argument
        } else if (strcmp(argv[i], "-cfl") == 0 && i + 1 < argc) {
            cfl = std::atof(argv[i + 1]);
            i++; // Skip the next argument
        } else if (strcmp(argv[i], "-save_freq") == 0 && i + 1 < argc) {
            save_freq = std::atof(argv[i + 1]);
            i++; // Skip the next argument
        } else if (strcmp(argv[i], "-scheme") == 0 && i + 1 < argc) {
            scheme =argv[i + 1];
            i++; // Skip the next argument
        } else {
            std::cerr << "Unknown or incomplete argument: " << argv[i] << std::endl;
            return 1;
        }
    }

    TwoDProblem twodproblem(nx, ny, Tf, cfl, save_freq, scheme);
    twodproblem.run();
    return 0;
}