#ifndef __QUADRATURE_H__
#define __QUADRATURE_H__

#include <string>
#include <cmath>
#include <math.h>
#include <limits>
#include <algorithm>
#include <vector>
#include<iostream>
#include<stdio.h>
#define SZERO 0.001

using namespace std;

// GaussLobatto quadraature rule

long double gamma(const unsigned int n);


long double JacobiP(const long double x,
                                      const int alpha,
                                        int beta,
                                      const unsigned int n);

std::vector<long double>
compute_quadrature_points(const unsigned int q,
                          const int alpha,
                          const int beta);



std::vector<long double>
compute_quadrature_weights(const std::vector<long double> &x,
                           const int alpha,
                           const int beta);


void get_points_weights(std::vector<long double>& q_p, std::vector<long double>& q_w, const unsigned int n,
		const double& left_boundary, const double& right_boundary);


#endif

