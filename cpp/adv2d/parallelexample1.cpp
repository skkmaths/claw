#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>  // Include for sin and cos functions

// Function to perform computation for non-parallel version
void perform_computation(std::vector<double>& results) {
    const int N = results.size();
    for (int i = 0; i < N; ++i) {
        results[i] = std::sin(i * 0.001) * std::cos(i * 0.002);
    }
}

int main() {
    const int N = 1000000000;  // Number of iterations
    int num_threads = 10;  // Specify the number of threads
    omp_set_num_threads(num_threads);

    // Vector to store results (optional, depending on what you do inside the loop)
    std::vector<double> results(N);

    // Start measuring time for parallel version
    auto start_parallel = std::chrono::steady_clock::now();

    // Parallelize the for loop using OpenMP
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        // Example computation inside the loop (can be replaced with actual work)
        results[i] = std::sin(i * 0.001) * std::cos(i * 0.002);
    }

    // End measuring time for parallel version
    auto end_parallel = std::chrono::steady_clock::now();

    // Calculate elapsed time for parallel version
    std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
    std::cout << "Elapsed time (Parallel): " << elapsed_seconds_parallel.count() << " seconds" << std::endl;

    // Start measuring time for non-parallel version
    auto start_sequential = std::chrono::steady_clock::now();

    // Perform the computation sequentially (non-parallel)
    perform_computation(results);

    // End measuring time for non-parallel version
    auto end_sequential = std::chrono::steady_clock::now();

    // Calculate elapsed time for non-parallel version
    std::chrono::duration<double> elapsed_seconds_sequential = end_sequential - start_sequential;
    std::cout << "Elapsed time (Sequential): " << elapsed_seconds_sequential.count() << " seconds" << std::endl;

    return 0;
}
