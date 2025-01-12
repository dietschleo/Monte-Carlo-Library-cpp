// main
#include <iostream>
#include <fstream>
#include <vector>
#include <optional> 
#include <map>
#include <cmath>
#include "functions.h"
#include "objects.h"
#include <tuple>
#include <string>
#include <random>
void test_asian_options() {
    // Parameters for the test
    double S = 100.0;        // Initial stock price
    double r = 0.05;         // Risk-free rate
    double T = 1.0;          // Maturity (1 year)
    double K = 100.0;        // Strike price
    double sigma = 0.2;      // Volatility
    int N = 10;             // Number of time steps
    int Nmc = 10;         // Number of Monte Carlo simulations
    int seed = 42;           // Random seed for reproducibility
/*
    // Step 1: Generate GBM trajectories
    auto underlying = multi_simmulations_SBM(T, N, sigma, S, r, Nmc, seed);

    // Step 2: Calculate Asian option prices
    auto results = asian_options(S, r, T, K, underlying);

    // Step 3: Print results
    std::cout << "Asian Option Prices:\n";
    std::cout << "Fixed Strike Call: " << results["fixed_strike_call"] << "\n";
    std::cout << "Fixed Strike Put: " << results["fixed_strike_put"] << "\n";
    std::cout << "Floating Strike Call: " << results["floating_strike_call"] << "\n";
    std::cout << "Floating Strike Put: " << results["floating_strike_put"] << "\n";
*/

   Simulation sim;
   AsianOption asian = sim.CreateAsianOption();
   asian.S = S;
   asian.r = r;
   asian.T = T;
   asian.K = K;
   asian.sigma = sigma;
   asian.rand.sigma = sigma;
   asian.rand.S = S;
   asian.rand.DayNumber = N;
   asian.rand.SimulationNumber = Nmc;


 auto prices = asian.Price();
    std::cout << "Asian Option Prices:\n";
    std::cout << "Fixed Strike Call: " << prices["fixed_strike_call"] << "\n";
    std::cout << "Fixed Strike Put: " << prices["fixed_strike_put"] << "\n";
    std::cout << "Floating Strike Call: " << prices["floating_strike_call"] << "\n";
    std::cout << "Floating Strike Put: " << prices["floating_strike_put"] << "\n";

auto greeks = asian.ExtractGreeks();

for (const auto& pair : greeks) {
    std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
}



}


// Assume all the required function declarations are included here, such as:
// multi_simmulations_SBM, asian_options, delta_gamma_extraction, vega_vomma_extraction
int main() {
double r;

test_asian_options();

    std::cout << "computation finished!";
    std::cin >> r;
    return 0;
}

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <string>

// Assuming all functions from your provided code are already defined.


