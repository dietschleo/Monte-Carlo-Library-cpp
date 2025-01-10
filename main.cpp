// main
#include <iostream>
#include <fstream>
#include <vector>
#include <optional> 
#include <map>
#include <cmath>
#include "functions.h"
#include <tuple>
#include <string>
#include <random>

// Assume all the required function declarations are included here, such as:
// multi_simmulations_SBM, asian_options, delta_gamma_extraction, vega_vomma_extraction

int main() {


    // Define inputs for the Asian option simulation
    std::string option_type = "asian_options";

    double S = 75.0;       // Initial stock price
    double t = 0.0; // Current time
    double K = 100.0;       // Strike price
    double T = 1.0;         // Time to maturity (in years)
    double sigma = 0.2;     // Volatility
    double r = 0.05;        // Risk-free interest rate
    double N = 252;            // Number of time steps
    int Nmc = 1000000;        // Number of Monte Carlo simulations
    double v_h = 0.01;      // Perturbation for volatility
    double s_h = 0.5;       // Perturbation for stock price
    int seed = 42;

    // Create a Simulation object
    Simulation simulation(seed);

    // Create a EuropeanOption using the simulation
    EuropeanOption europeanOption = simulation.createEuropeanOption(S, K, T, t, sigma, r);

    // Get the prices of the EuropeanOption
    std::map<std::string, double> prices = europeanOption.Price();

    // Print the results
    std::cout << "European Option Prices:" << std::endl;
    for (const auto& price : prices) {
        std::cout << price.first << ": " << price.second << std::endl;
    }



    std::cin >> r;
    return 0;
}

