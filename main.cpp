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

// Assume all the required function declarations are included here, such as:
// multi_simmulations_SBM, asian_options, delta_gamma_extraction, vega_vomma_extraction

int main() {



    double S = 75.0;       // Initial stock price
    double t = 0.0; // Current time
    double K = 100.0;       // Strike price
    double T = 1.0;         // Time to maturity (in years)
    double sigma = 0.2;     // Volatility
    double r = 0.05;        // Risk-free interest rate
    double N = 252;            // Number of time steps
    int Nmc = 10000;        // Number of Monte Carlo simulations
    double v_h = 0.01;      // Perturbation for volatility
    double s_h = 0.5;       // Perturbation for stock price
    int seed = 42;

    // Create a Simulation object
    Simulation simulation;

    // Create a EuropeanOption using the simulation
    EuropeanOption europeanOption = simulation.CreateEuropeanOption();

    // Get the prices of the EuropeanOption
    std::map<std::string, double> prices = europeanOption.Price();
    std::map<std::string, double> mcprices = europeanOption.MonteCarloPrice();


    std::cout << "call: " << prices["call"] << std::endl;

    std::cout << "mc call: " << mcprices["call"] << std::endl;

    std::cout << "put: " << prices["put"] << std::endl;

    std::cout << "mc put: " << mcprices["call"] << std::endl;



    std::cin >> r;
    return 0;
}

