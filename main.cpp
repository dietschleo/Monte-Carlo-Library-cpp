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
    double K = 100.0;       // Strike price
    double T = 1.0;         // Time to maturity (in years)
    double sigma = 0.2;     // Volatility
    double r = 0.05;        // Risk-free interest rate
    double N = 252;            // Number of time steps
    int Nmc = 1000000;        // Number of Monte Carlo simulations
    double v_h = 0.01;      // Perturbation for volatility
    double s_h = 0.5;       // Perturbation for stock price

    // Create the arguments vector (must align with the function's expected order)
    std::vector<double> arguments = {
        S,          // Stock price
        K,          // Strike price
        T,          // Time to maturity
        0.0,        // Placeholder for `t` (unused in this example)
        sigma,      // Volatility
        r,          // Risk-free rate
        0.0,        // Placeholder for `K2` (unused in this example)
        0.0,        // Placeholder for `T1` (unused in this example)
        0.0,        // Placeholder for `T2` (unused in this example)
        N // Number of time steps
    };

    std::cout << "S,fixed_strike_call_delta,fixed_strike_call_gamma,fixed_strike_call_price,fixed_strike_call_vega,fixed_strike_call_vomma,fixed_strike_put_delta,fixed_strike_put_gamma,fixed_strike_put_price,fixed_strike_put_vega,fixed_strike_put_vomma,floating_strike_call_delta,floating_strike_call_gamma,floating_strike_call_price,floating_strike_call_vega,floating_strike_call_vomma,floating_strike_put_delta,floating_strike_put_gamma,floating_strike_put_price,floating_strike_put_vega,floating_strike_put_vomma";
    for (int i=0; i<=50; ++i){
        arguments[0]=S+ static_cast<double>(i);
        // Call the Monte Carlo simulation function
        std::map<std::string, double> results = monte_carlo_simmulation(option_type, arguments, Nmc, v_h, s_h);

        // Display the results
        //std::cout << "Monte Carlo Simulation Results for Asian Options:\n" ;

        if (results.empty()) {
            std::cout << "Results map is empty!" << std::endl;
        }
        std::cout << arguments[0] << ",";
        for (const auto& result : results) {
            std::cout << result.second << ",";
        }
    std::cout << "\n";
    }

    std::cin >> r;

    return 0;
}

