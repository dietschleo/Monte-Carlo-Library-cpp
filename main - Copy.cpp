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
    double r;
    // Create a Simulation object
    Simulation simulation;

    // Create an AsianOption using the simulation
    AsianOption asian = simulation.CreateAsianOption();

    // Set initial parameters
    asian.K = 5;
    asian.S = 100;

    std::cout << "K " << asian.K << "T:" << asian.T << "t:" << asian.t << "sigma: " << asian.sigma << "day number:" << asian.rand.DayNumber << " sims: " << asian.rand.SimulationNumber << "\n";
    std::cout << "S Call_price Put_price call_delta put_delta call_gamma put_gamma call_vega put_vega call_vomma put_vomma\n";
    asian.Price();
    
    std::vector<std::vector<double>> matrix = multi_simmulations_SBM(1.0, 252, 0.2, 100, 0.05, 1000, 42);

    std::cout << "price1 =" << asian.prices["fixed_strike_call"];
    std::cout << "price2 =" << asian_options(100, 105, 1.0, 0.05, asian.rand.Znestedvect)["fixed_strike_call"];
    std::cout << "price3 =" << asian_options(100, 105, 1.0, 0.05, matrix)["fixed_strike_call"];
    /*for (int i = 0; i <= 50; i++) {
        asian.S = asian.S + 1; // Increment the underlying price
        asian.Greeks();
    } */

    std::cout << "\n S : " << asian.rand.S << "sigma : " << asian.rand.sigma ;

    std::cout << "computation finished!";
    std::cin >> r;
    return 0;
}

