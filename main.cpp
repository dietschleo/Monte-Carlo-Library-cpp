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

    // Create a EuropeanOption using the simulation
    EuropeanOption euro = simulation.CreateEuropeanOption();

    for (int i = 0; i <= 50; i++){

    euro.S = euro.S + 1;

    // Get the prices of the EuropeanOption
    std::map<std::string, double> prices = euro.Price();
    
    for (const auto& [key, value] : prices) {
        std::cout  << key << ": " << value << "\n" << std::endl;
    }

    prices = euro.Greeks();

    for (const auto& [key, value] : prices) {
        std::cout << "S=" << euro.S << " " << key << ": " << value << std::endl;
    }
    
    }

    std::cout << "computation finished!";
    std::cin >> r;
    return 0;
}

