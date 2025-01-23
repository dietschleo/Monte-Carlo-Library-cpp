#include <iostream>
#include <vector>
#include <map>
#include "Simulation.h"
#include "RandomNumber.h"

int main() {
    
    // Initialize Simulation
    Simulation sim;
    sim.SimulationNumber = 100000;

    // Enable multi-threading
    sim.n_threads = 4;

    // Create Asian options
    auto option = sim.CreateAsianOption();



    // Iterate over S values
    for (int i = 90; i <= 110; ++i) {
        option.S = i;

        // Multi-threaded pricing
        auto greeks = option.ExtractGreeks();
        std::cout << option.S << " : fixed strike call price: " << option.prices["fixed_strike_call"]
            << ", corresponding vega: " << option.greeks["fixed_strike_call_vega"]<< std::endl;

        // Clear memory and reset for next iteration
        option.Clear();
        option.rand.Reset();
    }

    std::cin >> sim.n_threads;
    return 0;
}
