// main
#include <iostream>
#include <fstream>
#include <vector>
#include <optional> 
#include <map>
#include <cmath>
#include "functions.h"
#include "objects.h"
#include "RandomNumber.h"
#include <tuple>
#include <string>
#include <random>
#include <chrono> // For timing


int main() {
    double r;
    auto start_time = std::chrono::high_resolution_clock::now();


    Simulation sim;
    sim.SimulationNumber = 1000;
    //sim.n_threads = 6;
    
    //sim.K = 5;
   auto option = sim.CreateAsianOption();
   auto optiono = sim.CreateAsianOption();
   
    option.n_threads = 1;
    optiono.n_threads = 6;
    std::cout << "S price delta gamma vega vomma \n";
    for (int i = 90; i <= 110; ++i){
        option.S = i;
        optiono.S = i;
        option.rand.S = i;
        optiono.rand.S = i;
        auto greeks = option.ExtractGreeks();
        std::cout << option.S << " single thread " << option.prices["fixed_strike_call"] << "\n";
        option.n_threads = 6;
        auto greekos = optiono.ExtractGreeks();
        std::cout << option.S << " multi  thread " << optiono.prices["fixed_strike_call"] << "\n";
        option.Clear();
        optiono.Clear();
        option.rand.Reset();
//        std::cout <<i << " " << option.prices["call"] << " " << greeks["call_delta"] << " " << greeks["call_gamma"] << " " << greeks["call_vega"] << " " << greeks["call_vomma"] << "\n";
//        std::cout <<i << " " << option.prices["call"] << " " << option.greeks["call_delta"] << " " << option.greeks["call_gamma"] << " " << option.greeks["call_vega"] << " " << option.greeks["call_vomma"] << "\n";

    }
    

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
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


