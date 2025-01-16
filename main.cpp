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
#include <chrono> // For timing


int main() {
    double r;
    auto start_time = std::chrono::high_resolution_clock::now();


    Simulation sim;
    sim.SimulationNumber = 100000;
    
    auto eur = sim.CreateEuropeanOption();
    std::cout << "eur.rand.SimulationNumber "<< eur.rand.SimulationNumber << std::endl;
    eur.Analytical(true);

    auto p = eur.MonteCarloPrice();
    auto p2= eur.AnalyticalPrice();
    std::cout << "mcPrice: " << p["call"] << "  P_analytical: " << p2["call"];

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Execution time multi thread: " << duration.count() << " ms" << std::endl;
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


