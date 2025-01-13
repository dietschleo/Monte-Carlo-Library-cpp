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

    auto start_time3 = std::chrono::high_resolution_clock::now();

    Simulation sim;
    AsianOption asian = sim.CreateAsianOption();
    asian.rand.SimulationNumber = 1000000;
/*
    auto end_time3 = std::chrono::high_resolution_clock::now();
    auto duration3 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time3 - start_time3);
    std::cout << "Execution time init: " << duration3.count() << " ms" << std::endl;

    auto start_time = std::chrono::high_resolution_clock::now();

    auto price = asian.Price();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cout << "Execution time single thread: " << duration.count() << " ms" << std::endl;
*/
    auto start_time2 = std::chrono::high_resolution_clock::now();

    auto prices = asian.TempPrice();

    auto end_time2 = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time2 - start_time2);
    std::cout << "Execution time multi thread: " << duration2.count() << " ms" << std::endl;

    



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


