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

//test_asian_options();
    Simulation sim;

    for(int i = 75; i <= 125; ++i){
    AmericanOption US = sim.CreateAmericanOption();
    US.S = i;
    US.rand.S = i;
    double price = US.CallPrice();
    std::cout << "S= "<< US.rand.S << "|  price: " << price << "\n";
    }


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


