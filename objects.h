#ifndef OBJECTS_H
#define OBJECTS_H

//import dependencies
#include "functions.h"
#include <map>


class EuropeanOption {
private:
    bool pricesCalculated;

public:
    double S, K, T, t, sigma, r;
    std::map<std::string, double> prices;

    //constructor
    EuropeanOption(double S, double K, double T, double t, double sigma, double r)
        : S(S), K(K), T(T), t(t), sigma(sigma), r(r), pricesCalculated(false) {}

    //getter function for prices:
    std::map<std::string, double> Price() {
        if (!pricesCalculated){
            prices = black_scholes_european(S, K, T, t, sigma, r);
            pricesCalculated = true;
        }
        return prices;
    }
};

class RandomNumber {
private:
    int seed, SimulationNumber, DayNumber; 
    double S, sigma, T;
    

public:
    std::vector<double> Zvector, Zvector2;
    std::vector<std::vector<double>> Znestedvect;

    std::vector<std::vector<double>> CreateRandomSeries(){
        Zvector = vector_std_dist(0.0, 1.0, SimulationNumber, seed);
        Zvector2 = vector_std_dist(0.0, 1.0, SimulationNumber, seed + 1);
        return {Zvector, Zvector2};
    }

    std::vector<std::vector<double>> CreateBrownianMotion(){
        Znestedvect = multi_simmulations_SBM(T, DayNumber, sigma, SimulationNumber, seed);
        return Znestedvect
    }

    RandomNumber(int seed, int SimulationNumber, int DayNumber, double S, double sigma, double T)
        : seed(seed), SimulationNumber(SimulationNumber), DayNumber(DayNumber), S(S), sigma(sigma), T(T){}
};

class Simulation{ //Factory class
private:
    int seed = 42; //instance of random number

    //define default simmulation values 
    double S = 100, K = 105, T = 1.0, t = 0.0, sigma = 0.2;
    double r = 0.05, K2 = 5, T2 = 2.0;

public: 
    int DayNumber = 252, SimulationNumber = 1000;

    explicit Simulation(){}

    RandomNumber CreateRandomNumber(std::mt19937 rng){
        RandomNumber rand(seed, SimulationNumber, DayNumber, S, sigma);
        seed = seed + 1; //change seed for future instance of RandomNumber
        return rand;
    }
    /*<
    AsianOption createAsianOption(double S, double K, double T, double sigma, double r){
        //placeholder at this point
        return AsianOption(S, K, T, sigma, r);
    }*/

    EuropeanOption CreateEuropeanOption(double S, double K, double T, double t, double sigma, double r){
        EuropeanOption euro(S, K, T, t, sigma, r); //initiate an instance of the EuropeanOption object
        return euro;
    }
};


#endif