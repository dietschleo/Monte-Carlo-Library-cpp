#ifndef OBJECTS_H
#define OBJECTS_H

//import dependencies
#include "functions.h"
#include <map>

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
        if (Znestedvect.size() == 0){
            Znestedvect = multi_simmulations_SBM(T, DayNumber, sigma, SimulationNumber, seed);
        }
        return Znestedvect;
    }

    RandomNumber(int seed, int SimulationNumber, int DayNumber, double S, double sigma, double T)
        : seed(seed), SimulationNumber(SimulationNumber), DayNumber(DayNumber), S(S), sigma(sigma), T(T){}
};

class EuropeanOption {
private:
    bool pricesCalculated;
    RandomNumber rand;
public:
    double S, K, T, t, sigma, r;
    std::map<std::string, double> prices;

    //constructor
    EuropeanOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r)
        : rand(rand), S(S), K(K), T(T), t(t), sigma(sigma), r(r), pricesCalculated(false) {}

    //getter function for prices:
    std::map<std::string, double> Price() {
        prices = black_scholes_european(S, K, T, t, sigma, r);
        return prices;
    }

    std::map<std::string, double> MonteCarloPrice() {
        if (!pricesCalculated){
            rand.CreateRandomSeries();
            prices = monte_carlo_european(S, K, T, t, sigma, r, rand.Zvector);
            pricesCalculated = true;
        }
        return prices;
    }
};



class Simulation{ //Factory class
private:
    int seed = 42; //instance of random number

    //define default simmulation values 
    double S = 100, K = 105, T = 1.0, t = 0.0, sigma = 0.2;
    double r = 0.05, K2 = 5, T2 = 2.0, h_v = 0.01, h_s = 0.2 ;

    RandomNumber rand;

public: 
    int DayNumber = 252, SimulationNumber = 1000;

    explicit Simulation()
        : rand(seed, SimulationNumber, DayNumber, S, sigma, T){} //constructor automatically creates instance rand

    //RandomNumber(int seed, int SimulationNumber, int DayNumber, double S, double sigma, double T)
    RandomNumber CreateRandomNumber(int seed){
        RandomNumber newrand(seed, SimulationNumber, DayNumber, S, sigma, T);
        seed = seed + 1; //change seed for future instance of RandomNumber
        return newrand;
    }
    /*<
    AsianOption createAsianOption(double S, double K, double T, double sigma, double r){
        //placeholder at this point
        return AsianOption(S, K, T, sigma, r);
    }*/

    EuropeanOption CreateEuropeanOption(){
        EuropeanOption euro(rand, S, K, T, t, sigma, r); //initiate an instance of the EuropeanOption object
        return euro;
    }
};


#endif