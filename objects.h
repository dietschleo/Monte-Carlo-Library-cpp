#ifndef OBJECTS_H
#define OBJECTS_H


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Class creation///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    bool ZvectComputed, ZnestedvectComputed;
    

public:
    std::vector<double> Zvector, Zvector2;
    std::vector<std::vector<double>> Znestedvect;

    std::vector<std::vector<double>> CreateRandomSeries(){
        if (!ZvectComputed){
            Zvector = vector_std_dist(0.0, 1.0, SimulationNumber, seed);
            Zvector2 = vector_std_dist(0.0, 1.0, SimulationNumber, seed + 1);
            ZvectComputed = true;
        }
    return {Zvector, Zvector2};
    }

    RandomNumber(int seed, int SimulationNumber, int DayNumber)
        : seed(seed), SimulationNumber(SimulationNumber), DayNumber(DayNumber), ZvectComputed(false), ZnestedvectComputed(false) {}
};

class Simulation{ //Factory class
private:
    int seed; //instance of random number

public: 
    int DayNumber = 252;
    int SimulationNumber = 1000; //Nmc
    explicit Simulation(int seed) : seed(seed) {}

    RandomNumber CreateRandomNumber(std::mt19937 rng){
        RandomNumber rand(seed, SimulationNumber, DayNumber);
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