#ifndef OBJECTS_H
#define OBJECTS_H

//import dependencies
#include "functions.h"
#include <map>

class RandomNumber {
private:
    int seed;    
    double S_hist, sigma_hist, r_hist; //control variables for updates
public:
    double S, sigma, T, r;
    int SimulationNumber, DayNumber;
    std::vector<double> Zvector, Zvector2;
    std::vector<std::vector<double>> Znestedvect;

    std::vector<std::vector<double>> CreateRandomSeries(){
        Zvector = vector_std_dist(0.0, sigma, SimulationNumber, seed);
        Zvector2 = vector_std_dist(0.0, sigma, SimulationNumber, seed + 1);
        return {Zvector, Zvector2};
    }

    std::vector<std::vector<double>> CreateBrownianMotion(){
        //check if the vector has never been computed or variables have been updated since last computation
        if (Znestedvect.size() == 0 || S != S_hist || sigma != sigma_hist || r != r_hist){
            Znestedvect = multi_simmulations_SBM(T, DayNumber, sigma, S, r, SimulationNumber, seed);
            S_hist = S;
            sigma_hist= sigma;
            r_hist=r;
        }
        return Znestedvect;
    }

    void Reset(){
        //method to clear memory of instance
        Zvector.clear();
        Zvector2.clear();
        Znestedvect.clear();
    }

    RandomNumber(int seed, int SimulationNumber, int DayNumber, double S, double r, double sigma, double T)
        : seed(seed), SimulationNumber(SimulationNumber), DayNumber(DayNumber), S(S), sigma(sigma), r(r), T(T){}
};

class EuropeanOption {
private:
    bool pricesCalculated, MCdefault = false;
    std::map<std::string, double> MCprice; //cache for MC calculated price
public:
    RandomNumber rand;
    double S, K, T, t, sigma, r, s_h, v_h;
    std::map<std::string, double> prices, greeks;

    //constructor
    EuropeanOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r,double v_h,double s_h)
        : rand(rand), S(S), K(K), T(T), t(t), sigma(sigma), r(r), v_h(v_h), s_h(s_h), pricesCalculated(false) {}

    //function for prices:

    std::map<std::string, double> Price() {
        if (!MCdefault){
            prices = AnalyticalPrice();
        } else {
            prices = MonteCarloPrice();
        }
        return prices;
    }
    std::map<std::string, double> AnalyticalPrice() {
        prices = black_scholes_european(S, K, T, t, sigma, r);
        return prices;
    }

    std::map<std::string, double> MonteCarloPrice() {
        if (MCprice.empty()){
            rand.CreateRandomSeries();
            prices = monte_carlo_european(S, K, T, t, sigma, r, rand.Zvector);
            MCprice = prices; 
            pricesCalculated = true;
        }
        return prices;
    }

    void Analytical(bool analytical){
        // set setting to use monte carlo by default
        if(!analytical){
            MCdefault = false;
        } else {
            MCdefault = true;
        }
    }

    std::map<std::string, double> ExtractGreeks(){
        if (prices.empty()){
            prices = Price();
        }
        //placeholders for current variables and prices
        std::map<std::string, double> temp_prices = prices;
        double temp_S = S;
        double temp_sigma = sigma;
        S = S + s_h;
        std::map<std::string, double> p_plus = Price();
        S = temp_S - s_h;
        std::map<std::string, double> p_minus = Price();
        sigma = sigma + v_h;
        std::map<std::string, double>  p_plus_v = Price();
        sigma = temp_sigma - v_h;
        std::map<std::string, double>  p_minus_v = Price();

        prices = temp_prices;
        std::map<std::string, double> results;
        for (auto& pair : prices) { //for each type of contract:
            results[pair.first + "_price"] = pair.second;
            std::map<std::string, double> greeks = delta_gamma_extraction(prices[pair.first], p_plus[pair.first], p_minus[pair.first], s_h);
            greeks["vega"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vega"];
            greeks["vomma"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vomma"];
            for (auto& duo : greeks) { //for each greek:
                std::string contract_greek = pair.first + "_" + duo.first;
                results[contract_greek] = duo.second;
                //std::cout << "Value for " << contract_greek << " is : " << results[contract_greek] << "\n";
            }
        }
        //reset original price and variable
        sigma = temp_sigma;
        S = temp_S;
        return results;
    }


};

class CompoundOption {
private:
    bool pricesCalculated;
//compound_option(double S, double K1, double K2, double T1, double T2, double sigma, double r, const std::vector<double>& z1, const std::vector<double>& z2)
public:
    RandomNumber rand;
    double S, K, K2, T, T2, t, sigma, r, s_h, v_h;
    std::map<std::string, double> prices;

    CompoundOption(RandomNumber rand, double S, double K, double K2, double T, double T2, double t, double sigma, double r,double v_h,double s_h)
    : rand(rand), S(S), K(K), K2(K2), T(T), T2(T2), t(t), sigma(sigma), r(r), v_h(v_h), s_h(s_h), pricesCalculated(false) {}

    std::map<std::string, double> Price(){
        if (!pricesCalculated){
            rand.CreateRandomSeries();
            prices = compound_option(S, K, K2, T, T2, sigma, r, rand.Zvector, rand.Zvector2);
        }
        return prices;
    }

    std::map<std::string, double> ExtractGreeks(){
        if (prices.empty()){
            prices = Price();
        }
        //placeholders for current variables and prices
        std::map<std::string, double> temp_prices = prices;
        double temp_S = S;
        double temp_sigma = sigma;
        S = S + s_h;
        std::map<std::string, double> p_plus = Price();
        S = temp_S - s_h;
        std::map<std::string, double> p_minus = Price();
        sigma = sigma + v_h;
        std::map<std::string, double>  p_plus_v = Price();
        sigma = temp_sigma - v_h;
        std::map<std::string, double>  p_minus_v = Price();

        prices = temp_prices;
        std::map<std::string, double> results;
        for (auto& pair : prices) { //for each type of contract:
            results[pair.first + "_price"] = pair.second;
            std::map<std::string, double> greeks = delta_gamma_extraction(prices[pair.first], p_plus[pair.first], p_minus[pair.first], s_h);
            greeks["vega"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vega"];
            greeks["vomma"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vomma"];
            for (auto& duo : greeks) { //for each greek:
                std::string contract_greek = pair.first + "_" + duo.first;
                results[contract_greek] = duo.second;
                //std::cout << "Value for " << contract_greek << " is : " << results[contract_greek] << "\n";
            }
        }
        //reset original price and variable
        
        sigma = temp_sigma;
        S = temp_S;
        return results;
    }
    
};

class AsianOption {
private:
    bool pricesCalculated;
    

public:
    RandomNumber rand;
    double S, K, T, t, sigma, r, v_h, s_h;
    std::map<std::string, double> prices, greeks;

    AsianOption(RandomNumber rand, double S, double K, double T, double t, double sigma, double r, double v_h, double s_h)
        : rand(rand), S(S), K(K), T(T), t(t), sigma(sigma), r(r), v_h(v_h), s_h(s_h), pricesCalculated(false) {}
    
    std::map<std::string, double> Price(){
        rand.S = S;
        rand.sigma = sigma;
        rand.CreateBrownianMotion();
        prices = asian_options(S, r, T, K, rand.Znestedvect);
    return prices;
    }

        std::map<std::string, double> ExtractGreeks(){
        //placeholders for current variables and prices
        if (prices.empty()){
            prices = Price();
        }
        std::map<std::string, double> temp_prices = prices;
        double temp_S = S;
        double temp_sigma = sigma;
        S = S + s_h;
        std::map<std::string, double> p_plus = Price();
        S = temp_S - s_h;
        std::map<std::string, double> p_minus = Price();
        sigma = sigma + v_h;
        std::map<std::string, double>  p_plus_v = Price();
        sigma = temp_sigma - v_h;
        std::map<std::string, double>  p_minus_v = Price();

        prices = temp_prices;
        std::map<std::string, double> results;
        for (auto& pair : prices) { //for each type of contract:
            results[pair.first + "_price"] = pair.second;
            std::map<std::string, double> greeks = delta_gamma_extraction(prices[pair.first], p_plus[pair.first], p_minus[pair.first], s_h);
            greeks["vega"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vega"];
            greeks["vomma"] = vega_vomma_extraction(prices[pair.first], p_plus_v[pair.first], p_minus_v[pair.first], v_h)["vomma"];
            for (auto& duo : greeks) { //for each greek:
                std::string contract_greek = pair.first + "_" + duo.first;
                results[contract_greek] = duo.second;
                std::cout << S << " Value for " << contract_greek << " is : " << results[contract_greek] << "\n";
            }
        }
        //reset original price and variable

        sigma = temp_sigma;
        S = temp_S;
        return results;
        }
};

class Simulation{ //Factory class

public:
    int seed = 42; //instance of random number

    //define default simmulation values 
    double S = 100, K = 105, T = 1.0, t = 0.0, sigma = 0.2;
    double r = 0.05, K2 = 5, T2 = 2.0, v_h = 0.01, s_h = 0.2 ;

    RandomNumber rand;


    int DayNumber = 252, SimulationNumber = 1000;

    explicit Simulation()
        : rand(seed, 1000, 252, S, r, sigma, T){} //constructor automatically creates instance rand w/ default values

    //int seed, int SimulationNumber, int DayNumber, double S, double r, double sigma, double T)
    RandomNumber CreateRandomNumber(int seed){
        RandomNumber newrand(seed, SimulationNumber, DayNumber, S, r, sigma, T);
        seed = seed + 1; //change seed for future instance of RandomNumber
        return newrand;
    }
    
    CompoundOption CreateCompoundOption(){
        CompoundOption compound(rand, S, K, K2, T, T2, t, sigma, r, v_h, s_h);
        return compound;
    }

    AsianOption CreateAsianOption(){
        AsianOption asian(rand, S, K, T, t, sigma, r, v_h, s_h);
        return asian;
    }

    EuropeanOption CreateEuropeanOption(){
        EuropeanOption euro(rand, S, K, T, t, sigma, r, v_h, s_h); //initiate an instance of the EuropeanOption object
        return euro;
    }
};


#endif