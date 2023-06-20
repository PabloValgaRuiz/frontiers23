#include "model.h"
#include "ABC.h"
#include <pcg/pcg_random.hpp>

#include <random>
#include <fstream>
#include <iostream>
#include <tuple>
#include <string>
#include <future>
#include <mutex>
#include <math.h>

std::mutex mutex;

double readSigma(const std::string& path, const std::string& country_code){
    double sigma;
    std::ifstream file_sigma(path + "households.txt");
    std::string country_temp;
    std::string result;

    while(file_sigma.peek() != EOF){
        
        std::getline (file_sigma,country_temp, ',');
        std::getline (file_sigma, result, '\n');
        
        if(country_temp == country_code){
            sigma = std::stod(result);
            std::lock_guard<std::mutex> lock(mutex);
            std::cout << country_code << " sigma: " << sigma << std::endl;
            break;
        }
    }

    if(file_sigma.eof()){
        std::lock_guard<std::mutex> lock(mutex);
        std::cerr << "Couldn't find country " << country_code << std::endl;
    }
    return sigma;
}

double readPopulation(const std::string& path, const std::string& country_code){
    double pop;
    std::ifstream file_pop(path + "countries_population.txt");
    std::string country_temp;
    std::string result;

    while(file_pop.peek() != EOF){
        
        std::getline (file_pop,country_temp, ',');
        std::getline (file_pop, result, '\n');
        
        if(country_temp == country_code){
            pop = std::stod(result);
            std::lock_guard<std::mutex> lock(mutex);
            std::cout << country_code << " pop: " << pop << std::endl;
            break;
        }
    }

    if(file_pop.eof()){
        std::lock_guard<std::mutex> lock(mutex);
        std::cerr << "Couldn't find country " << country_code << std::endl;
    }

    return pop;
}

std::pair<double, double> readContacts(const std::string& path, const std::string& country){
    double k_active{}, k_passive{};

    std::ifstream file_countries(path + "countries_list.txt");
    std::string country_code;
    std::string country_name;
    while(file_countries.peek() != EOF){
        
        std::getline (file_countries,country_code, ',');
        std::getline (file_countries,country_name, '\n');
        
        if(country_code == country){
            break;
        }

    }

    std::ifstream file_active(path + "k_average/k_active.txt");
    std::string country_temp;
    std::string result;

    while(file_active.peek() != EOF){
        
        std::getline (file_active,country_temp, ',');
        std::getline (file_active, result, '\n');
        
        if(country_temp == country_name){
            k_active = std::stod(result);
            break;
        }

    }

    std::ifstream file_passive(path + "k_average/k_passive.txt");
    while(file_passive.peek() != EOF){
        
        std::getline (file_passive,country_temp, ',');
        std::getline (file_passive, result, '\n');
        
        if(country_temp == country_name){
            k_passive = std::stod(result);
            break;
        }
    }

    if(file_passive.eof()){
        std::lock_guard<std::mutex> lock(mutex);
        std::cerr << "Couldn't find country " << country << std::endl;
    }
    std::lock_guard<std::mutex> lock(mutex);
    std::cout << "k_active = " << k_active << ", k_passive = " << k_passive << std::endl;
    return {k_active, k_passive};
}

std::vector<double> readDead(const std::string& path) {
    std::ifstream file(path);

    if(!file) std::cerr << "Cannot open dead file" << std::endl;

    std::string date;
    double tempdead;

    std::vector<double> dead;

    while (file >> date >> tempdead) {
        dead.push_back(tempdead);
    }

    return dead;
}

std::vector<double> readMobility(const std::string& path) {
    std::ifstream file(path);

    if(!file) std::cerr << "Cannot open mobility file" << std::endl;

    std::string date;
    double tempMob;

    std::vector<double> Mob;

    while (file >> date >> tempMob) {
        Mob.push_back(tempMob);
    }
    return Mob;
}



 MarkovSEIRPDParameters buildMarkov(const FixedParamsABC& fixedParamsABC, const VariableParamsABC& variableParamsABC){
    return MarkovSEIRPDParameters{
        fixedParamsABC.p,
        fixedParamsABC["I_init"],
        (int)fixedParamsABC["N"],
        fixedParamsABC["k_active"],
        fixedParamsABC["k_passive"],
        fixedParamsABC["sigma"],
        fixedParamsABC["nu"],
        fixedParamsABC["mu"],

        variableParamsABC["beta"],
        variableParamsABC["permeability"],
        fixedParamsABC["IFR"],
        variableParamsABC["xi"],
        int(variableParamsABC["init_days"] + fixedParamsABC.p.size())
    };
}

std::vector<ResultsABC> ABC(int n_simulations, int n_top, PrioriParametersABC priori, const std::string& country) {

    pcg64 rng(pcg_extras::seed_seq_from<std::random_device>{});

    constexpr int chunk_size = 1e4;

    std::vector<ResultsABC> results; results.reserve(chunk_size);

    std::string dataPath{ "../../data/" };

    FixedParamsABC fixed{dataPath, country};

    auto D_obs = readDead(dataPath + "dead/dead_" + country + ".txt");

    for (int t = 0; t < n_simulations; t++) {

        VariableParamsABC variable(priori, rng); //Generate random variable set
        auto params = buildMarkov(fixed, variable);
        MarkovSEIRPD mkModel{params};

        auto D_sim = mkModel.iterate();
        if(variable["init_days"] >= 0)
            D_sim.erase(D_sim.begin(), D_sim.begin() + (int)floor(variable["init_days"]));
        else
            D_sim.insert(D_sim.begin(), -(int)floor(variable["init_days"]), 0); //Add zeros at the beginning (for the negative init_days)

        double value = 0;
        for (int i = priori["delay"].max(); i < D_obs.size(); i++) {
            value += log(std::abs(D_obs[i] - D_sim[i - (int)(variable["delay"])]) + 1);
        }

        results.push_back({ variable, D_sim, value });

        if (!(t % chunk_size) && t) {

            std::partial_sort(results.begin(), results.begin() + n_top, results.end(), [](const ResultsABC& a, const ResultsABC& b) {return a.value < b.value; });
            results.erase(results.begin() + n_top, results.end());
            std::lock_guard<std::mutex> lock(mutex);
            if(!((t * 100) % n_simulations))
                std::cout << 100 * (double)t / n_simulations << "%" << std::endl;
        }

    }

    std::partial_sort(results.begin(), results.begin() + n_top, results.end(), [](const ResultsABC& a, const ResultsABC& b) {return a.value < b.value; });
    results.erase(results.begin() + n_top, results.end());

    return results;
}


FixedParamsABC::FixedParamsABC(const std::string dataPath, const std::string country_code){
    int sigma = readSigma(dataPath, country_code);
    int pop = readPopulation(dataPath, country_code);
    auto ks = readContacts(dataPath, country_code);

    p = readMobility(dataPath + "mobility/mobility_" + country_code + ".txt");

    map["I_init"] = 1.0/pop;
    map["N"] = int(pop);

    map["k_active"] = ks.first;
    map["k_passive"] = ks.second;

    map["sigma"] = sigma;
    map["nu"] = 1.0 / 5.2;
    map["mu"] = 1.0 / 4.2;

    map["IFR"] = 0.01;
    // map["beta"] = 0.08;
    // map["xi"] = 0.10;
}

PrioriParametersABC::PrioriParametersABC(){
    map["beta"] = std::uniform_real_distribution<double>{0.01, 0.20};
    map["permeability"] = std::uniform_real_distribution<double>{0, 1};
    // map["IFR"] = std::uniform_real_distribution<double>{0.004, 0.024};
    map["xi"] = std::uniform_real_distribution<double>{1.0 / 20, 1.0 / 6};
    map["init_days"] = std::uniform_real_distribution<double>{0, 100};
    map["delay"] = std::uniform_real_distribution<double>{2, 20};
}