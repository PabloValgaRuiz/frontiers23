#pragma once

#include <vector>
#include <random>
#include <pcg/pcg_random.hpp>
#include <unordered_map>

struct FixedParamsABC{
    std::vector<double> p{};

    std::unordered_map<std::string,double> map;

    FixedParamsABC(const std::string dataPath, const std::string country_code);



    inline double& operator[](const std::string& param){
        return map.at(param);
    }
    inline double operator[](const std::string& param) const{
        return map.at(param);
    }
};


struct PrioriParametersABC{

    std::unordered_map<std::string, std::uniform_real_distribution<double>> map;

    PrioriParametersABC();

    inline std::uniform_real_distribution<double>& operator[](const std::string& param){
        return map.at(param);
    }

};

struct VariableParamsABC{
    std::unordered_map<std::string,double> map;

    VariableParamsABC(PrioriParametersABC priori, pcg64& rng){
        for(auto& param : priori.map)
            map[param.first] = param.second(rng);
    };

    inline double& operator[](const std::string& param){
        return map.at(param);
    }
    inline double operator[](const std::string& param) const{
        return map.at(param);
    }
};

struct ResultsABC {
    VariableParamsABC params;
    std::vector<double> D;
    double value;
};


std::vector<ResultsABC> ABC(int n_simulations, int n_top, PrioriParametersABC priori, const std::string& country);