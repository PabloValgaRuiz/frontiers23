#include "model.h"
#include <math.h>
#include <iostream>
#define PROFILING 1
#include "benchmark.hpp"

MarkovSEIRPD::MarkovSEIRPD(const Params& _params) : params{_params}{
    if(_params.T > _params.p.size()){
        params.p.insert(params.p.begin(), params.T - params.p.size(), 1.0);
        params.p_residential.insert(params.p_residential.begin(), params.T - params.p_residential.size(), 1.0);
    }
    I = _params.I_init * _params.N;
    S_h = 0;
    S = _params.N - I;
    S_total = S;
    E = 0;
    Pd = 0;
    D = 0;
    R = 0;
}

std::vector<double> calculateCombinatorials(int n){
    std::vector<double> combinatorials(n);
    combinatorials[0] = 1;
    for(int i = 1; i < n; i++){
        combinatorials[i] = combinatorials[i-1] * (n - i + 1) / i;
    }
    return combinatorials;
}

std::vector<double> MarkovSEIRPD::iterate()
{
    PROFILE_FUNCTION();

    std::vector<double> daily_dead(params.T, 0);

    double Pactivo;
    double Ppasivo;
    double Pconfinado;
    double sh = 1;
    double Pcontagio;

    double k_active = params.k_active;
    double k_passive = params.k_passive;

    double sigma = round(params.sigma);

    auto probs_I_equals_i = calculateCombinatorials(params.sigma - 1);

    int t_init = 0;
    if(params.T < params.p.size()){
        t_init = params.p.size() - params.T;
    }

    for(int t = 0; t < params.T; t++){

        double pt = params.p[t + t_init] <= 1 ? params.p[t + t_init] : 1;
        double p_residential = params.p_residential[t + t_init];

        k_active = params.k_active * pt;
        k_passive = params.k_passive * p_residential;

        double rho = I/params.N;

        Pactivo = 1 - pow((1 - params.beta * rho), k_active);
        Ppasivo = 1 - pow((1 - params.beta * rho), k_passive);
        Pconfinado = 0;
        if (params.sigma > 2){
            for(int i = 1; i < params.sigma - 1; i++){
                Pconfinado += probs_I_equals_i[i] * pow(rho, i) * pow(1 - rho, params.sigma - 1 - i) * (1 - pow(1 - params.beta * i / (params.sigma - 1), k_passive));
            }
        }

        S_total = params.N - I - E - Pd - D - R;
        double S_active = S_total * pt;
        double S_inactive = S_total * (1 - pt) * params.permeability;
        double S_confined = S_total * (1 - pt) * (1 - params.permeability);
        
        D = params.xi * Pd + D;
        Pd = params.mu * params.IFR * I + (1 - params.xi) * Pd;
        R = params.mu * (1 - params.IFR) * I + R;
        I = params.nu * E + (1 - params.mu) * I;

        E = S_active * Pactivo + S_inactive * Ppasivo + S_confined * Pconfinado + (1 - params.nu) * E;

        daily_dead[t] = params.xi * Pd;
    }

    return daily_dead;
}
