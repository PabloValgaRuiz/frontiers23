#include "model.h"
#include <math.h>
#include <iostream>
#define PROFILING 1
#include "benchmark.hpp"

MarkovSEIRPD::MarkovSEIRPD(const Params& _params) : params{_params}{

    params.p.insert(params.p.begin(), params.T - params.p.size(), 1.0);
    
    I = _params.I_init * _params.N;
    S_h = 0;
    S = _params.N - I;
    S_total = S;
    E = 0;
    Pd = 0;
    D = 0;
    R = 0;
}

std::vector<double> MarkovSEIRPD::iterate()
{
    PROFILE_FUNCTION();

    std::vector<double> daily_dead(params.T);

    double Pactivo;
    double Pconfinado;
    double sh;
    double Pcontagio;

    for(int t = 0; t < params.T; t++){

        double pt = params.p[t] <= 1 ? params.p[t] : 1;

        double rho = I/params.N;
        if(rho < 0.0001)
        {
            Pactivo = params.k_active * params.beta * rho;
            Pconfinado = params.k_passive * params.beta * rho;
            sh = 1 - (params.sigma - 1) * rho;
        }
        else {
            Pactivo = 1 - pow((1 - params.beta * rho), params.k_active);
            Pconfinado = 1 - pow((1 - params.beta * rho),params.k_passive);
            sh = pow((1 - rho),(params.sigma - 1));
        }
        Pcontagio = pt * Pactivo + (1-pt)*(1 - sh * (1 - params.permeability)) * Pconfinado;
        
        S_total = params.N - I - E - Pd - D - R;
        
        S_h = S_total * (1 - pt) * sh * (1 - params.permeability);
        
        D = params.xi * Pd + D;
        Pd = params.mu * params.IFR * I + (1 - params.xi) * Pd;
        R = params.mu * (1 - params.IFR) * I + R;
        I = params.nu * E + (1 - params.mu) * I;
        E = S_total * Pcontagio + (1 - params.nu) * E;

        daily_dead[t] = params.xi * Pd;
    }

    return daily_dead;
}
