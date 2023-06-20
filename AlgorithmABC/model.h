#pragma once

#include <vector>


struct MarkovSEIRPDParameters{

    std::vector<double> p{};
    double I_init{};
    int N{};
    double k_active{}, k_passive{}, sigma{}, nu{}, mu{};


    double beta{}, permeability{}, IFR{}, xi{};
    int T{};
};

class MarkovSEIRPD{

private:
    using Params = MarkovSEIRPDParameters;
    
    Params params{};

    double I{}, S_h{}, S{}, S_total{}, E{}, Pd{}, D{}, R{};

public:

    MarkovSEIRPD(const Params& _params);

    std::vector<double> iterate();

};