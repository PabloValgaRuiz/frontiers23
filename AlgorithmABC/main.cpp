#include "model.h"
#include "ABC.h"
#include "ThreadPool.hpp"
#define PROFILING 1
#include "benchmark.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <filesystem>

#include "pcg/pcg_random.hpp"


/*This will use as many threads as the machine has available*/
static ThreadPool pool{std::thread::hardware_concurrency()};
/*To change the number of threads, comment the line above and uncomment the line below with the number of concurrent threads you want*/
//static ThreadPool pool{4};

/*The path where the results will be saved*/
constexpr std::string_view outPath{"../out/results/"};



template<typename Distribution, typename Type>
Distribution CompareResult(const Distribution& priori, Type result){
    Distribution new_priori = priori;
    if(priori.min() > result) new_priori = Distribution{result, priori.max()};
    if(priori.max() < result) new_priori = Distribution{priori.min(), result};
    return new_priori;
}

std::vector<ResultsABC> iterateABC(int max_chosen, const std::string& country){

    int n_top = 300;//When choosing the new priori distributions, take the 300 highest values

    PrioriParametersABC priori{};
    std::vector<ResultsABC> results;

    for(int i = 0; i < 2; i++){
        results = ABC(2e6, n_top, priori, country);

        for(auto& priori_param : priori.map)
            priori_param.second = std::uniform_real_distribution<double>(results[0].params[priori_param.first], results[0].params[priori_param.first]);

        for(int j = 0; j < results.size(); j++){
            for(auto& priori_param : priori.map)
                priori_param.second = CompareResult(priori_param.second, results[j].params[priori_param.first]);
        }

        for(const auto & param_name : {"init_days", "delay"}){
            priori[param_name] = std::uniform_real_distribution<double>((int)(priori[param_name].min()), (int)(priori[param_name].max())+1);
        }

        for(auto& priori_param : priori.map)
            std::cout << priori_param.first + "_priori: " << priori_param.second.min() << " - " << priori_param.second.max() << std::endl;
    }

    results = ABC(1e7, max_chosen, priori, country);

    return results;

}

int main() {
Instrumentor::Get().BeginSession("Session Name");

    std::filesystem::create_directory(outPath);

    int max_chosen = 1000;

    std::vector<std::future<void>> futures;

    const std::unordered_set<std::string> countries = {"AR", "AT", "BD", "BE", "BO", "BG", "CA", "CL", "CO", "EG", "FR", "GR", "DE", "GT",
                                                       "HN", "HU", "ID", "IQ", "IE", "IL", "IT", "KW", "LU", "MY", "MX", "MA", "NG", "PK",
                                                       "PA", "PH", "PL", "PT", "RO", "RU", "SA", "ZA", "ES", "CH", "TR", "US", "GB", "UA"};

    for(auto country : countries) futures.push_back(std::move(pool.enqueue([=]{
    
        auto begin = std::chrono::high_resolution_clock::now();

        auto results = iterateABC(max_chosen, country);

        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - begin).count() / 1000.0 << std::endl;
        
        std::ofstream fileValues(std::string{outPath} + "values_" + country + ".txt");
        std::ofstream fileResults(std::string{outPath} + "results_" + country + ".txt");
        std::ofstream fileDead(std::string{outPath} + "dead_" + country + ".txt");

        // fileValues << "beta,permeability,IFR,xi,init_days,delay,value\n";
        fileValues << "beta,permeability,xi,init_days,delay,value\n";

        fileResults << "iter,D_sim\n";
        fileDead << "mean,variance\n";

        for (int i = 0; i < max_chosen; i++) {
            const auto& result = results[i];

            // for(auto& param_name : {"beta","permeability","IFR","xi","init_days","delay"}){
            for(auto& param_name : {"beta","permeability","xi","init_days","delay"}){
                fileValues << result.params[param_name] << ",";
            }
            fileValues << result.value << "\n";

            for (int j = 0; j < result.D.size(); j++) {
                fileResults << i << "," << result.D[j] << "\n";
            }
        }
        
        int size = results[0].D.size();
        std::vector<double> mean(size);
        std::vector<double> mean2(size);

        for (int j = 0; j < size; j++) {
            for (int i = 0; i < max_chosen; i++) {
                mean[j] += results[i].D[j];
                mean2[j] += results[i].D[j] * results[i].D[j];
            }
            mean2[j] = (mean2[j] - mean[j] * mean[j] / max_chosen) / max_chosen;
            mean[j] /= max_chosen;

        }

        for (int i = 0; i < size; i++)
            fileDead << mean[i] << "," << mean2[i] << "\n";

    })));

    for(auto& future : futures) future.wait();
    futures.clear();

Instrumentor::Get().EndSession();
    return 0;
}