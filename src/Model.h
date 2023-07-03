#ifndef MODEL_H
#define MODEL_H

#include "../util/includes.h"

//templated alias of vec
template<typename T>
using vec = std::vector<T>;
using string = std::string;

template <typename T>
struct ModelData
{
    vec<vec<T>> values;  // Store the values of all variables at each repetition
};

class Model {
public:
    Model() = default;
    Model(const vec<string>& vnames, const vec<double>& rates,
        const vec<double>& inits, const vec<vec<double>>& tmat,
        const vec<std::function<double(const vec<double>&, int)>>& propensity,
        const string& debug);

    Model(const Model& other);

    ~Model();

    Model& operator=(const Model& o);

    std::tuple<vec<int>, vec<vec<vec<double>>>, int> getStats();
    std::tuple<double, int, vec<double>> run(const string& method = "SSA", int tmax = 10, int reps = 1);

private:
    constexpr static int pvSize = 1000;

    vec<string> vn;
    vec<double> rates;
    vec<double> inits;
    vec<std::function<double(const vec<double>&, int)>> pv;
    vec<vec<double>> tm;
    int nvars;
    vec<vec<vec<double>>> res;
    vec<int> time;
    vec<vec<vec<double>>> series;
    int steps;
    string debug;

    double pvFunc(int index, const vec<double>& ini);

    double randomDouble();

    int roulette(const vec<double>& pv, double limite);
};

#endif //MODEL_H
