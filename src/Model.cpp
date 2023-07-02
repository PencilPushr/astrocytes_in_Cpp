#include "model.h"

#include <random>
#include <numeric>


Model::Model(const vec<string>& vnames, const vec<double>& rates,
    const vec<double>& inits, const vec<vec<double>>& tmat,
    const vec<std::function<double(const vec<double>&, int)>>& propensity,
    const string& debug)
    : vn(vnames), rates(rates), inits(inits), pv(propensity), tm(tmat),
    nvars(inits.size()), debug(debug) {}

//This feels braindead but I can't figure out a better solution
//if a better solution is found this can be remodeled, or better yet, simply don't have this at all.
Model::Model(const Model& other)
    : vn(other.vn), rates(other.rates), inits(other.inits), pv(other.pv), tm(other.tm),
    nvars(other.nvars), time(other.time), series(other.series), steps(other.steps),
    debug(other.debug)
{
    // Copy the contents of 'res' vector
    res.resize(other.res.size());
    for (std::size_t i = 0; i < other.res.size(); i++) {
        res[i].resize(other.res[i].size());
        for (std::size_t j = 0; j < other.res[i].size(); j++) {
            res[i][j].resize(other.res[i][j].size());
            for (std::size_t k = 0; k < other.res[i][j].size(); k++) {
                res[i][j][k] = other.res[i][j][k];
            }
        }
    }
}

std::tuple<vec<int>, vec<vec<vec<double>>>, int> Model::getStats() {
    return std::make_tuple(time, series, steps);
}

std::tuple<double, int, vec<double>> Model::run(const string& method, int tmax, int reps) {
    res.resize(tmax, vec<vec<double>>(nvars, vec<double>(reps, 0.0)));
    time.resize(tmax);
    series.resize(tmax, vec<vec<double>>(nvars, vec<double>(reps, 0.0)));

    for (int i = 0; i < tmax; i++) {
        time[i] = i;
    }

    if (method == "SSA") {
        vec<double> pv(pvSize, 0.0);

        for (int i = 0; i < reps; i++) {
            vec<double> ini = inits;
            vec<double> r = rates;
            int l = pv.size();

            double tc = 0.0;
            steps = 0;
            res[0][0][i] = ini[0];
            double a0 = 1.0;

            for (int j = 0; j < l; j++) {
                pv[j] = pvFunc(j, ini);
            }

            a0 = std::accumulate(pv.begin(), pv.end(), 0.0);
            double tau = (-1.0 / a0) * log(randomDouble());

            while (tc <= tmax) {
                double limite = a0 * randomDouble();
                int event = roulette(pv, limite);
                // Rest of the code...
                // Update tc, steps, res, etc.
            }
        }
    }
    else if (method == "SSAct") {
        // Implementation for "SSAct" method...
    }

    return std::make_tuple(tau, event, pv2);
}

double Model::pvFunc(int index, const vec<double>& ini) {
    return pv[index](ini, index);
}

double Model::randomDouble() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<double> dis(0.0, 1.0);
    return dis(gen);
}

int Model::roulette(const vec<double>& pv, double limite) {
    double a0 = std::accumulate(pv.begin(), pv.end(), 0.0);
    double aux = 0.0;
    int k = 0;
    while (k < pv.size() && aux < limite) {
        aux += pv[k];
        k++;
    }
    k--;
    return k;
}

Model& Model::operator=(const Model& other) {
    if (this == &other) {
        return *this;
    }

    vn = other.vn;
    rates = other.rates;
    inits = other.inits;
    pv = other.pv;
    tm = other.tm;
    nvars = other.nvars;
    res = other.res;
    time = other.time;
    series = other.series;
    steps = other.steps;
    debug = other.debug;

    return *this;
}