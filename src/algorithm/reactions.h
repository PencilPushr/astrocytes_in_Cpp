// reactions.h
#pragma once

#include <vector>
#include <cmath>
#include <array>

namespace FormulaVariables {
    enum Variables {
        v1 = 0,     // 0v1
        Y = 1,      // 1Y
        vin = 2,    // 2vin
        VM2 = 3,    // 3VM2
        Ca = 4,     // 4C
        n = 5,      // 5n
        K2 = 6,     // 6K2
        VM3 = 7,    // 7VM3
        kout = 8,   // 8kout
        S = 9,      // 9S
        kf = 10,    // 10kf
        kp = 11,    // 11kp
        kdeg = 12,  // 12kdeg
        vp = 13,    // 13vp
        kcaaa = 14, // 14kcaaa
        kcai = 15,  // 15kcai
        kip3 = 16,  // 16kip3
        Ia = 17,    // 17Z
        q = 18,     // 18q
        W = 19,     // 19W
        A = 20,     // 20A
        kia = 21,   // 21kia
        D = 22,     // 22D
        l = 23,     // 23l
        K = 24,     // 24K
        ka = 25,    // 25ka
        m = 26,     // 26m
        phh = 27,   // 27phh
        phl = 28,   // 28phl
        plh = 29    // 29plh
    };
};

class Reactions {
public:
    Reactions(const std::vector<double>& ini, int number_of_rec, int number_of_var);

    // Core reaction calculations
    double calculate_sigma0(int i);  // Reaction 0
    double calculate_sigma1(int i);  // Reaction 1
    double calculate_sigma2(int i);  // Reaction 2
    double calculate_sigma3(int i);  // Reaction 3
    double calculate_sigma4(int i);  // Reaction 4
    double calculate_sigma5(int i);  // Reaction 5
    double calculate_sigma6(int i);  // Reaction 6
    double calculate_sigma7(int i);  // Reaction 7
    double calculate_sigma8(int i);  // Reaction 8
    double calculate_sigma9(int i);  // Reaction 9 (Receptor reactions)

    // Diffusion reactions
    double calculate_diffusion_x_positive(int i);  // C_x,y,z -> C_x+1,y,z
    double calculate_diffusion_x_negative(int i);  // C_x,y,z -> C_x-1,y,z
    double calculate_diffusion_y_positive(int i);  // C_x,y,z -> C_x,y+1,z
    double calculate_diffusion_y_negative(int i);  // C_x,y,z -> C_x,y-1,z
    double calculate_diffusion_z_positive(int i);  // C_x,y,z -> C_x,y,z+1
    double calculate_diffusion_z_negative(int i);  // C_x,y,z -> C_x,y,z-1

    // Helper functions
    std::vector<double> get_all_propensities(int i);

private:
    double rec(int i, int x);
    double calculate_base_diffusion(int i, int target_idx);

private:
    const std::vector<double>& ini;
    int number_of_rec;
    int number_of_var;
    static constexpr double ALPHA = 0.01;
};