// reactions.cpp
#include "reactions.h"

Reactions::Reactions(const std::vector<double>& ini, int number_of_rec, int number_of_var)
    : ini(ini), number_of_rec(number_of_rec), number_of_var(number_of_var) 
{}

double Reactions::calculate_sigma0(int i) 
{
    // ini[(i/numberofrec)*numberofvar+0]*ini[(i/numberofrec)*numberofvar+1]
    return rec(i, FormulaVariables::v1) * rec(i, FormulaVariables::Y);
}

double Reactions::calculate_sigma1(int i) 
{
    return rec(i, FormulaVariables::vin);
}

double Reactions::calculate_sigma2(int i) 
{
    using namespace FormulaVariables;

    double VM3_val = rec(i, VM3);
    double n_val = rec(i, n);
    double Ca_val = rec(i, Ca);
    double kcaaa_val = rec(i, kcaaa);
    double kcai_val = rec(i, kcai);
    double Ia_val = rec(i, Ia);
    double m_val = rec(i, m);
    double kip3_val = rec(i, kip3);
    double S_val = rec(i, S);

    // Calculate power terms
    double Ca_n = std::pow(Ca_val, n_val);
    double kcaaa_n = std::pow(kcaaa_val, n_val);
    double kcai_n = std::pow(kcai_val, n_val);
    double Ia_m = std::pow(Ia_val, m_val);
    double kip3_m = std::pow(kip3_val, m_val);

    // Calculate the complex formula parts
    double term1 = 4 * VM3_val * (kcaaa_n * Ca_n);
    double term2 = (Ca_n + kcaaa_n) * (Ca_n + kcai_n);
    double term3 = Ia_m / (Ia_m + kip3_m);

    return (term1 / term2) * term3 * (S_val - Ca_val);
}

double Reactions::calculate_sigma3(int i) 
{
    using namespace FormulaVariables;

    double VM2_val = rec(i, VM2);
    double Ca_val = rec(i, Ca);
    double K2_val = rec(i, K2);

    double Ca_squared = Ca_val * Ca_val;
    double K2_squared = K2_val * K2_val;

    return VM2_val * (Ca_squared / (Ca_squared + K2_squared));
}

double Reactions::calculate_sigma4(int i) {
    using namespace FormulaVariables;
    return rec(i, kf) * rec(i, S);
}

double Reactions::calculate_sigma5(int i) 
{
    using namespace FormulaVariables;
    return rec(i, kf) * rec(i, Ca);
}

double Reactions::calculate_sigma6(int i) 
{
    using namespace FormulaVariables;
    return rec(i, kout) * rec(i, Ca);
}

double Reactions::calculate_sigma7(int i) 
{
    using namespace FormulaVariables;

    double vp_val = rec(i, vp);
    double vin_val = rec(i, vin);
    double Ca_val = rec(i, Ca);
    double kp_val = rec(i, kp);

    double vin_squared = vin_val * vin_val;
    double Ca_squared = Ca_val * Ca_val;
    double kp_squared = kp_val * kp_val;

    return vp_val * (vin_squared / (Ca_squared + kp_squared));
}

double Reactions::calculate_sigma8(int i) 
{
    using namespace FormulaVariables;
    return rec(i, kdeg) * rec(i, Ia);
}

double Reactions::calculate_sigma9(int i) 
{
    using namespace FormulaVariables;

    double ka_val = rec(i, ka);
    double Ca_val = rec(i, Ca);
    double q_val = rec(i, q);
    double K_val = rec(i, K);
    double W_val = rec(i, W);
    double A_val = rec(i, A);

    double Ca_q = std::pow(Ca_val, q_val);
    double K_q = std::pow(K_val, q_val);

    return ka_val * (Ca_q / (K_q + Ca_q)) * (W_val - 0 * A_val);
}

double Reactions::calculate_base_diffusion(int i, int target_idx) 
{
    using namespace FormulaVariables;

    double D_val = rec(i, D);
    double l_val = rec(i, l);
    double Ca_current = rec(i, Ca);
    double Ca_target = ini[target_idx * number_of_var + Ca];

    if (Ca_current <= Ca_target) {
        return 0.0;
    }

    return (D_val / l_val) * std::fabs(Ca_current - Ca_target);
}

double Reactions::calculate_diffusion_x_positive(int i) 
{
    int target_idx = i + 1;
    if (target_idx >= ini.size() / number_of_var) return 0.0;

    double base = calculate_base_diffusion(i, target_idx);
    return base * rec(i, FormulaVariables::phh);
}

double Reactions::calculate_diffusion_x_negative(int i) 
{
    int target_idx = i - 1;
    if (target_idx < 0) return 0.0;

    double base = calculate_base_diffusion(i, target_idx);
    return base * rec(i, FormulaVariables::phh);
}

double Reactions::calculate_diffusion_y_positive(int i) 
{
    // Need grid width for y direction calculation
    int grid_width = std::sqrt(ini.size() / number_of_var);
    int target_idx = i + grid_width;
    if (target_idx >= ini.size() / number_of_var) return 0.0;

    double base = calculate_base_diffusion(i, target_idx);
    return base * rec(i, FormulaVariables::phh);
}

double Reactions::calculate_diffusion_y_negative(int i) 
{
    int grid_width = std::sqrt(ini.size() / number_of_var);
    int target_idx = i - grid_width;
    if (target_idx < 0) return 0.0;

    double base = calculate_base_diffusion(i, target_idx);
    return base * rec(i, FormulaVariables::phh);
}

double Reactions::calculate_diffusion_z_positive(int i) 
{
    int grid_area = ini.size() / (number_of_var * 3); // Assuming 3 layers
    int target_idx = i + grid_area;
    if (target_idx >= ini.size() / number_of_var) return 0.0;

    double base = calculate_base_diffusion(i, target_idx);
    return base * rec(i, FormulaVariables::phh);
}

double Reactions::calculate_diffusion_z_negative(int i) 
{
    int grid_area = ini.size() / (number_of_var * 3); // Assuming 3 layers
    int target_idx = i - grid_area;
    if (target_idx < 0) return 0.0;

    double base = calculate_base_diffusion(i, target_idx);
    return base * rec(i, FormulaVariables::phh);
}

double Reactions::rec(int i, int x) 
{
    return ini[(i / number_of_rec) * number_of_var + x];
}

std::vector<double> Reactions::get_all_propensities(int i) 
{
    std::vector<double> props;
    props.reserve(number_of_rec);

    // Add all reaction propensities
    props.push_back(calculate_sigma0(i));
    props.push_back(calculate_sigma1(i));
    props.push_back(calculate_sigma2(i));
    props.push_back(calculate_sigma3(i));
    props.push_back(calculate_sigma4(i));
    props.push_back(calculate_sigma5(i));
    props.push_back(calculate_sigma6(i));
    props.push_back(calculate_sigma7(i));
    props.push_back(calculate_sigma8(i));
    props.push_back(calculate_sigma9(i));

    // Add diffusion reactions
    props.push_back(calculate_diffusion_x_positive(i));
    props.push_back(calculate_diffusion_x_negative(i));
    props.push_back(calculate_diffusion_y_positive(i));
    props.push_back(calculate_diffusion_y_negative(i));
    props.push_back(calculate_diffusion_z_positive(i));
    props.push_back(calculate_diffusion_z_negative(i));

    return props;
}