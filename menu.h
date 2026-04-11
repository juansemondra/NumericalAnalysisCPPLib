#ifndef MENU_H
#define MENU_H

    void helper_function();
    double call_bisection();
    double call_fixed_point();
    double call_fake_position();
    double call_newton_raphson();
    double call_secant_method();

    void call_regressive_substitution();
    void call_gaussian_elimination();
    void call_lu_substitution();
    void call_gauss_seidel();

    void call_inferior_sums();
    void call_superior_sums();
    void call_trapezoidal_rule();

    void print_menu();
    void check_error(double value);

#endif