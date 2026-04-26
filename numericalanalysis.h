#ifndef NUMERICALANALYSIS_H
#define NUMERICALANALYSIS_H

#include <map>
#include <string>
#include <vector>
#include <iostream>

namespace NumericalAnalysis {

    class Function {
    private:
        std::map            <std::string, double>               coeff;
        bool valid_key      (const std::string& key)            const;
        double evaluate_key (const std::string& key, double x)  const;

    public:
        Function();
        double  evaluate                (double x) const;
        double  derivate_evaluate       (double x) const;
        float   get                     (const std::string& key) const;
        void    update                  (const std::string& key, float val);
        void    add                     (const std::string& key, float val);
        void    extract_expression      (const std::string& expression);
        void    print                   () const;
    };

    class Matrix {
    private:
        std::vector<std::vector<double>> data;
        int rows;
        int columns;
    public:
        Matrix                          ();
        Matrix                          (const std::string& filename);
        Matrix                          (int rows, int columns);
        Matrix                          (const std::vector<std::vector<double>>& data);
        void    set                     (int row, int column, double value);
        double  get                     (int row, int column) const;
        void    print                   () const;
        void    add                     (const Matrix& other);
        void    subtract                (const Matrix& other);
        void    multiply                (const Matrix& other);
        void    divide                  (const Matrix& other);
        void    transpose               ();
        void    inverse                 ();
        double  determinant             ();
        int     rank                    ();
        void    read_from_file          (const std::string& filename);
        void    write_to_file           (const std::string& filename) const;
        int     getRows                 () const;
        int     getCols                 () const;
    };
    
    static double eval_arg          (const std::string &arg, double x);
    static double eval_arg_deriv    (const std::string &arg, double x);

    bool evaluate_tolerance (double xn, double xnp1, double tolerance);

    double bisection        (Function func, double point_a, double point_b, double tolerance, int iterations);
    double fixed_point      (Function func, double initial_point, double tolerance, int iterations);
    double fake_position    (Function func, double point_a, double point_b, double tolerance, int iterations);
    double newton_raphson   (Function func, double initial_point, double tolerance, int iterations);
    double secant_method    (Function func, double point_a, double point_b, double tolerance, int iterations);


    // Funciones segundo corte

    Matrix regressive_substitution(Matrix matrix);
    Matrix gaussian_elimination_step(Matrix matrix);
    Matrix gaussian_elimination_with_regressive_substitution(Matrix matrix);
    void   lu_factorization(Matrix A, Matrix& L, Matrix& U);
    Matrix lu_substitution(Matrix matrix);
    Matrix gauss_seidel(Matrix matrix, Matrix initial, double tolerance, int iterations);

    // Funciones segundo porte parte 2
    double inferior_sums(Function func, double a, double b, int n);
    double superior_sums(Function func, double a, double b, int n);
    double trapezoidal_rule(Function func, double a, double b, int n);
    double simpson_rule(Function func, double a, double b, int n);
}

#endif 