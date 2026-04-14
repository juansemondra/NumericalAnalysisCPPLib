#include "menu.h"
#include "numericalanalysis.h"
#include <iostream>
#include <iomanip>
#include <limits>
#include <stdexcept>

void helper_function(){
    std::cout << "helper" << std::endl;
}

static void flush_cin()
{
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

template <typename T>
static T read_value(const std::string &prompt)
{
    T value;
    while (true)
    {
        std::cout << prompt;
        if (std::cin >> value) { flush_cin(); return value; }
        std::cout << "  Entrada inválida, intenta de nuevo.\n";
        std::cin.clear();
        flush_cin();
    }
}

static NumericalAnalysis::Function read_function()
{
    while (true)
    {
        std::cout << "Ingresa la función: ";
        std::string expression;
        std::getline(std::cin, expression);
        try
        {
            NumericalAnalysis::Function func;
            func.extract_expression(expression);
            return func;
        }
        catch (const std::invalid_argument &e)
        {
            std::cout << "  " << e.what() << "\nIntenta de nuevo.\n";
        }
    }
}

double call_bisection()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();

    double point_a = read_value<double>("Ingrese el punto a: ");
    double point_b = read_value<double>("Ingrese el punto b: ");

    if (func.evaluate(point_a) * func.evaluate(point_b) > 0)
    {
        std::cerr << "f(a) y f(b) deben tener signos opuestos para garantizar una raíz en [a, b].\n";
        return -1;
    }

    double tolerance = 0;
    while (tolerance <= 0)
    {
        tolerance = read_value<double>("Ingrese la tolerancia (> 0): ");
        if (tolerance <= 0) std::cout << "  La tolerancia debe ser un valor positivo.\n";
    }

    int iterations = 0;
    while (iterations <= 0)
    {
        iterations = read_value<int>("Ingrese el número de iteraciones (> 0): ");
        if (iterations <= 0) std::cout << "  El número de iteraciones debe ser positivo.\n";
    }

    return NumericalAnalysis::bisection(func, point_a, point_b, tolerance, iterations);
}

double call_fixed_point()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();

    double initial_point = read_value<double>("Ingrese el punto inicial: ");

    double tolerance = 0;
    while (tolerance <= 0)
    {
        tolerance = read_value<double>("Ingrese la tolerancia (> 0): ");
        if (tolerance <= 0) std::cout << "  La tolerancia debe ser un valor positivo.\n";
    }

    int iterations = 0;
    while (iterations <= 0)
    {
        iterations = read_value<int>("Ingrese el número de iteraciones (> 0): ");
        if (iterations <= 0) std::cout << "  El número de iteraciones debe ser positivo.\n";
    }

    return NumericalAnalysis::fixed_point(func, initial_point, tolerance, iterations);
}

double call_fake_position()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();

    double point_a = read_value<double>("Ingrese el punto a: ");
    double point_b = read_value<double>("Ingrese el punto b: ");

    if (func.evaluate(point_a) * func.evaluate(point_b) > 0)
    {
        std::cerr << "f(a) y f(b) deben tener signos opuestos para que el método de posición falsa funcione.\n";
        return -1;
    }

    double tolerance = 0;
    while (tolerance <= 0)
    {
        tolerance = read_value<double>("Ingrese la tolerancia (> 0): ");
        if (tolerance <= 0) std::cout << "  La tolerancia debe ser un valor positivo.\n";
    }

    int iterations = 0;
    while (iterations <= 0)
    {
        iterations = read_value<int>("Ingrese el número de iteraciones (> 0): ");
        if (iterations <= 0) std::cout << "  El número de iteraciones debe ser positivo.\n";
    }

    return NumericalAnalysis::fake_position(func, point_a, point_b, tolerance, iterations);
}

double call_newton_raphson()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();

    double initial_point = read_value<double>("Ingrese el punto inicial: ");

    double tolerance = 0;
    while (tolerance <= 0)
    {
        tolerance = read_value<double>("Ingrese la tolerancia (> 0): ");
        if (tolerance <= 0) std::cout << "  La tolerancia debe ser un valor positivo.\n";
    }

    int iterations = 0;
    while (iterations <= 0)
    {
        iterations = read_value<int>("Ingrese el número de iteraciones (> 0): ");
        if (iterations <= 0) std::cout << "  El número de iteraciones debe ser positivo.\n";
    }

    return NumericalAnalysis::newton_raphson(func, initial_point, tolerance, iterations);
}

double call_secant_method()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();

    double point_a = read_value<double>("Ingrese el punto a (x_{n-1}): ");
    double point_b = read_value<double>("Ingrese el punto b (x_n): ");

    while (point_a == point_b)
    {
        std::cout << "  Los dos puntos iniciales deben ser distintos.\n";
        point_b = read_value<double>("Ingrese el punto b (x_n): ");
    }

    double tolerance = 0;
    while (tolerance <= 0)
    {
        tolerance = read_value<double>("Ingrese la tolerancia (> 0): ");
        if (tolerance <= 0) std::cout << "  La tolerancia debe ser un valor positivo.\n";
    }

    int iterations = 0;
    while (iterations <= 0)
    {
        iterations = read_value<int>("Ingrese el número de iteraciones (> 0): ");
        if (iterations <= 0) std::cout << "  El número de iteraciones debe ser positivo.\n";
    }

    return NumericalAnalysis::secant_method(func, point_a, point_b, tolerance, iterations);
}

static std::string read_path(const std::string &prompt)
{
    std::string path;
    std::cout << prompt;
    std::getline(std::cin, path);
    return path;
}

static void print_solution(const NumericalAnalysis::Matrix& result)
{
    int n = result.getRows();
    std::cout << "\nVector solución:\n";
    for (int i = 0; i < n; i++)
        std::cout << "  x_" << i + 1 << " = "
                  << std::fixed << std::setprecision(6)
                  << result.get(i, 0) << "\n";
    std::cout << "\n";
}

static NumericalAnalysis::Matrix read_augmented_matrix()
{
    std::string filename = read_path("Ruta del archivo de la matriz: ");
    NumericalAnalysis::Matrix m(filename);
    int r = m.getRows();
    int c = m.getCols();

    if (r == 0)
    {
        std::cerr << "No se pudo leer la matriz.\n";
        return NumericalAnalysis::Matrix();
    }

    if (c == r + 1)
    {
        std::cout << "\nMatriz aumentada [A|b] (" << r << "x" << c << "):\n";
        m.print();
        return m;
    }

    if (c == r)
    {
        std::cout << "\nMatriz cuadrada A (" << r << "x" << c << "):\n";
        m.print();
        std::cout << "\nIngrese el vector b (" << r << " valores):\n";

        NumericalAnalysis::Matrix aug(r, r + 1);
        for (int i = 0; i < r; i++)
            for (int j = 0; j < r; j++)
                aug.set(i, j, m.get(i, j));

        for (int i = 0; i < r; i++)
            aug.set(i, r, read_value<double>(
                "  b_" + std::to_string(i + 1) + " = "));

        std::cout << "\nMatriz aumentada [A|b]:\n";
        aug.print();
        return aug;
    }

    std::cerr << "\nError: Se leyó una matriz de " << r << "x" << c
              << ". Se esperaba cuadrada (" << r << "x" << r
              << ") o aumentada (" << r << "x" << (r + 1) << ").\n\n";
    return NumericalAnalysis::Matrix();
}

void call_regressive_substitution()
{
    std::cin.ignore();
    NumericalAnalysis::Matrix matrix = read_augmented_matrix();
    if (matrix.getRows() == 0) return;

    NumericalAnalysis::Matrix result =
        NumericalAnalysis::regressive_substitution(matrix);
    print_solution(result);
}

void call_gaussian_elimination()
{
    std::cin.ignore();
    NumericalAnalysis::Matrix matrix = read_augmented_matrix();
    if (matrix.getRows() == 0) return;

    NumericalAnalysis::Matrix result =
        NumericalAnalysis::gaussian_elimination_with_regressive_substitution(matrix);
    print_solution(result);
}

void call_lu_substitution()
{
    std::cin.ignore();
    NumericalAnalysis::Matrix matrix = read_augmented_matrix();
    if (matrix.getRows() == 0) return;

    NumericalAnalysis::Matrix result =
        NumericalAnalysis::lu_substitution(matrix);
    print_solution(result);
}

void call_gauss_seidel()
{
    std::cin.ignore();
    NumericalAnalysis::Matrix matrix = read_augmented_matrix();
    if (matrix.getRows() == 0) return;

    int n = matrix.getRows();

    std::cout << "\nIngrese el vector inicial (" << n << " valores):\n";
    NumericalAnalysis::Matrix initial(n, 1);
    for (int i = 0; i < n; i++)
        initial.set(i, 0, read_value<double>(
            "  x0_" + std::to_string(i + 1) + " = "));

    double tolerance = 0;
    while (tolerance <= 0)
    {
        tolerance = read_value<double>("Ingrese la tolerancia (> 0): ");
        if (tolerance <= 0)
            std::cout << "  La tolerancia debe ser un valor positivo.\n";
    }

    int iterations = 0;
    while (iterations <= 0)
    {
        iterations = read_value<int>(
            "Ingrese el número máximo de iteraciones (> 0): ");
        if (iterations <= 0)
            std::cout << "  El número de iteraciones debe ser positivo.\n";
    }

    NumericalAnalysis::Matrix result =
        NumericalAnalysis::gauss_seidel(matrix, initial, tolerance, iterations);
    print_solution(result);
}

static void read_integration_params(double &a, double &b, int &n)
{
    a = read_value<double>("Ingrese el límite inferior a: ");
    b = read_value<double>("Ingrese el límite superior b: ");
    while (b <= a)
    {
        std::cout << "  b debe ser mayor que a.\n";
        b = read_value<double>("Ingrese el límite superior b: ");
    }
    n = 0;
    while (n <= 0)
    {
        n = read_value<int>("Ingrese el número de subintervalos n (> 0): ");
        if (n <= 0) std::cout << "  n debe ser un entero positivo.\n";
    }
}

void call_inferior_sums()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();
    func.print();

    double a, b;
    int n;
    read_integration_params(a, b, n);

    double result = NumericalAnalysis::inferior_sums(func, a, b, n);
    std::cout << "\nSuma inferior L(A) = " << std::fixed
              << std::setprecision(8) << result << "\n\n";
}

void call_superior_sums()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();
    func.print();

    double a, b;
    int n;
    read_integration_params(a, b, n);

    double result = NumericalAnalysis::superior_sums(func, a, b, n);
    std::cout << "\nSuma superior U(A) = " << std::fixed
              << std::setprecision(8) << result << "\n\n";
}

void call_trapezoidal_rule()
{
    std::cin.ignore();
    NumericalAnalysis::Function func = read_function();
    func.print();

    double a, b;
    int n;
    read_integration_params(a, b, n);

    double result = NumericalAnalysis::trapezoidal_rule(func, a, b, n);
    std::cout << "\nRegla del trapecio T(h) = " << std::fixed
              << std::setprecision(8) << result << "\n\n";
}

void print_menu(){
    std::cout << "\n========================================\n";
    std::cout << " Análisis Numérico\n";
    std::cout << "========================================\n";
    std::cout << "--- Primer Corte: Ceros de Funciones ---\n";
    std::cout << "  1. Método de la bisección\n";
    std::cout << "  2. Método del punto fijo\n";
    std::cout << "  3. Método de la posición falsa\n";
    std::cout << "  4. Método de Newton-Raphson\n";
    std::cout << "  5. Método de la secante\n";
    std::cout << "--- Segundo Corte: Sistemas Lineales ---\n";
    std::cout << "  6. Sustitución regresiva\n";
    std::cout << "  7. Eliminación gaussiana\n";
    std::cout << "  8. Factorización LU\n";
    std::cout << "  9. Método de Gauss-Seidel\n";
    std::cout << "--- Integración Numérica --------------\n";
    std::cout << " 10. Sumas inferiores\n";
    std::cout << " 11. Sumas superiores\n";
    std::cout << " 12. Regla del trapecio\n";
    std::cout << "----------------------------------------\n";
    std::cout << "  0. Salir\n";
    std::cout << "========================================\n";
    std::cout << "Opción: ";
}

void check_error(double value){
    if (value == -1) std::cout << "No se pudo encontrar resultado con la tolerancia propuesta." << std::endl;
    else std::cout << "El resultado de la operación es: " << value << std::endl;
}