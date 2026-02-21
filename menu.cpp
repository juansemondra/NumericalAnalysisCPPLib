#include "menu.h"
#include "numericalanalysis.h"
#include <iostream>
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

void print_menu(){
    
    std::cout << std::endl << "Bienvendo al sistema de análisis numérico." << std::endl;
    std::cout << "Opciones: " << std::endl;
    std::cout << "1. Método de la bisección " << std::endl;
    std::cout << "2. Método del punto fijo " << std::endl;
    std::cout << "3. Método de la posición falsa " << std::endl;
    std::cout << "4. Método de la Newton-Raphson " << std::endl;
    std::cout << "5. Método de la secante " << std::endl;
    std::cout << "0. Salir " << std::endl;
}

void check_error(double value){
    if (value == -1) std::cout << "No se pudo encontrar resultado con la tolerancia propuesta." << std::endl;
    else std::cout << "El resultado de la operación es: " << value << std::endl;
}