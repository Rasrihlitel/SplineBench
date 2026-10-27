#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <omp.h>
#include "splineFunc.h"

int main() {
    std::cout << std::fixed << std::setprecision(6); 

    //Функция y = x^3 на отрезке [0, 3]
    std::vector<double> x = { 0.0, 1.0, 2.0, 3.0 };
    std::vector<double> y = { 0.0, 1.0, 8.0, 27.0 };

    std::cout << "Starting points: \n";
    for (int i = 0; i < x.size(); i++) {
        std::cout << "Point " << i << " (" << x[i] << ", " << y[i] << ")\n";
    }
    std::cout << "\n";

    //Построение сплайна
    CoeffSpline coef;
    try {
        coef = calculateCoeff(x, y);
        std::cout << "The spline is built\n\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error in calculation: " << e.what() << "\n";
        return 1;
    }

    //Вывод полученных коэффициентов
    std::cout << "Spline coefficients\n";
    for (size_t i = 0; i < x.size(); ++i) {
        std::cout << "Section " << i << ":\n";
        std::cout << " a = " << coef.a[i] << "\n";
        std::cout << " b = " << coef.b[i] << "\n";
        std::cout << " c = " << coef.c[i] << "\n";
        std::cout << " d = " << coef.d[i] << "\n";
    }
    std::cout << "\n";

    //Проверка в узлах интерполяции
    std::cout << "Check in the nodes\n";
    for (size_t i = 0; i < x.size(); ++i) {
        double val = evaluateSpline(coef, x, x[i]);
        std::cout << "S(" << x[i] << ") = " << val
            << " | Error rate: " << std::abs(val - y[i]) << "\n";
    }
    std::cout << "\n";

    //Проверяем в промежуточных точках
    std::cout << "Check at intermediate points\n";
    std::vector<double> test_points = { 0.5, 1.5, 2.5 };
    for (double pt : test_points) {
        double val = evaluateSpline(coef, x, pt);
        double exact_val = pt * pt * pt; 
        std::cout << "S(" << pt << ") = " << val
            << " | correct x^3 = " << exact_val
            << " | Error rate: " << std::abs(val - exact_val) << "\n";
    }
    std::cout << "\n";


    //Заполнение таблицы данными
    exportSplineToCSV(coef, x, 200);

    return 0;
}