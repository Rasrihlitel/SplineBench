#pragma once
#ifndef SPLINEFUNC_H
#define SPLINEFUNC_H
#include <vector>

struct CoeffSpline {
	std::vector<double> a; //Коэффициент a
	std::vector<double> b; //Коэффициент b
	std::vector<double> c; //Коэффициент c
	std::vector<double> d; //Коэффициент d
};

void sweepMethod(const std::vector<double>& A, //Решение трехдиагональной матрицы
	const std::vector<double>& B,
	const std::vector<double>& C,
	const std::vector<double>& F,
	std::vector<double>& x);

CoeffSpline calculateCoeff(const std::vector<double>& x, const std::vector<double>& y); //Вычисление коэффициентов сплайна

double evaluateSpline(const CoeffSpline& coef, const std::vector<double>& x, double x_val); //Вычисление значения сплайна в точке

#endif //CUBICSPLINE_H