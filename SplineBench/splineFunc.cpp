#include "splineFunc.h"
#include <stdexcept>

void sweepMethod( //Решение трехдиагональной матрицы
	const std::vector<double>& A, 
	const std::vector<double>& B,
	const std::vector<double>& C,
	const std::vector<double>& F,
	std::vector<double>& x
) {
	int size = A.size();
	std::vector<double> alpha(size);
	std::vector<double> beta(size);
	alpha[0] = B[0] / C[0];
	beta[0] = F[0] / C[0];

	//Подсчет коэффициентов x[i] = beta[i] - alpha[i]x[i+1]
	for (int i = 1; i < size - 1; i++) {
		double denom = C[i] - alpha[i - 1] * A[i];
		alpha[i] = B[i] / denom;
		beta[i] = (F[i] - A[i] * beta[i - 1]) / denom;
	}

	//Подсчет резуьтирующего вектора - обратный ход
	x[size - 1] = (F[size - 1] - A[size - 1] * beta[size - 2]) / (C[size - 1] - A[size - 1] * alpha[size - 2]);
	for (int i = size - 2; i >= 0; i--) {
		x[i] = beta[i] - alpha[i] * x[i + 1];
	}
}

CoeffSpline calculateCoeff(const std::vector<double>& x, const std::vector<double>& y) { //Вычисление коэффициентов сплайна
	int size = x.size();

	//Проверка корректности входных данных
	if (size != y.size()) {
		throw std::invalid_argument("Размеры массивов x и y должны совпадать.");
	}
	if (size < 2) {
		throw std::invalid_argument("Нужно как минимум 2 точки.");
	}
	for (int i = 0; i < size - 1; i++) {
		if (x[i + 1] <= x[i]) {
			throw std::invalid_argument("Значения x должны быть строго возрастающими!");
		}
	}

	//Обозначение коэффициентов
	CoeffSpline coef;
	coef.a = y;
	coef.b.resize(size);
	coef.c.resize(size);
	coef.d.resize(size);
	std::vector<double> h(size - 1);

	//Подсчет h[i] = x[i] - x[i-1] (как в методичке, но в программе - смещение)
	for (int i = 0; i < size - 1; i++) {
		h[i] = x[i + 1] - x[i];
	}

	//Подсчет коэффициентов системы A[i]c[i-1] + C[i]c[i] + B[i]c[i+1] = F[i]
	std::vector<double> A(size, 0.0), C(size, 0.0), B(size, 0.0), F(size, 0.0);
	C[0] = 1.0;
	C[size - 1] = 1.0;
	for (int i = 1; i < size - 1; i++) {
		A[i] = h[i - 1];
		C[i] = 2.0 * (h[i - 1] + h[i]);
		B[i] = h[i];
		F[i] = 6.0 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1]);
	}

	//Решение системы
	sweepMethod(A, B, C, F, coef.c);

	//Вычисление коэффициентов
	for (int i = 1; i < size; i++) {
		coef.d[i] = (coef.c[i] - coef.c[i - 1]) / h[i - 1];
		coef.b[i] = h[i - 1] * coef.c[i] / 2.0 - h[i - 1] * h[i - 1] * coef.d[i] / 6.0 + (y[i] - y[i - 1]) / h[i - 1];
	}
	coef.d[0] = coef.d[1];
	coef.b[0] = coef.b[1];

	return coef;
}

double evaluateSpline(const CoeffSpline& coef, const std::vector<double>& x, double x_val) { //Вычисление значения сплайна в точке
	int size = x.size();
	if (x_val < x.front() || x_val > x.back()) {
		throw std::out_of_range("Точка лежит вне области интерполяции!");
	}
	//Бинарный поиск, так как быстрее
	int left = 0, right = size - 1;
	while (left < right - 1) {
		int mid = left + (right - left) / 2;
		if (x[mid] <= x_val) {
			left = mid;
		}
		else {
			right = mid;
		}
	}

	//Вычисление в точке
	double dx = x_val - x[right];
	return coef.a[right] + coef.b[right] * dx + (coef.c[right] / 2.0) * dx * dx + (coef.d[right] / 6.0) * dx * dx * dx;
}