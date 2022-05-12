#include <iostream>
#include <vector>

/*   Только в фунцкии 
	 * n - строк матрицы
	 * a - над главной (нумеруется: [0;n-2])
	 * b - главная диагональ (нумеруется: [0;n-1])
	 * c - под главной (нумеруется: [1;n-1])
	 * d - столбец
	 * x - решение
*/

void solveTriagonalSlae(int n, std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& x) {

	double m;

	for (int i = 1; i < n; i++) {
		m = c[i] / b[i - 1];
		b[i] = b[i] - m * a[i - 1];
		d[i] = d[i] - m * d[i - 1];
	}

	x[n - 1] = d[n - 1] / b[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		x[i] = (d[i] - a[i] * x[i + 1]) / b[i];
	}
}


int main() {

	setlocale(LC_ALL, "Russian");

	std::cout << "Количество строк матрицы: " << std::endl;
	int n;
	std::cin >> n;
	
	double num;

	std::cout << "Верхняя диагональ: " << std::endl;
	std::vector<double> a;
	for (int i = 0; i < n-1; i++) {
		std::cin >> num;
		a.push_back(num);
	}
	
	std::cout << "Главная диагональ: " << std::endl;
	std::vector<double> b;
	for (int i = 0; i < n; i++) {
		std::cin >> num;
		b.push_back(num);
	}
	b[n - 1] = 0;
	std::cout << "Нижняя диагональ: " << std::endl;
	std::vector<double> c(n);
	for (int i = 1; i < n; i++) {
		std::cin >> num;
		c.push_back(num);
	}
	c[0] = 0;
	std::cout << "Столбец " << std::endl;
	std::vector<double> d;
	for (int i = 0; i < n; i++) {
		std::cin >> num;
		d.push_back(num);
	}
	
	std::vector<double>x(n);
	solveTriagonalSlae(n, a, b, c, d, x);
	
	std::cout << "Ответ " << std::endl;
	for (int i = 0; i < n; i++) {
		std::cout << x[i] << " ";
	}

	return 0;
}