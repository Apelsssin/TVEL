#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
#include <cmath>
#include <SFML/Graphics.hpp>
using namespace sf;
/*
Допущения:
1. Длина много больше диаметра. Теплопередачей по высоте пренебрегаем. 
2. Топливная таблетка однородна
3. Коэф. теплопроводности постоянны
*/
class TVEL {
	const double fuel_el_step = 0.01275; // Шаг размешения ТВЭЛов
	double S;
	double P;
	long double Pe;
	const double pi = 3.14159;
	const double delta_shell = 0.000685; //Толщина оболочки ТВЭЛа [м]
	const double delta_He = 0.000065; //Толщина газового зазора [м]
	const double outer_radius_fuel_tablet = 0.0038 ; //Наружний радиус топливной таблетки [м]
	const double inner_radius_fuel_tablet = 0.00115; //Внутренний радиус топливной таблетки [м]
	double coolant_temp_in_K = 571.35;
	double coolant_temp_out_K = 602.05;
	double coolant_temp_avg_K;
	double q_v = 1000000000; //Обьемное энерговыделение топлива [Вт/м3]
	std::vector <std::pair <double, double>> v_of_temp_K_r; //Вектор расперделения температур по радиусу в Цельсиях
	const double r_step = 0.0001; //Шаг построения графика по радиусу [м]
	double Q_1 = 0.00048066;// 88000 / 163 / 312 / 3600; // Расход через сечение 1 ТВЭЛ [м/с]
	double Cp = 4190; // Теплоносителя
	double coolant_velocity;
	double coolant_rho = 701.63;
	double coolant_lambda = 0.539; // [Вт/м*К]
	double fuel_lambda = 3.22;
	double He_lambda = 0.275; 
	double shell_lambda = 19.3;
	double K_r_max = 1;
	double K_z = 1;
	long double N_t;
	long double l = 0.0000001;
	size_t n_assembly = 1; //Число ТВС
	size_t n_fuel_el = 1; //Число ТВЭЛ
	double delta_t;

public:
	TVEL() { 
		N_t = 62923; //3200 * 1000000 / 163 / 312; //[Вт] на 1 ТВЭЛ
		S = sqrt(3) / 2 * fuel_el_step * fuel_el_step - pi * (outer_radius_fuel_tablet + delta_He + delta_shell) * (outer_radius_fuel_tablet + delta_He + delta_shell); // 1 ТВЭЛ
		P = 2 * pi * (outer_radius_fuel_tablet + delta_He + delta_shell);// 1 ТВЭЛ
		coolant_velocity = Q_1 / S;
		Pe = coolant_rho * Cp * coolant_velocity * l / coolant_lambda;
		delta_t = sqrt(1 / Pe) * l;
		coolant_temp_avg_K = (coolant_temp_in_K + coolant_temp_out_K) / 2;
		double rho_coolant_avg = 691.1;
		double mu_coolant_avg = 88.3 * 1000000;
		double nu_coolant_avg = mu_coolant_avg / rho_coolant_avg;
		double d_g_assembly = 4 * S / P;
		double coolant_alpha = nu_coolant_avg * coolant_lambda / d_g_assembly;
		
		//Все разности температур зависят только от z
		double temp_between_max_and_outer_radius_fuel_tablet = q_v * outer_radius_fuel_tablet * outer_radius_fuel_tablet / 4 / fuel_lambda * K_r_max * K_z;
		double temp_between_outer_radius_fuel_tablet_and_inner_shell = delta_He * N_t / He_lambda * n_assembly / ((outer_radius_fuel_tablet) * 2) / pi / n_fuel_el * K_z * K_r_max;
		double temp_between_inner_and_outer_shell = delta_shell * N_t * n_assembly / shell_lambda / ((outer_radius_fuel_tablet + delta_He) * 2) / pi / n_fuel_el * K_z * K_r_max;
		double temp_between_outer_shell_and_coolant = N_t * n_assembly * K_z * K_r_max / coolant_alpha / pi / 2 / (outer_radius_fuel_tablet + delta_He + delta_shell);
		
		double temp_outer_shell_K = coolant_temp_avg_K + temp_between_outer_shell_and_coolant;
		double temp_inner_shell_K = temp_outer_shell_K + temp_between_inner_and_outer_shell;
		double temp_outer_radius_fuel_tablet_K = temp_inner_shell_K + temp_between_outer_radius_fuel_tablet_and_inner_shell;
		double temp_max_K = temp_outer_radius_fuel_tablet_K + temp_between_max_and_outer_radius_fuel_tablet;

		auto push = [&](double r, double temp) {
			std::pair <double, double> var(r, temp);
			v_of_temp_K_r.push_back(var);
			};
		push(0, temp_max_K);
		push(outer_radius_fuel_tablet, temp_outer_radius_fuel_tablet_K);
		push(outer_radius_fuel_tablet + delta_He, temp_inner_shell_K);
		push(outer_radius_fuel_tablet + delta_He + delta_shell, temp_outer_shell_K);
		push(fuel_el_step / 2, coolant_temp_avg_K);

		//Температурное распределение оболочки
		double pred = outer_radius_fuel_tablet + delta_He + delta_shell;
		for (double i = outer_radius_fuel_tablet + delta_He + r_step; i < pred; i += r_step) {
			double temp_var = temp_inner_shell_K - (temp_between_inner_and_outer_shell) / log(1 + delta_shell / (outer_radius_fuel_tablet + delta_He)) * log(i / (outer_radius_fuel_tablet + delta_He));
			push(i, temp_var);
		}
		//Температурное распределение таблетки
		for (double i = 0; i < outer_radius_fuel_tablet; i += r_step) {
			double temp_var = temp_outer_radius_fuel_tablet_K + q_v / 4 / fuel_lambda * (outer_radius_fuel_tablet * outer_radius_fuel_tablet - i * i);
			push(i, temp_var);
		}

	 }
	void print() {
		std::cout << "Значение радиуса" << "\t" << "Значение температуры в К и С" << std::endl;
		for (auto it : v_of_temp_K_r) {
			std::cout << it.first << "\t\t" << it.second << "\t\t" << it.second - 273.15 << '\n';
		}
	};
	auto get_size() {
		return v_of_temp_K_r.size();
	}
	auto get_first(int i) {
		return v_of_temp_K_r[i].first;
	}
	auto get_second(int i) {
		return v_of_temp_K_r[i].second;
	}
	void sort() { //Сортировка по значению радиуса по возрастанию
		std::stable_sort(v_of_temp_K_r.begin(), v_of_temp_K_r.end(),
			[&](std::pair <double, double> pair1, std::pair <double, double> pair2) {
				return (pair1.first < pair2.first);
			});
	}
	TVEL(const TVEL&) = delete;
	void operator= (const TVEL&) = delete;
};

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Rus");
	TVEL t;
	t.sort();
	//t.print();

	srand(time(NULL));
	RenderWindow window(VideoMode(1500, 1000), L"Titul", Style::Default);
	int y0 = 900;
	int x0 = 100;
	VertexArray lines(LineStrip, t.get_size() + 6);
	lines[0].position = sf::Vector2f(x0 - 10, y0);
	lines[1].position = sf::Vector2f(x0 + 1500, y0);
	lines[2].position = sf::Vector2f(x0, y0);
	lines[3].position = sf::Vector2f(x0, y0 - 10); 
	lines[4].position = sf::Vector2f(x0, y0 + 10);
	lines[5].position = sf::Vector2f(x0, y0);
	for (long int i = 0; i < t.get_size() - 1; i += 1) {
		lines[i + 6].position = Vector2f(x0 + t.get_first(i) * 200000, y0 - t.get_second(i) / 5);
		lines[i + 7].position = Vector2f(x0 + t.get_first(i + 1) * 200000,y0 - t.get_second(i + 1) / 5);
	}

	while (window.isOpen())
	{
		Event event;
		while (window.pollEvent(event))
		{
			if (event.type == Event::Closed) window.close();
		}
		window.clear(Color::Transparent);
		window.draw(lines); 
		window.display();
	}
	return 0;
}