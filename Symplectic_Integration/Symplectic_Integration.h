#ifndef RUNGEKUTTA4THORDER_H

#include <iostream>
#include <vector>

long double norm(std::vector<long double>);

void print_vector(std::vector<long double> vect);

void get_acceleration(std::vector<long double>&,
	std::vector<long double>,
	std::vector<long double>,
	std::vector<long double>,
	long double,
	long double);

void update_kvector_vel(std::vector<long double>&,
	std::vector<long double>,
	std::vector<long double>,
	std::vector<long double>,
	std::vector<long double>,
	long double,
	long double,
	long double,
	long double);

void RungeKutta4thOrder(std::vector<std::vector<long double>>&,
	std::vector<std::vector<long double>>&,
	std::vector<std::vector<long double>>&,
	std::vector<long double>,
	std::vector<long double>,
	long double,
	long double,
	long double,
	long int);

void update_kvector_pos(std::vector<long double>&,
	std::vector<long double>,
	std::vector<long double>,
	long double,
	long double);


#endif // !RUNGEKUTTA4THORDER_H
