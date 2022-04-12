#include "Symplectic_Integration.h"
#include <iostream>
#include <vector>

void print_vector(std::vector<long double> vect)
{
	for (int i = 0; i < vect.size(); i++)
		std::cout << vect[i] << " ";
}

void get_acceleration(std::vector<long double>& acc,
	std::vector<long double> vel,
	std::vector<long double> E,
	std::vector<long double> B,
	long double q,
	long double m)
{

	acc[0] = q * (E[0] + (vel[1] * B[2] - vel[2] * B[1])) / m;
	acc[1] = q * (E[1] + (vel[2] * B[0] - vel[0] * B[2])) / m;
	acc[2] = q * (E[2] + (vel[0] * B[1] - vel[1] * B[0])) / m;

}

void update_kvector_vel(std::vector<long double>& kvector,
	std::vector<long double> kupdate,
	std::vector<long double> vel,
	std::vector<long double> E,
	std::vector<long double> B,
	long double weight,
	long double q,
	long double m,
	long double h)
{

	kvector[0] = h * q * (E[0] + ((vel[1] + kupdate[1] * weight) * B[2] - (vel[2] + kupdate[2] * weight) * B[1])) / m;
	kvector[1] = h * q * (E[1] + ((vel[2] + kupdate[2] * weight) * B[0] - (vel[0] + kupdate[0] * weight) * B[2])) / m;
	kvector[2] = h * q * (E[2] + ((vel[0] + kupdate[0] * weight) * B[1] - (vel[1] + kupdate[1] * weight) * B[0])) / m;

}

void update_kvector_pos(std::vector<long double>& kvector,
	std::vector<long double> kupdate,
	std::vector<long double> vel,
	long double weight,
	long double h)
{

	kvector[0] = h * (vel[0] + kupdate[0]);
	kvector[1] = h * (vel[1] + kupdate[1]);
	kvector[2] = h * (vel[2] + kupdate[2]);

}