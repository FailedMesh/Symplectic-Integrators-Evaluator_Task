#include "Symplectic_Integration.h"
#include <iostream>
#include <vector>


void RungeKutta4thOrder(std::vector<std::vector<long double>>& accelerations,
	std::vector<std::vector<long double>>& velocities,
	std::vector<std::vector<long double>>& positions,
	std::vector<long double> E,
	std::vector<long double> B,
	long double time_step,
	long double q,
	long double m,
	long int max_time)
{

	int num_turns = 0;
	long int time = 1;
	std::vector<long double> k_vel_row(3, 0);
	std::vector<std::vector<long double>> k(4, k_vel_row);

	while (time < max_time)
	{
		update_kvector_vel(k[0], { 0, 0, 0 }, velocities[time - 1], E, B, 0, q, m, time_step);
		update_kvector_vel(k[1], k[0], velocities[time - 1], E, B, 0.5, q, m, time_step);
		update_kvector_vel(k[2], k[1], velocities[time - 1], E, B, 0.5, q, m, time_step);
		update_kvector_vel(k[3], k[2], velocities[time - 1], E, B, 1, q, m, time_step);

		velocities[time][0] = velocities[time - 1][0] + (k[0][0] / 6 + k[1][0] / 3 + k[2][0] / 3 + k[3][0] / 6);
		velocities[time][1] = velocities[time - 1][1] + (k[0][1] / 6 + k[1][1] / 3 + k[2][1] / 3 + k[3][1] / 6);
		velocities[time][2] = velocities[time - 1][2] + (k[0][2] / 6 + k[1][2] / 3 + k[2][2] / 3 + k[3][2] / 6);

		update_kvector_pos(k[0], { 0, 0, 0 }, velocities[time - 1], 0, time_step);
		update_kvector_pos(k[1], k[0], velocities[time - 1], 0.5, time_step);
		update_kvector_pos(k[2], k[1], velocities[time - 1], 0.5, time_step);
		update_kvector_pos(k[3], k[2], velocities[time - 1], 1, time_step);

		positions[time][0] = positions[time - 1][0] + (k[0][0] / 6 + k[1][0] / 3 + k[2][0] / 3 + k[3][0] / 6);
		positions[time][1] = positions[time - 1][1] + (k[0][1] / 6 + k[1][1] / 3 + k[2][1] / 3 + k[3][1] / 6);
		positions[time][2] = positions[time - 1][2] + (k[0][2] / 6 + k[1][2] / 3 + k[2][2] / 3 + k[3][2] / 6);

		system("CLS");
		std::cout << "Finished time step = " << time << "\nVelocity vector = ";
		print_vector(velocities[time]);
		std::cout << "\nPosition vector = ";
		print_vector(positions[time]);

		//Check if a turn is completed:
		if (positions[time][2] > 0 && positions[time - 1][2] < 0)
			num_turns++;

		time = time + 1;

	}

	std::cout << "\n" << num_turns;

}