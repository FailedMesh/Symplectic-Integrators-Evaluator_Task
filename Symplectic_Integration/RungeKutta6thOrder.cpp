#include "Symplectic_Integration.h"


void RungeKutta6thOrder(std::vector<std::vector<long double>>& accelerations,
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
	std::vector<std::vector<long double>> k(7, k_vel_row);
	std::vector<long double> kupdate(3, 0);

	while (time < max_time)
	{
		update_kvector_vel(k[0], { 0, 0, 0 }, velocities[time - 1], E, B, 0, q, m, time_step);
		
		kupdate[0] = k[0][0];
		kupdate[1] = k[0][1];
		kupdate[2] = k[0][2];

		update_kvector_vel(k[1], kupdate, velocities[time - 1], E, B, 1, q, m, time_step);

		kupdate[0] = (3 * k[0][0] + k[1][0]) / 8;
		kupdate[1] = (3 * k[0][1] + k[1][1]) / 8;
		kupdate[2] = (3 * k[0][2] + k[1][2]) / 8;

		update_kvector_vel(k[2], kupdate, velocities[time - 1], E, B, 1, q, m, time_step);

		kupdate[0] = (8 * k[0][0] + 2 * k[1][0] + 8 * k[2][0]) / 27;
		kupdate[1] = (8 * k[0][1] + 2 * k[1][1] + 8 * k[2][1]) / 27;
		kupdate[2] = (8 * k[0][2] + 2 * k[1][2] + 8 * k[2][2]) / 27;

		update_kvector_vel(k[3], kupdate, velocities[time - 1], E, B, 1, q, m, time_step);

		kupdate[0] = (3 * (3 * sqrt(21) - 7) * k[0][0] - 8 * (7 - sqrt(21) * k[1][0] + 48 * (7 - sqrt(21)) * k[2][0] - 3 * (21 - sqrt(21)) * k[3][0])) / 392;
		kupdate[1] = (3 * (3 * sqrt(21) - 7) * k[0][1] - 8 * (7 - sqrt(21) * k[1][1] + 48 * (7 - sqrt(21)) * k[2][1] - 3 * (21 - sqrt(21)) * k[3][1])) / 392;
		kupdate[2] = (3 * (3 * sqrt(21) - 7) * k[0][2] - 8 * (7 - sqrt(21) * k[1][2] + 48 * (7 - sqrt(21)) * k[2][2] - 3 * (21 - sqrt(21)) * k[3][2])) / 392;

		update_kvector_vel(k[4], kupdate, velocities[time - 1], E, B, 1, q, m, time_step);

		kupdate[0] = (-5 * (231 + 51 * sqrt(21)) * k[0][0] - 40 * (7 + sqrt(21)) * k[1][0] - 320 * sqrt(21) * k[2][0] + 3 * (21 + 121 * sqrt(21)) * k[3][0] + 392 * (6 + sqrt(21)) * k[4][0]) / 1960;
		kupdate[1] = (-5 * (231 + 51 * sqrt(21)) * k[0][1] - 40 * (7 + sqrt(21)) * k[1][1] - 320 * sqrt(21) * k[2][1] + 3 * (21 + 121 * sqrt(21)) * k[3][1] + 392 * (6 + sqrt(21)) * k[4][1]) / 1960;
		kupdate[2] = (-5 * (231 + 51 * sqrt(21)) * k[0][2] - 40 * (7 + sqrt(21)) * k[1][2] - 320 * sqrt(21) * k[2][2] + 3 * (21 + 121 * sqrt(21)) * k[3][2] + 392 * (6 + sqrt(21)) * k[4][2]) / 1960;

		update_kvector_vel(k[5], kupdate, velocities[time - 1], E, B, 1, q, m, time_step);

		kupdate[0] = (15 * (22 + 7 * sqrt(21)) * k[0][0] + 120 * k[1][0] + 40 * (7 * sqrt(21) - 5) * k[2][0] - 63 * (3 * sqrt(21) - 2) * k[3][0] - 14 * (49 + 9 * sqrt(21)) * k[4][0] + 70 * (7 - sqrt(21)) * k[5][0]) / 180;
		kupdate[1] = (15 * (22 + 7 * sqrt(21)) * k[0][1] + 120 * k[1][1] + 40 * (7 * sqrt(21) - 5) * k[2][1] - 63 * (3 * sqrt(21) - 2) * k[3][1] - 14 * (49 + 9 * sqrt(21)) * k[4][1] + 70 * (7 - sqrt(21)) * k[5][1]) / 180;
		kupdate[2] = (15 * (22 + 7 * sqrt(21)) * k[0][2] + 120 * k[1][2] + 40 * (7 * sqrt(21) - 5) * k[2][2] - 63 * (3 * sqrt(21) - 2) * k[3][2] - 14 * (49 + 9 * sqrt(21)) * k[4][2] + 70 * (7 - sqrt(21)) * k[5][2]) / 180;

		update_kvector_vel(k[5], kupdate, velocities[time - 1], E, B, 1, q, m, time_step);

		velocities[time][0] = velocities[time - 1][0] + (9 * k[0][0] + 64 * k[2][0] / 3 + 49 * k[4][0] / 3 + 49 * k[5][0] + 9 * k[6][0]) / 180;
		velocities[time][1] = velocities[time - 1][1] + (9 * k[0][1] + 64 * k[2][1] / 3 + 49 * k[4][1] / 3 + 49 * k[5][1] + 9 * k[6][1]) / 180;
		velocities[time][2] = velocities[time - 1][2] + (9 * k[0][2] + 64 * k[2][2] / 3 + 49 * k[4][2] / 3 + 49 * k[5][2] + 9 * k[6][2]) / 180;

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