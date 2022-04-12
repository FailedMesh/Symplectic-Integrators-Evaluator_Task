#include "Symplectic_Integration.h"

int main()
{
	std::vector<long double> Electric_Field{ 0, 0, 0 };
	std::vector<long double> Magnetic_Field{ 0, 0, 100000 };

	//Constants:
	long double c = 3 * pow(10, 8);
	long double m = 9.109 * pow(10, -31);
	long double q = -1.60218 * pow(10, -19);

	//Timestep:
	long double time_length = pow(10, -13);
	long double tStep = pow(10, -16);
	long int num_iterations = time_length / tStep;

	std::cout << "Number of timesteps = " << num_iterations;

	std::vector<long double> position_row(3, 0);
	std::vector<long double> velocity_row(3, 0);
	std::vector<long double> acceleration_row(3, 0);

	std::cout << "\nAllocating memory before starting integration...";

	std::vector<std::vector<long double>> positions(num_iterations, position_row);
	std::vector<std::vector<long double>> velocities(num_iterations, velocity_row);
	std::vector<std::vector<long double>> accelerations(num_iterations, acceleration_row);

	std::cout << "\nAllocation done";

	positions[0] = { 0, 0, 0 };
	velocities[0] = { 0, 0.9 * c, 0 };
	get_acceleration(accelerations[0], velocities[0], Electric_Field, Magnetic_Field, q, m);
	//accelerations[0] = curr_acc;

	//int max_turns = 2;

	RungeKutta4thOrder(accelerations, velocities, positions, Electric_Field, Magnetic_Field, tStep, q, m, num_iterations);

	return 0;
}