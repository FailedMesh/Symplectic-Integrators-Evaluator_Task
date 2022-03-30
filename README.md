# Symplectic-Integrators-Evaluator_Task
## An implementation of the Symplectic Integration Task for GSOC-2022

The Force equation of an electron is integrated by using various symplectic integration methods to generate its trajectory, radius and momentum scaling factor with
each turn. The motion of the electron can be controlled and visualized by changing the values of the electric and magnetic fields.

### Methods implemented in this repository:
1) Dormand Prince 5th Order
2) Runge Kuta 4th Order
3) Runge Kutta 6th Order
4) Leapfrog Method

Each of these methods generated the same trajectory for the electron. But because of the difference in integration method, the momentum scaling and time steps required for a specific number of turns tend to change. The value of time step and maximum number of turns can be changed in the code to see how they affect the momentum scaling, keeping in mind the time complexity of the code

Few examples of graphs generated:

![Leapfrog Trajectory](https://user-images.githubusercontent.com/68490344/160836058-6bff0a77-467c-4f9c-b2cf-aefe30dbf73f.png)
![Leapfrog Momentum](https://user-images.githubusercontent.com/68490344/160836092-337a9bc9-627f-4fd2-a2fd-cc7bcf9cc21c.png)

Graphs are generated automatically by running the codes.


