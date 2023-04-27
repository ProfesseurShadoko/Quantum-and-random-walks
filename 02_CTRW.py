from QW.qw import *

# creating a continous time random walk in 2 dimensions over a grid of size NxN
N=40

# initializing the particle on position (N//2,N//2)
particle:Particle = tensor(Particle(N,N//2),Particle(N,N//2))
print(repr(particle))
print(f"Magnitude on {N//2,N//2} :",particle[N//2,N//2])
print(f"Norm of the particle (should always be 1) :",particle.norm())

# creating walk operator with coupling between direct neighbour positions
adjacency_matrix:Operator = tensor(Operator(N),pow=2)
for position in adjacency_matrix.positions(): # all positions of the grid
    for neighbour in particle.neighbours(*position):
        adjacency_matrix[position,neighbour] = 1

laplacian = adjacency_matrix.laplacian()

# defining the time step of the walk
dt = adjacency_matrix.timescale() * 0.01

# performing the walk (CTW stands for continous time walk)
walk = CTW(particle,laplacian,dt)
walk.solve(500)
walk.run()

