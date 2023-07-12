from QW.qw import *

#################################
# 2 DIMENSIONAL GRID OF SIZE N² #
#################################
N=40


########################################################
# INITIALIZATION OF THE PARTICLE ON POSITION (N/2,N/2) #
########################################################

particle:Particle = tensor(Particle(N,N//2),Particle(N,N//2))
print(repr(particle))
print(f"Magnitude on {N//2,N//2} :",particle[N//2,N//2])
print(f"Norm of the particle (should always be 1) :",particle.norm())


#######################################
# INITIALIZATION OF THE WALK OPERATOR #
#######################################
adjacency_matrix:Operator = tensor(Operator(N),pow=2)
for position in adjacency_matrix.positions(): # all positions of the grid
    for neighbour in particle.neighbours(*position):
        adjacency_matrix[position,neighbour] = 1

laplacian = adjacency_matrix.laplacian()

# defining the time step of the walk
dt = adjacency_matrix.timescale() * 0.01


############################################
# PERFORMING THE CONTINOUS-TIME-WALK (CTW) #
############################################
walk = CTW(particle,laplacian,dt)
walk.solve(500)
walk.run()

