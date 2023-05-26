from QW.qw import *

# creating a continous time random walk in 3 dimensions over a grid of size NxNxN
N=11

particle:QParticle = tensor(
    QParticle(N),
    pow=3
)
particle = particle.uniform()
print(repr(particle))

# creating walk operator with coupling between direct neighbour positions
hamiltonian:Operator = Operator([N,N,N]) #equivalent to tensor(Operator(N),pow=3))
for position in hamiltonian.positions():
    for neighbour in particle.neighbours(*position):
        hamiltonian[position,neighbour] = 1

# adding a target to the walk
target = (5,7,3)
for neighbour in particle.neighbours(*target):
    hamiltonian[target,neighbour] = 2
    hamiltonian[neighbour,target] = 2

# defining the time step of the walk (adding 1j factor for resolving Schrödinger's equation)
dt = hamiltonian.timescale()* 1j * 0.01

# performing the walk (CTW stands for continous time walk)
walk = CTW(particle,hamiltonian,dt)
walk.solve(3000)
walk.run()

walk.follow_target(target)