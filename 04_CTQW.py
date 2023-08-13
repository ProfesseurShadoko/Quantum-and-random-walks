from QW.qw import *

####################################
# 3 DIMENSIONAL GRID OF SIZE NxNxN #
####################################
N=10
d=3

##################################
# INITIALIZATION OF THE PARTICLE #
##################################
particle:QParticle = tensor(
    QParticle(N,N//2),
    pow=d
)
particle = particle.gaussian(particle.mean(),2)
print(particle.std())
print(repr(particle))

#####################################
# INITIALIZATION OF THE HAMILTONIAN #
#####################################
# creating walk operator with coupling between direct neighbour positions
hamiltonian:Operator = Operator([N]*d) #equivalent to tensor(Operator(N),pow=3))
for position in hamiltonian.positions():
    for neighbour in particle.neighbours(*position):
        hamiltonian[position,neighbour] = 1

# defining the time step of the walk (adding 1j factor for resolving Schrödinger's equation)
dt = hamiltonian.timescale()* 1j * 0.01

# performing the walk (CTW stands for continous time walk)
walk = CTW(particle,hamiltonian,dt)
walk.solve(1000)
walk.run()
