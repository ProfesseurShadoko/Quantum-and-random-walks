from QW.qw import *

d=2
N=11
std=4
target=tuple(np.random.randint(N//3,2*N//3) for _ in range(d))

gammas=[0,1.2,1.4,2.0] #gamma-1D = gammas[1] = 1.2


# initializing particle with gaussian distribution around the center of the grid
particle:QParticle = tensor(QParticle(N),pow=d)
particle = particle.gaussian(
    tuple(N//2 for _ in range(d)),
    std
)

# initializing search hamiltonian
hamiltonian:Operator = tensor(Operator(N),pow=d)

for pos in hamiltonian.positions():
    for neighbour in hamiltonian.neighbours(*pos):
        hamiltonian[pos,neighbour] = 1
        

for neighbour in hamiltonian.neighbours(*target):
    hamiltonian[target,neighbour] = gammas[d]
    hamiltonian[neighbour,target] = gammas[d]
    
    
# definig the walk as a quantum walk + choosing the right timescale
dt = hamiltonian.timescale()*0.01*1j
walk = CTW(particle,hamiltonian,dt)

# solving and running
walk.solve(1000)
walk.run(
    keep_scale=False,
    speed_up=3
)