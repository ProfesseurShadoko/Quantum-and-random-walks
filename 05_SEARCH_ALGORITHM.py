from QW.qw import *

###############
# WALK PARAMS #
###############
d=2
N=20
std=2
target=(13,13) # this has to be a tuple (even in 1D), for a walk in 1D, target = (13,)
use_uniform_particle = True


################################
# OPTIMIZED GAMMA FOR 1D,2D,3D #
################################
gammas=[
    np.nan,
    1.15,
    1.4,
    2.0
]


##################################
# INITIALIZATION OF THE PARTICLE #
##################################
particle:QParticle = tensor(QParticle(N),pow=d)
particle = particle.gaussian(
    tuple(N//2 for _ in range(d)),
    std
)
if use_uniform_particle:
    particle = particle.uniform()


#####################################
# INITIALIZATION OF THE HAMILTONIAN #
#####################################
hamiltonian:Operator = tensor(Operator(N),pow=d)

for pos in hamiltonian.positions():
    for neighbour in hamiltonian.neighbours(*pos):
        hamiltonian[pos,neighbour] = 1
        

for neighbour in hamiltonian.neighbours(*target):
    hamiltonian[target,neighbour] = gammas[d]
    hamiltonian[neighbour,target] = gammas[d]
    
# definig the walk as a quantum walk (1j factor) + choosing the right timescale
dt = hamiltonian.timescale()*0.01*1j
walk = CTW(particle,hamiltonian,dt)


######################
# COMPUTING THE WALK #
######################
print(f"Continous-Grover :",f"- {N=}",f"- target=({','.join(map(str,target))})",f"- gamma={gammas[d]}",f"- std of initial particle : {std}",sep="\n\t")
walk.solve(3000)
walk.run(
    keep_scale=False,
    speed_up=0.5
)

# plotting evolution of target site
walk.follow_target(target)
