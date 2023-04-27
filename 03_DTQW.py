from QW.qw import *

# creating a quantum walk in 1 dimension over a grid of size N
N=101

# initializing the quantum particle on position N//2 with spin +
particle = tensor(
    QParticle(N,N//2),
    QParticle(2,0) #|+>
)
print(repr(particle))

# coin operator
C = Operator(2)
C.array = np.array([
    [1,1],
    [1,-1]
])

I = Operator(N).eye()
C = tensor(I,C)

# translation operator
T = tensor(Operator(N),Operator(2))
for x in range(1,N-1):
    left_minus = tensor(QParticle(N,x-1),QParticle(2,1))    # |x-1,->
    right_plus = tensor(QParticle(N,x+1),QParticle(2,0))    # |x+1,+>
    minus = tensor(QParticle(N,x),QParticle(2,1))           # |x,->
    plus = tensor(QParticle(N,x),QParticle(2,0))            # |x,+>
    
    T+=QParticle.dag_product(left_minus,minus)                  # |x-1,-><x,-|
    T+=QParticle.dag_product(right_plus,plus)                   # |x+1,+><x,+|

U = T*C

# performing walk
walk = DTW(particle,U)
walk.solve(N//2)

# removing spin information from particle
particle:Particle = walk.particle.project_on_space(remove_index=1)

# plotting results
X = range(0,N,2)
Y = particle.distribution()[::2]
plt.plot(X,Y)
plt.show()