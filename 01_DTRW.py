from QW.qw import *

# creating a random walk in 1 dimension over a grid of size N
N=101

# initializing the particle on position N//2
particle = Particle(N,N//2)
print(repr(particle))

# crating walk operator
operator = Operator(N)
for i in range(N-1):
    operator[i,i+1] = 0.5
    operator[i+1,i] = 0.5

# performing the walk (DTW stands for discrete time walk)
walk = DTW(particle,operator)
walk.solve(N//2)

# getting the result of the walk
last_particle = walk.particle

print(f"Mean position : {last_particle.mean()}")
print(f"Variance : {last_particle.variance()}")
print(f"Std :  {last_particle.std()}")

# plotting the result of the walk (only even positions)
X = range(0,N,2)
Y = walk.particle.distribution()[::2]
plt.plot(X,Y)
plt.show()

