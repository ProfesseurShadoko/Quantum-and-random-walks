from QW.qw import *

################################
# 1 DIMENSIONAL GRID OF SIZE N #
################################
N=101


##################################################
# INITIALIZATION OF THE PARTICLE ON POSITION N/2 #
##################################################
particle = Particle(N,N//2)
print(repr(particle))


#######################################
# INITIALIZATION OF THE WALK OPERATOR #
#######################################
operator = Operator(N)
for i in range(N-1):
    operator[i,i+1] = 0.5
    operator[i+1,i] = 0.5


#######################
# PERFORMING THE WALK #
#######################
walk = DTW(particle,operator)
walk.solve(N//2)

last_particle = walk.particle
print(f"Mean position : {last_particle.mean()}")
print(f"Variance : {last_particle.variance()}")
print(f"Std :  {last_particle.std()}")

X = range(0,N,2) # plotting the result of the walk (only even positions)
Y = walk.particle.distribution()[::2]
plt.plot(X,Y)
plt.show()

