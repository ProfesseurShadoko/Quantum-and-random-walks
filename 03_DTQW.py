from QW.qw import *
# see PDF\Rapport.pdf in order to know what we are doing with the spin here


################################
# 1 DIMENSIONAL GRID OF SIZE N #
################################
N=101


######################################################
# INITIALIZATION OF THE PARTICLE ON POSITION (N/2,+) #
######################################################
particle = tensor(
    QParticle(N,N//2), #QParticle just changes the way the probability is computes (abs(_)² for quantum particle, just the coefficient of the array for classical particle)
    QParticle(2,0) #|+>
)
print(repr(particle))

#######################################
# INITIALIZATION OF THE COIN OPERATOR #
#######################################
C = Operator(2)
C.array = np.array([
    [1,1],
    [1,-1]
])

I = Operator(N).eye()
C = tensor(I,C)

##############################################
# INITIALIZATION OF THE TRANSLATION OPERATOR #
##############################################
T = tensor(Operator(N),Operator(2))
for x in range(1,N-1):
    left_minus = tensor(QParticle(N,x-1),QParticle(2,1))    # |x-1,->
    right_plus = tensor(QParticle(N,x+1),QParticle(2,0))    # |x+1,+>
    minus = tensor(QParticle(N,x),QParticle(2,1))           # |x,->
    plus = tensor(QParticle(N,x),QParticle(2,0))            # |x,+>
    
    T+=QParticle.dag_product(left_minus,minus)                  # |x-1,-><x,-|
    T+=QParticle.dag_product(right_plus,plus)                   # |x+1,+><x,+|


###
# #
###
U = T*C

#######################
# PERFORMING THE WALK #
#######################
walk = DTW(particle,U)
walk.solve(N//2)

# removing spin information from particle
particle:Particle = walk.particle.project_on_space(remove_index=1)

X = range(0,N,2)
Y = particle.distribution()[::2]
plt.plot(X,Y)
plt.show()