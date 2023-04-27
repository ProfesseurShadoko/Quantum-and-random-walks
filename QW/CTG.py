from .walk import CountinousTimeWalk as CTW
from .tensor import tensor
from .operators import Operator
import numpy as np
import matplotlib.pyplot as plt
from .particle import QParticle

"""Les paramètres à étudier :
- N => la taille du réseau
- d => la dimension de la marche
- t => la position de la cible (proche du bord ou non ?)
- g => gamma
- tau => le temps idéal où faire la mesure pour un maximum de résultat
- pmax => la probabilité maximale obtenue (et la présence de fantômes)
"""



class CTG(CTW):
    """Cette classe correspond à la réalisation d'une CTQW ciblée."""
    speed = 0.01
    period_limit = 0.6
    
    
    def __init__(self,target:tuple,N:int=21,gamma:float=1.5,steps:int=1000) -> None:
        self.gamma = gamma
        self.d = len(target)
        self.N = N
        self.target = target
        self.dt = 1j*CTG.speed*2*np.pi
        
        super().__init__(self._get_gamma_particle(),self._get_gamma_hamiltonian(),self.dt)
        
        self.solve(steps)
    
    def _get_gamma_particle(self)->QParticle:
        particle = tensor(QParticle(self.N),pow=self.d)
        particle = particle.uniform()
        particle.normalize()
        return particle
    
    def _get_gamma_hamiltonian(self)->Operator:
        ham:Operator = tensor(Operator(self.N),pow=self.d)
        for pos in ham.positions():
            for neighbourg in ham.neighbours(*pos):
                ham[pos,neighbourg]=1
                
        for neighbour in ham.neighbours(*self.target):
            ham[self.target,neighbour]=self.gamma
            ham[neighbour,self.target]=self.gamma
        self.dt = 1j * ham.timescale() * CTG.speed
        return ham
    
    def strongest_site(self,time:int)->float:
        """It is hoped that te strongest site will be the target, and that the probiability will be high
        Args:
            particle (Particle)

        Returns:
            float:probability
        """
        return np.max(np.abs(self.get_particle(time).array))**2
    
    def second_strongest_site(self,time:int)->float:
        particle = self.get_particle(time)
        abs_vals = np.abs(particle.array)
        
        #remove the target site
        abs_vals[self.particle.to_index(*self.target)] = 0
        
        #find new maximum value
        return np.max(np.max(abs_vals))**2
        
    def target_prob(self)->np.ndarray:
        """returns the probability of the target over time
        Returns:
            np.ndarray
        """
        return [self.get_particle(t).amplitude(*self.target) for t in range(len(self.particle_over_time))]
    
    def first_peak(self)->int:
        """computes the time for the first peek of the walk

        Returns:
            float: T where P(w,T) is local mixima
        """
        P = self.target_prob()
        T = self.find_first_peak(P)
        return T
    
    def plot_target_evolution(self,ax:plt.Axes=None,show_second_strongest:bool=False,plot_topt_line:bool=False)->None:
        
        create_ax = (ax==None)
        
        if create_ax:
            fig = plt.figure()
            ax = plt.subplot(111)
            
        P = self.target_prob()
        P2 = [self.second_strongest_site(time) for time in range(len(P))]
        T = CTG.find_first_peak(P)
        X = [t for t in range(len(P))]
        
        ax.plot(X,P,label=f"P(x={self.target})")
        
        if show_second_strongest:
            ax.plot(X,P2,color="grey",label=f"Fantôme")
        
        
        if T!=0 and plot_topt_line:
            ax.axvline(x=T,color="r",linestyle=":",label=f"Premier pic : {T} (10e-2 u)")
        
        if create_ax:
            plt.title("Probabilité du site cible")
            plt.xlabel(f"Temps ({CTG.speed:.1e}*u)")
            plt.ylabel("Probabilité")
            plt.legend()
            plt.show()
    
    @staticmethod
    def find_first_peak(signal:np.ndarray)->int:
        """doesn't exactly finds the period of the signal, but the index
        of the first maximum, that corresponds approximatively to the 
        period of the signal

        Args:
            signal (np.ndarray)

        Returns:
            float: index of the first maximum
        """
        high_signal = signal > np.min(signal) + CTG.period_limit * (np.max(signal) - np.min(signal))
        
        #compute connex components of high signal
        connexe_components = []
        region = []
        for i, b in enumerate(high_signal):
            if b:
                region.append(i)
            else:
                if region==[]:
                    continue
                else:
                    connexe_components.append(region)
                    region=[]
                    
        #get the average position of each region
        avg_positions = [np.mean(region) for region in connexe_components]
       
        if len(avg_positions)==1:
            return np.argmax(signal)
        if len(avg_positions)==0:
            print("weird, no position...")
            return 0
        #distance between pics
        distances = np.diff(avg_positions)
        T = np.max(distances)
        
        #let's make sure that the first maximum falls on a period, because it is the first maximum that is most important
        first_max = np.argmax(signal[:int(T)])
        return first_max
        
    
    @staticmethod
    def nearest_square(root:int)->int:
        """
        find the nearest integer n so that n**2 is the closest to root**3
        """
        return int(round(root**(3/2)))
    
    @staticmethod
    def random_target(d:int,N:int)->tuple:
        """
        Args:
            d (int): 1D,2D,3D
            N (int): size of the grid

        Returns:
            tuple: random position, not too close to the border of the grid
        """
        out=[]
        for i in range(d):
            out.append(np.random.randint(N//3,2*N//3))
        return tuple(out)

class CTG_n(CTG):
    
    def _get_gamma_hamiltonian(self)->Operator:
        ham:Operator = tensor(Operator(self.N),pow=self.d)
        for pos in ham.positions():
            for neighbourg in ham.neighbours(*pos):
                ham[pos,neighbourg]=1
        
        self.dt = 1j * ham.timescale() * CTG.speed
             
        for neighbour in ham.neighbours(*self.target):
            ham[self.target,neighbour]=self.gamma
            ham[neighbour,self.target]=self.gamma
        
        return ham

class CTG_t(CTG):
    
    def _get_gamma_hamiltonian(self)->Operator:
        ham:Operator = tensor(Operator(self.N),pow=self.d)
        for pos in ham.positions():
            for neighbourg in ham.neighbours(*pos):
                ham[pos,neighbourg]=1
        
        self.dt = 1j * ham.timescale() * CTG.speed
             
        ham[self.target,self.target]=self.gamma
        
        return ham
    
