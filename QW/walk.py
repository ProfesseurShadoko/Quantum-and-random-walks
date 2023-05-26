

from .tensor import Tensor
from .particle import Particle
from .operators import Operator
from .particle_plotter import ParticlePlotter

from scipy import linalg as alg
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

class Walk:
    
    log:bool = True
    
    def __init__(self,particle:Particle,operator:Operator):
        self.particle_over_time=[particle]
        self.operator = operator
    
    @property
    def particle(self)->Particle:
        return self.particle_over_time[-1]
    
    @particle.setter
    def particle(self,part:Tensor):
        self.particle = [part]
    
    def __len__(self)->int:
        return len(self.particle_over_time)
    
    def add_step(self)->None:
        self.particle_over_time.append(self.operator * self.particle)
    
    def get_particle(self,time:int)->Particle:
        return self.particle_over_time[time]
    
    def solve(self,length:int=100):
        for i in range(length):
            if self.log:
                print(f"\r[log] : solving... {i/length:.2%}",end="",flush=True)
            self.add_step()
        if self.log:
            print("\r[log] : solving... 100%     ")
    
    def plot(self,i:int=-1,ax:plt.Axes=None):
        particle = self.get_particle(i)
        if ax==None:
            if len(particle.dim)==1:
                _ax = plt.subplot(1,1,1)
            else:
                _ax = plt.subplot(1,1,1,projection="3d")
        else:
            _ax = ax
        ParticlePlotter.plot(particle,_ax)
        
        if ax==None:
            plt.show()
        
    
    def run(self,keep_scale:bool=False,speed_up:int=1)->None:
        fig = plt.figure()
        
        if len(self.particle.dim)==1:
            ax = fig.add_subplot(1,1,1)
        else:
            ax = fig.add_subplot(1,1,1,projection="3d")
        
        def animate(t:int):
            t = int(t*speed_up)
            if t >= len(self):
                return
            fig.suptitle(f"Time : {t} (u)")
            ax.clear()
            
            d = len(self.particle.dim)
            
            if keep_scale:
                if d==1:
                    ax.set_ylim(0,1)
                if d==2:
                    ax.set_zlim(0,1)
            
            if d==1:
                ax.set_xlabel("x")
                ax.set_ylabel("p")
            if d==2:
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                ax.set_zlabel("p")
            if d==3:
                ax.set_xlabel("x")
                ax.set_ylabel("y")
                ax.set_zlabel("z")
                
            
            self.plot(t,ax)
        
        ani = anim.FuncAnimation(fig,animate,interval=10)
        plt.show()
        
    # analysis tools for search walk
    
    def follow_target(self,target:tuple,ax:plt.Axes=None,plot_ghost:bool=True,plot_topt:bool=True)->None:
        """plots the evolution of the amplitude of the targeted site

        Args:
            target (tuple): coordinates of the site that should be looked at
            ax (plt.Axes, optional): ax on which to plot. Defaults to None (figure is created  automatically and plt.show() is called).
            plot_ghost (bool, optional): ghost is the maximal amplitude over the lattice except for the target. Defaults to True.
            plot_topt (bool, optional): finds t_opt and plots a vertical red line. Defaults to True.
        """
        
        create_ax = (ax==None)
        
        if create_ax:
            fig = plt.figure()
            ax = plt.subplot(111)
            
        if type(target)==int:
            target=(target,)
        
        # probability of target over time
        X = [t for t in range(len(self.particle_over_time))]
        P = [self.get_particle(t).amplitude(*target) for t in X]
        T = self._find_first_peak(P)
        
        # lets find the ghost = maximal amplitude over the lattice with target removed
        P2=[]
        for t in X:
            abs_vals = np.abs(self.get_particle(t).array)
            abs_vals[self.particle.to_index(*target)] = 0
            P2.append(np.max(abs_vals)**2)
        
        ax.plot(X,P,label=f"P(x={target})")
        if plot_ghost:
            ax.plot(X,P2,color="grey",label=f"Fantôme")
        
        if T!=0 and plot_topt:
            ax.axvline(x=T,color="r",linestyle=":",label=f"Premier pic : {T} (10e-2 u)")
        
        if create_ax:
            plt.title("Probabilité du site cible au cours du temps")
            plt.xlabel(f"Temps ({10e-2:.1e}*u)")
            plt.ylabel("Probabilité")
            plt.legend()
            plt.show()
        
    def get_topt(self,target:tuple):
        
        if type(target)==int:
            target=(target,)
    
        P = [self.get_particle(t).amplitude(*target) for t in range(len(self.particle_over_time))]
        return self._find_first_peak(P)
    
    @staticmethod
    def _find_first_peak(signal:np.ndarray,treshold:float=0.6)->int:
        
        # signal above the treshold
        high_signal = signal > np.min(signal) + treshold * (np.max(signal) - np.min(signal))
        
        # compute connex components of high signal
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
                    
        # get the average position of each region
        avg_positions = [np.mean(region) for region in connexe_components]
       
        if len(avg_positions)==1:
            return np.argmax(signal)
        if len(avg_positions)==0:
            print("weird, no position...")
            return 0
        
        # distance between pics
        distances = np.diff(avg_positions)
        T = np.max(distances)
        
        # let's make sure that the first maximum falls on a period, because it is the first maximum that is most important
        first_max = np.argmax(signal[:int(T)])
        return first_max
        
    

class CountinousTimeWalk(Walk):
    
    def __init__(self, particle: Particle, H: Operator,dt:complex):
        """
        Args:
            particle (Particle): start particle
            operator (Operator): hamiltonian (evolution operator = exp(time_step * H))
            dt (complex): solving equation dΨ/dt = H*Ψ. You can find dt by using H.timescale() and dividing it by the discretization
            
        Perform the walk using evolutonary operator exp(dt*H)
        """
        self.H = H
        operator = H.copy()
        
        if self.log:
            print("[log] : computing evolutionnary operator... ",end="",flush=True)
        
        operator.array = alg.expm(dt * H.array)
        
        if self.log:
            print("OK")
            
        super().__init__(particle, operator)




        