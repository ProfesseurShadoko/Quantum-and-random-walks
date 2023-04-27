

from .tensor import Tensor
from .particle import Particle
from .operators import Operator
from .particle_plotter import ParticlePlotter

from scipy import linalg as alg
import matplotlib.pyplot as plt
import matplotlib.animation as anim

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




        