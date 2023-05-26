from .particle import Particle

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

class ParticlePlotter:
    
    use_cbar = False
    
    _cbar = None
    _warned = False
    
    @staticmethod
    def plot_1D(particle:Particle,ax:plt.Axes):
        assert len(particle.dim)==1, f"Cannot plot in 1D with dimensions too high ({particle.dim_to_str()})"
        dist = particle.distribution()
        
        ax.plot(
            [pos[0] for pos in particle.positions()],
            [dist[pos] for pos in particle.positions()],
        )

    @staticmethod
    def plot_2D_surface(particle:Particle,ax:plt.Axes):
        
        assert len(particle.dim)==2, f"Cannot plot in 2D with incompatible dimensions ({particle.dim_to_str()})"
        X = np.arange(particle.dim[0])
        Y = np.arange(particle.dim[1])
        X,Y = np.meshgrid(X,Y)
        Z=particle.distribution()
        
        ax.contourf(
            X,Y,Z,cmap="plasma"
        )
        
        ax.set_xlabel("x-position")
        ax.set_ylabel("y-position")
        plt.title("Particle distriution over time (u)")
        plt.colorbar()
       
    @staticmethod 
    def plot_2D(particle:Particle,ax:Axes3D):
        
        assert len(particle.dim)==2, f"Cannot plot in 2D with incompatible dimensions ({particle.dim_to_str()})"
        X = np.arange(particle.dim[0])
        Y = np.arange(particle.dim[1])
        X,Y = np.meshgrid(Y,X)
        Z=particle.distribution()
        
        ax.plot_surface(X,Y,Z,cmap="plasma")
    
    @staticmethod   
    def plot_3D(particle:Particle,ax:Axes3D):
        dist = particle.distribution()
        X = [pos[0] for pos in particle.positions()]
        Y = [pos[1] for pos in particle.positions()]
        Z = [pos[2] for pos in particle.positions()]
        VAL = dist
        
        m=np.max(VAL)
        RADIUS = 10*(0.1 + VAL / m * 3)
        
        p1=ax.scatter(X,Y,Z,c=VAL,cmap="plasma",s=RADIUS)
        
        if not ParticlePlotter._warned:
            ParticlePlotter._warned = True
            print("[log] : WARNING when plotting in 3D, you might get an error with the colorbar duplicating itself. This is why the colorbar is disabled. You can reable it with :\n\tParticlePlotter.use_cbar=True")
        
        if ParticlePlotter.use_cbar:
            try:
                ParticlePlotter._cbar.remove()
            except:
                pass
            ParticlePlotter._cbar = plt.colorbar(p1,label="probability")
        
    @staticmethod
    def plot(particle:Particle,ax:Axes3D):
        if len(particle.dim)==1:
            ParticlePlotter.plot_1D(particle,ax)
        if len(particle.dim)==2:
            ParticlePlotter.plot_2D(particle,ax)
        if len(particle.dim)==3:
            ParticlePlotter.plot_3D(particle,ax)
    

    
    
    
        