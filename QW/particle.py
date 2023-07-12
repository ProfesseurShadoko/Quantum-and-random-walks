import numpy as np
from .tensor import tensor,Tensor


class Particle(Tensor):
    
    def __init__(self,dim:int,fock:int=None):
        super().__init__(dim)
        
        if fock != None:
            self[fock]=1
    
    def fock(self,*args:int)->np.ndarray:
        """Returns a particle array with a 1 at position corresponding to tensor coordinates"""
        assert len(args)==len(self.dim) and all([arg<dim for arg,dim in zip(args,self.dim)]),\
            f"Invalid arguments ({';'.join((str(arg) for arg in self.args))})\
            in tensor space {'x'.join((dim for dim in self.dim))}"
            
        out = self.zeros()
        n = self.to_index(*args)
        out[n] = 1
        return out
    
    def zeros(self)->np.ndarray:
        return np.zeros(self.space_dim,dtype=complex)
    
    def uniform(self):
        """
        Returns:
            Particle : particle filled with ones and then normalized
        
        """
        out = self.copy()
        out.array[None]=1
        out.normalize()
        return out
    
    def gaussian(self,mean:tuple,std:float):
        """
        Args:
            mean (float)
            std (float)
            
        Returns:
            Particle : particle filled with gaussian distribution depending on distance to mean, and then normalized
        """
        out = self.copy()
        
        for position in self.positions():
            out[position]=np.exp(-sum((x-y)**2 for x,y in zip(position,mean)) / (2*std**2))
        out.normalize()
        return out
    
    def __getitem__(self,idx:tuple)->complex:
        if isinstance(idx,int):
            idx = (idx,)
        return self.array[self.to_index(*idx)]
    
    def __setitem__(self,idx:tuple,__v:complex)->None:
        if type(idx)==int:
            idx = (idx,)
            
        self.array[self.to_index(*idx)] = __v
    
    def amplitude(self,*args)->float: #overwritten when using quantum objects : needs to be squared
        """probability amplitude of a tensor state represented by the coordinates in input.
        When looping over the positions for exemple, run p = tensor.amplitude(*position) for each position.

        Returns:
            float: usually the absolute value of the array
        """
        return abs(self[args])
    
    def distribution(self)->np.ndarray:
        """Probability (amplitude) distribution. Computed using self.amplitude().

        Returns:
            np.ndarray: distribution[x,y,z] = P(x,y,z)
            
        Usable when plotting with meshgrid(X,Y) for dimensions > 2
        """
        
        if hasattr(self,"_dist"):
            return self._dist
            
        dist = np.zeros(self.dim)
        for state in self.positions():
            dist[state] = self.amplitude(*state)
        
        self._dist = dist
        self._flatten_dist = dist.flatten()
        return dist
    
    def mean(self)->np.ndarray:
        """Computes distribution again and deduces it's mean 

        Returns:
            np.ndarray: mean position of the particle
        """
        dist = self.distribution()
        return sum(dist[pos]*np.array(pos) for pos in self.positions()) / self.norm()
    
    def variance(self)->float:
        e = self.mean()
        dist = self.distribution()
        
        def distance_squared(pos1,pos2):
            return sum(
                (x1-x2)**2 for x1,x2 in zip(pos1,pos2)
            )
        
        
        return sum(
            distance_squared(pos,e)*dist[pos] for pos in self.positions()
        )
    
    def std(self)->float:
        return np.sqrt(self.variance())
    
    def norm(self)->float:
        """
        Returns:
            float: sum of the amplitudes of every position
        """
        return sum(self.amplitude(*pos) for pos in self.positions())
    
    def normalize(self)->None:
        """in place modification
        """
        self.array/=self.norm()
    
    def normalized(self):
        """returns new normalized object
        """
        return self / self.norm()
    
    def mesure(self)->tuple:
        """
        Returns:
            tuple: random position chosen according to the distribution probability
        """
        self.distribution()
        flat_index = np.random.choice(len(self._flatten_dist),p=self._flatten_dist)
        return self.to_tuple(flat_index)
    
    def project_on_space(self,remove_index:list)->'Particle':
        """Let's say we study the position of the particle, but we introduce an additionnal degree of freedom (like the spin).
        We want to access the amplitude of |n> but we have two amplitudes |n,+> & |n,->.
        This function allows the Particle to be projected to a smaller tensor space, by adding the amplitudes, and resetting
        all phase (because |n,+> does not interfer with |n,->)

        Args:
            index (int): number of tensor space to remove
            
        Returns:
            Particle: particle in projected space (space with tensor_space{index} removed) as Particle !!! (even if origin particle was QParticle)
            In other terms, returns the probability distribution represented as a Particle.
        """
        _n = self.dim[remove_index]
        new_dim = [d for d in self.dim]
        new_dim.pop(remove_index)
        
        projected = Particle(new_dim)
        
        for pos in self.positions():            
            new_pos = list(pos)
            new_pos.pop(remove_index)
            new_pos = tuple(new_pos)
            
            #carefull, there is a correct way to add depending on the metric
            """
            Simply adding doesn't work :
            Start :
            |n> :: a*|n,+> + b*|n,->
            P(|n>) = ||n>|² = a²+b²
            Projected : 
            |n> :: (a+b)|n>
            P(|n>)=(a+b)**2=a²+b²+2ab
            
            This is because of the Quantum metrics. Let's start over :
            |n>  :: a*|n,+> + b*|n,->
            -->
            |n> :: f(a,b)|n>
            Who is f ?
            P(|n>) = amplitude(f(a,b)) = amplitude(a)+amplitude(b)
            f(a,b) = reverse_amplitude(amplitude(a)+amplitude(b))
            
            """
            projected[new_pos]+=self.amplitude(*pos)
        
        #
        
        return projected
            
    
    
class QParticle(Particle):
    """inherits from Particle, is a tensor. Only difference is the amplitude,
    wich is the square of the norm (quantum physics)
    """
    
    def amplitude(self, *args)->float:
        return super().amplitude(*args)**2
    
    def normalize(self) -> None:
        self.array/=np.sqrt(self.norm())
    
    def normalized(self)->'QParticle':
        return self / np.sqrt(self.norm())
