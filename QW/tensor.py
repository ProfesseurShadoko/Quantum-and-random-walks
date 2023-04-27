from .wrappers import abstract

import numpy as np
import itertools


class SpaceMissMatch(Exception):
    """Raised when trying to add / mulitply two objects that don't correspond to the same hilbrtian space"""
    
    def __init__(self, dims1, dims2,msg:str="") -> None:
        dim1 = "x".join((str(dim) for dim in dims1))
        dim2 = "x".join((str(dim) for dim in dims2))
        
        super().__init__(f"Dimensions do not match ({dims1} <> {dims2}). "+msg)
        
        
class Tensor:
    """A Tensor is a numpy array that keeps track of the ⊗ as a list of dimensions (.dim attribute)
    """
    
    def __init__(self,dim:int)->None:
        """Intitializes a particle, setting it's wavefunction to zeros

        Args:
            dim (int or list): dimension of the space (of list of the dimensions
            of the tensor space)
        """
        if type(dim)==list:
            self.dim = dim
        elif type(dim) == int:
            self.dim = [dim]
        else:
            self.dim = list(dim)
        
        self.array = self.zeros()
    
    @property
    def space_dim(self)->int:
        """Dimension of the Hilbertian Space that is a result of the tensor products"""
        out=1
        for dim in self.dim:
            out*=dim
        return out
    
    @abstract
    def zeros(self)->np.ndarray:
        """Depending on if the tensor is an operator or a particle, it should return a vector of zeros or a matrix of zeros.
        Should be reimplemented in children classes."""
        #return np.zeros((self.space_dim,self.space_dim),dtype=complex)
        #return np.zeros(self.space_dim,dtype=complex)
        pass
    
    def __str__(self):
        return str(self.to_array())
        
    def __repr__(self) -> str:
        typ = "wave_function" if self.is_particle() else "operator"
        return f"< {type(self).__name__} (tensor-{typ}) in space ({self.dim_to_str()}) >"
    
    def to_array(self)->np.ndarray:
        return self.array
    
    def dim_to_str(self)->str:
        """10x10x10 for example"""
        return  'x'.join((str(dim) for dim in self.dim))
        
    def __add__(self,__o):
        """returns addition as same type as the first operator"""
        if not isinstance(__o,Tensor):
            raise TypeError
        
        if self.dim != __o.dim:
            raise SpaceMissMatch(self.dim,__o.dim,"Unable to make addition.")
        
        if self.is_particle() != __o.is_particle():
            raise TypeError("Cannot add operator to wave_function")
        
        out = self.copy()
        out.array += __o.array
        
        return out
    
    def __mul__(self,__o):
        if isinstance(__o,Tensor):
            
            #print("Attempting multiplication between : ")
            #print(f"\t- {repr(self)} --> is_particle : {self.is_particle()}")
            #print(f"\t- {repr(__o)}  --> is_particle : {__o.is_particle()}")
        
            if self.dim != __o.dim:
                raise SpaceMissMatch(self.dim,__o.dim,"Unable to make multiplication.") #here, dim is the size of the hilbertian space, not the np.ndarray
            
            if self.is_particle() and __o.is_particle():
                raise TypeError("Cannot multiply two particles") #add projector ?
            
            if not self.is_particle() and not __o.is_particle():
                out = self.copy()
                out.array = self.array.dot(__o.array)
                return out
            
            if self.is_particle() and not __o.is_particle():
                out = self.copy()
                out.array = __o.array.dot(self.array) #WRONG ORDER HERE, IGNOR IT ?
                return out
            
            if not self.is_particle() and __o.is_particle():
                out = __o.copy()
                out.array = self.array.dot(__o.array)
                return out
            
        
        else: #float multiplication
            out = self.copy()
            out.array = self.array * __o
            return out
    
    def __rmul__(self,__o):#called by float multiplication
        """only used is case of multiplication of the tensor by a floating number. Simply calls __mul__ with the right order."""
        if not isinstance(__o,Tensor):
            return self.__mul__(__o)
        
    
    def __truediv__(self,__o):
        return self.__mul__(1/__o)
    
    def __sub__(self,__o):
        return self.__add__(-1*__o)
    
    def __neg__(self)->'Tensor':
        return self.__mul__(-1)
    
    @abstract
    def __getitem__(self,idx:tuple):
        pass
    
    @abstract
    def __setitem__(self,idx:tuple):
        pass
    
    def __tensor__(self,__o):
        """Allows to apply the ⊗ operation (through the tensor() function)"""
        if not isinstance(__o,Tensor):
            raise TypeError
        
        out = self.copy()
        out.dim = self.dim + __o.dim
        out.array = np.kron(self.array,__o.array)
        return out
    
    def to_index(self,*args)->int:
        """Translates the position (potential multi-dimensionnal tuple) in argument to and integer (the corresponding position on the array that the tensor represents)

        Returns:
            int: array[to_index(x,y,z)] = tensor[x,y,z]
        """
        assert len(args)==len(self.dim) and all([arg<dim for arg,dim in zip(args,self.dim)]),\
            f"Invalid arguments ({';'.join((str(arg) for arg in args))}) in tensor space {self.dim_to_str()}"
    
        out=0
        dim_fact=1
        
        for arg,dim in zip(reversed(args),reversed(self.dim)):
            out+=arg*dim_fact
            dim_fact*=dim
        return out
    
    def to_tuple(self,idx:int)->tuple:
        """Translates and index (the position on the array that the tensor represents) to the position (potential multi-dimensionnal tuple)

        Args:
            idx (int)
        
        Returns:
            tuple
        """
        out=[]
        dim_factor=self.space_dim
        for dim in self.dim:
            dim_factor//=dim
            out.append(
                idx//dim_factor
            )
            idx = idx%dim_factor
        return tuple(out)
    

        
    def positions(self)->itertools.product:
        """

        Returns:
            Iterable: allows to iterate over all the positions of the tensor space.
        
        For {1,2,3}⊗{a,b} => {
            (1,a),
            (1,b),
            (2,a),
            (2,b),
            (3,a),
            (3,b)
        }
        """
        ranges = [range(dim) for dim in self.dim]
        return itertools.product(*ranges)
    
    def neighbours(self,*coords:int)->list:
        """

        Returns:
            list: allows to iterate over all the neighbours of the positions in the tensor space.
        
        For SPACE = {1,2,3,4,5}⊗{1,2,3,4,5} and coords = 3,3 => [
            (3,4),
            (3,2),
            (2,3),
            (4,3)
        ]
        """
        assert len(coords)==len(self.dim) and all([arg<dim for arg,dim in zip(coords,self.dim)]),\
            f"Invalid arguments ({';'.join((str(arg) for arg in coords))}) in tensor space {self.dim_to_str()}"
        
        def is_ok(new_coords:tuple)->bool:
            """check if we do not go outside the lattice"""
            return all((x>=0 and x<dim for x,dim in zip(new_coords,self.dim)))
        
        out = []
        for i in range(len(self.dim)):
            point_up=[]
            point_down = []
            for k in range(len(self.dim)):
                up = 1 if i==k else 0
                point_up.append(coords[k]+up)
                point_down.append(coords[k]-up)
                
            if is_ok(point_up):
                out.append(tuple(point_up))
            if is_ok(point_down):
                out.append(tuple(point_down))
                  
        return out
    
    def is_particle(self)->bool:
        """when applying basic operations on tensors, it is important to keep track of the identity of the used objects
        We check if an object is a vector or an operator by checking the shape of the np array. 
        For a state, it is a 1d array. For an operator, it is a 2d matrix.
        
        Returns:
            bool: len(array.shape)==2
        """
        d = len(self.array.shape)
        return d==1
        

    def copy(self):
        """Performs a copy of the object (copying the list of dimensions and an array)

        Args:
            dtype (type, optional): type of the object in output. Defaults to None => output object has the same type as the copied object
        Return
        """
        out = type(self).__new__(type(self))
            
        out.dim = [dim for dim in self.dim]
        out.array = self.array.copy()
        return out
    
    @staticmethod
    def dag_product(p1,p2)->'Tensor':
        """

        Args:
            p1 (Tensor): target of the projection
            p2 (Tensor): source of the projection

        Returns:
            Tensor (operator): |p1><p2|
        """
        assert p1.is_particle() and p2.is_particle(), "Dag can only be computed between two particles"
        out = p1.copy()
        out.array = np.outer(p1.array,p2.array)
        return out
    
    
    
        

def tensor(*args:Tensor,pow:int=1) -> Tensor:
    """Operates tensor product ⊗ (the order is important !)

    Args:
        *args : Tensor
        power : int => puts the result of the tensor product of args to the power pow

    Returns:
        Tensor: A⊗B⊗C⊗...
    """
    out = args[0].copy()
    
    for i,arg in enumerate(args):
        if i==0:
            continue
        out = out.__tensor__(arg)
    
    if pow != 1:
        pow_list = [out.copy() for _ in range(pow)]
        return tensor(*pow_list)
    
    return out



