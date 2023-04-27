import numpy as np
from .tensor import Tensor


class Operator(Tensor):
    
    def zeros(self)->np.ndarray:
        return np.zeros((self.space_dim,self.space_dim),dtype=complex)
    
    def eye(self):
        """
        Returns:
            Operator: identity operator
        """
        out = Operator(self.dim)
        out.array = np.eye(self.space_dim,dtype=complex)
        return out
 
    def __getitem__(self,coords:tuple):
        
        if isinstance(coords[0],int) and isinstance(coords[1],int):
            coords = ((coords[0],),(coords[1],))
            
        assert isinstance(coords[0],tuple), "getitem and setitem methods require 2 tuples : opeator[pos2,pos1]"
        
        pos1 = coords[0]
        pos2 = coords[1]
        
        return self.array[self.to_index(*pos1),self.to_index(*pos2)]
    
    def __setitem__(self,coords:tuple,__v:complex)->None:
        
        if isinstance(coords[0],int) and isinstance(coords[1],int):
            coords = ((coords[0],),(coords[1],))
            
        assert isinstance(coords[0],tuple), "getitem and setitem methods require 2 tuples : opeator[pos1,pos2]"
        
        pos1 = coords[0]
        pos2 = coords[1]
        
        self.array[self.to_index(*pos1),self.to_index(*pos2)] = __v
        return
    
    
    def timescale(self)->float:
        """we don't know how fast the simulation is going to change.
        That's why we search or the highest transition rate in the hamiltonian
        and use it's inverse as a timescale.

        Returns:
            float: PI/2 * 1/max_rate
        """
        w_0 = self._find_max_ham_coeff()
        return 2*np.pi/w_0#T = 1/f = 2pi/w_0
    
    def _find_max_ham_coeff(self)->float:
        return np.max(np.abs(
            self.array - np.diag(np.diag(self.array))
        ))
    
    def laplacian(self):
        """
        Returns:
            Operator: Laplacian operator corresponding to the adjacency matrix
            (diagonal terms are adjusted so that the sum of each column / row is 0)
        """
        out = self.copy()
        out.array -= np.diag(np.diag(out.array)) #remove diagonal
        out.array -= np.diag(np.sum(out.array,axis=0)) #subtract row sums
        
        return out

    
    def herm_conj(self):
        """
        Returns:
            Operator: Hermitian conjugate of the operator
        """
        out = self.copy()
        out.array = np.transpose(np.conj(out.array))
        return out
