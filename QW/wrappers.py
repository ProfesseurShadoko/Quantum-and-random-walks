

from functools import wraps

class AbstractFunction(Exception):
    
    def __init__(self, func) -> None:
        super().__init__(f"Cannot call abstract function {func.__name__}. Must be defined in inherited class !")
        
        
def abstract(func):
    
    @wraps(func)
    def wrapper(self,*args,**kwargs):
        raise AbstractFunction(func)
    
    return wrapper

