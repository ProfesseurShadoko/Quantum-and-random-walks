
from time import perf_counter,sleep
import atexit

class Timer:
    
    def __init__(self) -> None:
        self.start_time = None
        self.pause_time = None
    
    def __str__(self) -> str:        
        return f"Timer : {self.get()}s"
    
    def __repr__(self) -> str:
        status = "RUNNING" if self.pause_time is None else "PAUSED"
        return f"< Timer : {self.get()}s ({status}) >"
    
    def start(self)->None:
        if self.start_time is None:
            self.start_time = perf_counter()
    
    def pause(self):
        if self.pause_time is None:            
            self.pause_time = perf_counter()
    
    def stop(self):
        self.pause()
    
    def resume(self):
        if self.pause_time is not None:
            self.start_time += perf_counter() - self.pause_time
            self.pause_time = None
    
    def reset(self):
        self.__init__()
    
    def restart(self):
        self.reset()
        self.start()
    
    def get(self)->int:
        if self.start_time is None:
            return 0
        if self.pause_time is None:
            return int(perf_counter() - self.start_time)
        if self.pause_time is not None:
            return int(self.pause_time - self.start_time)
    
    def format(self)->str:
        time = self.get()
        
        h = str(time // 3600)
        m = str((time % 3600) // 60)
        s = str(time%60)
        
        if len(h) == 1:
            h = "0" + h
        if len(m) == 1:
            m = "0" + m
        if len(s) == 1:
            s = "0" + s
        
        return f"{h}:{m}:{s}"
    
    def display(self):
        print(f"\n\nExecution time --> {self.format()}\n")
        

    
            