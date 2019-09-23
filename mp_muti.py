from multiprocessing import Process
from multiprocessing import cpu_count
import numpy as np
class multical:
    def multi(A,B,Ncpu=cpu_count()):
        shapeA=A.shape[0]
        Part=shapeA//Ncpu
        lPart=shapeA%Ncpu
        for i in range(Ncpu-1):
            p=Process(target=A[i*Part:i*Part+Part].dot(B))
            p_list.append(p)
            p.start()
        [ap.join() for ap in p_list]

