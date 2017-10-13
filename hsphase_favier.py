''' Phasing Half-Sib families using the WCSP approach of Favier(2011)

Test data corresponds to phase information:

T1: 1 * * 0 1 * 1
T2: * * * 0 1 * 0
T3: 0 * * * 0 * 0

Informative pairs of markers:

T1: (0,3),(3,4),(4,6)
T2: (3,4),(4,6)
T3: (0,4),(4,6)

Genetic distance is assumed to be 0.1 cM between each marker pair
Can be instanciated by a @rawphase parameter

'''
from __future__ import print_function
import numpy as np
import Numberjack as nj

test_phase={'T1':[1  ,'*','*',0  ,1,'*',1],
            'T2':['*','*','*',0  ,1,'*',0],
            'T3':[0  ,'*','*','*',0,'*',0],
            'T4':['*','*','*','*','*','*','*'] }

class PhaseData(object):
    ''' 
Data to test phase inference.
'''
    def __init__(self,rawphase):
        self.rawphase=rawphase
        ## info_pairs stores pairs of inform. markers as keys and [N+,N-] as values.
        self.info_pairs={}
        ## For each offspring
        for k,v in self.rawphase.items():
            ## get informative markers ...
            infomk=[i for i,x in enumerate(v) if x!='*']
            if len(infomk)<2:
                continue
            ## and corresponding successive pairs
            mypairs=zip(infomk[:-1],infomk[1:])
            ## fir each pair
            for p in mypairs:
                ## get phase information
                npair=int(v[p[0]]==v[p[1]]) ## + or - phase
                ## adds it to the pool
                try:
                    Nvec=self.info_pairs[p]
                except KeyError:
                    self.info_pairs[p]=[0,0]
                    Nvec=self.info_pairs[p]
                Nvec[1-npair]+=1
        ## informative markers are all mks seen in pairs
        info_mk=[]
        for k in self.info_pairs.keys():
            info_mk+=list(k)
        self.info_mk=np.unique(info_mk)
    @staticmethod
    def recrate(i,j,d=0.1):
        '''
        returns rec. rate between marker indices i and j
        assuming each adj. markers are separeted by d (default 0.1) cM
        '''
        dtot=abs(i-j)*d
        return 0.5*(1-exp(-d/50))

if __name__=='__main__':
    T=PhaseData(test_phase)
    print(*T.info_pairs.items())
