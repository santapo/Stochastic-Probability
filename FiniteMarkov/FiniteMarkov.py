#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 03:39:34 2020

@author: santapo
"""
import numpy as np
import math

class FiniteMarkov:
    
    def __init__(self,n,P):
        self.n = n
        self.P = P
        
    def transitiveClosure(self,n,P):
        
        reach = np.copy(P) 
        for k in range(n): 
            for i in range(n): 
                for j in range(n): 
                    reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
        return reach
    
    def period_i(self,i,P):
        d = []
        for k in range(2,100):
            P_tmp = np.linalg.matrix_power(P,k)
            if (P_tmp[i,i]>0):
                d.append(k)
                
        tmp = d[0]  
        for j in range(1,len(d)):
            tmp2=math.gcd(tmp,d[j])
            period=tmp2
        return period
    
    def isIrreducible(self,n,P):
        flag=1
        reach = transitiveClosure(n+1,P)
        for i in range(n):
            for j in range(n):
                if (reach[i][j] == 0):
                    flag = 0
                    break
        return flag
    
    def isRegular(self,n,P):
        flag = 0
        if not(isIrreducible(n,P)):
            return flag
        if (period_i(0,P)==1):
            flag = 1
            return flag
        return flag
    
    def computePi(self,P):
        if P.shape == (1,1):
            return 1
        
        size = P.shape[0]
        dP = P.T - np.identity(size)
        A = np.append(dP,np.ones(size)).reshape(size+1,size)
        b = np.append(np.zeros(size),1)
    
        return np.linalg.solve(A.T.dot(A),A.T.dot(b))

    def Ergodic(self,n,P):
        if not isRegular(n,P):
            print("not Ergodic Markov chain")
            return
        return computePi(P)