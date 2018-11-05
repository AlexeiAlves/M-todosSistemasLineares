
# coding: utf-8

# In[2]:


import math
import numpy as np
def substituicoes_retroativas(A, b):
    n = len(A)
    x = n * [0] 
    for i in range(0,n):
        for j in range(0,i+1):
            x[n-1-i]=x[n-1-i] + A[n-1-i][n-1-i+j]*x[n-1-i+j]
            
        x[n-1-i]=(b[n-1-i] - x[n-1-i])/(A[n-1-i][n-1-i])
    return x


# In[3]:


def gauss(A,b):
    n = len(A)
    for i in range(0,n-1):
        for j in range(i+1,n):
            b[j] += b[i] * -A[j][i]/A[i][i]
            A[j][:] = [a + b for a,b in zip (A[j][:],list(map(lambda x : x * (-A[j][i])/(A[i][i]),A[i][:])))]
            
    x = substituicoes_retroativas(A, b)
    #print(A)
    return x
    


# In[4]:


def gauss_det(A,b):
    n = len(A)
    det = A[0][0]
    for i in range(0,n-1):
        for j in range(i+1,n):
            b[j] += b[i] * -A[j][i]/A[i][i]
            A[j][:] = [a + b for a,b in zip (A[j][:],list(map(lambda x : x * (-A[j][i])/(A[i][i]),A[i][:])))]
        det = det*A[i+1][i+1]
    x = substituicoes_retroativas(A, b)
    return(x , det)


# In[5]:


def gauss_jordan(A,b):
    n=len(A)
    i=np.identity(n)
    A=np.concatenate((A,i),axis=1)
    for i in range(0,n):           
        b[i] = b[i]/A[i][i]
        A[i][:]=list(map(lambda x:x/A[i][i],A[i][:]))
        if(i==n-1):
            for k in range(0,i):
                    b[k]+=b[i] * -A[k][i]/A[i][i]
                    A[k][:]=[a+b for a,b in zip (A[k][:],list(map(lambda x:x*-A[k][i],A[i][:])))]
                    
        for j in range(i,n-1):
            if(i != n-1):
                b[j+1] += b[i] * -A[j+1][i]
                A[j+1][:] = [a + b for a,b in zip (A[j+1][:],list(map(lambda x : x * -A[j+1][i],A[i][:])))]
                    
            if(j==i and i>0) :
                for h in range(0,i):
                    
                    b[h]+=b[i] * -A[h][i]
                    A[h][:]=[a+b for a,b in zip (A[h][:],list(map(lambda x:x*-A[h][i],A[i][:])))]
                               
    return print(A[0:n,n:2*n],'\n\n', b)
        
    
            
        
                
            
        

