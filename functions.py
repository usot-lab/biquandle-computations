import numpy as np
import itertools
import dataset as data
from sympy import symbols, Matrix, zeros

def all_arrays(j, p):
    # Generates all ixj matrices over Z_p
    for flat_vals in itertools.product(range(p), repeat=j):
        yield np.array(flat_vals).reshape(j)
    return flat_vals

def inv(n,p):
    # Computes the inverse of an element n in Z_p
        for i in range(p):
            if (i*n)%p==1:
                return i

def alexander_bq(n, t, r):
    # Generates an Alexander biquandle over Z_n for the given invertible elements t and r.
    alex_over = np.zeros((n,n), dtype=('int'))
    alex_under = np.zeros((n,n), dtype=('int'))
    for i in range(n):
        for j in range(n):
            alex_under[i][j]=(t*i+(r-t)*j)%n
            alex_over[i][j]=(r*i)%n
    return alex_under,alex_over

def is_biquandle(under_mat, over_mat):
    # Checks if the given under and over operation matrices define a biquandle on Z_n
    over_mat = np.array(over_mat)
    under_mat = np.array(under_mat)
    n = over_mat.shape[0]
    for x in range(n):
        if over_mat[x, x] != under_mat[x, x]:
            return False
            
    seen1 = set()
    for a, b in itertools.product(range(n), repeat=2):
        s = (over_mat[b, a], under_mat[a, b])
        if s in seen1:
            return False
        seen1.add(s)
    for b in range(n):
        seen2 = set()
        for a in range(n):
            s = over_mat[a, b]
            if s in seen2:
                return False
            seen2.add(s)
    for b in range(n):
        seen3 = set()
        for a in range(n):
            s = under_mat[a, b]
            if s in seen3:
                return False
            seen3.add(s)

    for x, y, z in itertools.product(range(n), repeat=3):
        left1 = under_mat[under_mat[x, y], under_mat[z, y]]
        right1 = under_mat[under_mat[x, z], over_mat[y, z]]
        if left1 != right1:
            return False

        left2 = over_mat[under_mat[x, y], under_mat[z, y]]
        right2 = under_mat[over_mat[x, z], over_mat[y, z]]
        if left2 != right2:
            return False

        left3 = over_mat[over_mat[x, y], over_mat[z, y]]
        right3 = over_mat[over_mat[x, z], under_mat[y, z]]
        if left3 != right3:
            return False

    return True

def matrix_converter(matrix):
    # Converts a python matrix to a standard matrix.
    n = len(matrix)
    new_matrix = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==0 and j==0:
                new_matrix[n-1][n-1] = matrix[0][0]
            if i==0 and j!=0:
                new_matrix[n-1][j-1] = matrix[0][j]
            if i!=0 and j==0:
                new_matrix[i-1][n-1] = matrix[i][0]
            if i!=0 and j!=0:
                new_matrix[i-1][j-1] = matrix[i][j]
    return new_matrix

def is_biquandle_bracket(under_matrix, over_matrix, p, A_m, B_m, V_m, C_m, D_m, U_m):
    # For a given biquandle (with under and over matrices), and a given prime p, checks if the coefficient matrix [A_m B_m V_m C_m D_m U_m] defines a biquandle virtual bracket over Z_p.
    n = len(over_matrix)
    def inverse(n):
        return(hp.inv(n,p))
    def A(x,y):
        return A_m[x][y]
    def B(x,y):
        return B_m[x][y]
    def V(x,y):
        return V_m[x][y]
    def C(x,y):
        return C_m[x][y]
    def D(x,y):
        return D_m[x][y]
    def U(x,y):
        return U_m[x][y]
    def o(x,y):
        return over_matrix[x][y]
    def u(x,y):
        return under_matrix[x][y]
        
    delta = 0
    for i in range(n):
        t=0
        for j in range(n):
            if A(i,j)!=0 and C(i,j)!=0:
                delta = (-inv(A(i,j),p)*B(i,j)-inv(C(i,j),p)*D(i,j))%p
                t=1
                break
            elif B(i,j)!=0 and D(i,j)!=0:
                delta = (-A(i,j)*inv(B(i,j),p)-C(i,j)*inv(D(i,j),p))%p
                t=1
                break
        if t==1:
            break
    if delta==0:
        return (False, 0, 0)
    omega = (delta*A(0,0)+B(0,0)+V(0,0))%p
    if delta==0 or omega==0:
        return (False, 0, 0)

    # Checks if the given coefficient matrix satisfies the biquandle virtual axioms 1-23.
    all_good = True
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if inv(omega,p)!= (delta*C(z,z)+D(z,z)+U(z,z))%p: 
                    all_good = False
                    break
                if 1!= (A(x,y)*C(x,y)+V(x,y)*U(x,y))%p: #2.1
                    all_good = False
                    break
                if 1!= (B(x,y)*D(x,y)+V(x,y)*U(x,y))%p: #2.2
                    all_good = False
                    break
                if 0!= (A(x,y)*U(x,y)+V(x,y)*C(x,y))%p: #2.3
                    all_good = False
                    break
                if 0!= (B(x,y)*U(x,y)+V(x,y)*D(x,y))%p: #2.4
                    all_good = False
                    break
                if 0!= (delta*B(x,y)*D(x,y)+A(x,y)*D(x,y)+B(x,y)*C(x,y))%p: #2.5
                    all_good = False
                    break
                if 0!= (delta*A(x,y)*C(x,y)+A(x,y)*D(x,y)+B(x,y)*C(x,y))%p: #2.6
                    all_good = False
                    break
                if (A(x,y)*A(u(x,y),o(z,y))*A(y,z)+V(x,y)*A(u(x,y),o(z,y))*V(y,z))%p!=\
                (A(o(y,x),o(z,x))*A(x,z)*A(u(x,z),u(y,z)) + V(o(y,x),o(z,x))*A(x,z)*V(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*A(u(x,y),o(z,y))*B(y,z)+B(x,y)*A(u(x,y),o(z,y))*A(y,z)+\
                delta*B(x,y)*A(u(x,y),o(z,y))*B(y,z)++B(x,y)*A(u(x,y),o(z,y))*V(y,z)+\
                B(x,y)*B(u(x,y),o(z,y))*B(y,z)++B(x,y)*V(u(x,y),o(z,y))*B(y,z)+\
                +V(x,y)*A(u(x,y),o(z,y))*B(y,z))%p != (A(o(y,x),o(z,x))*B(x,z)*A(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*B(u(x,y),o(z,y))*A(y,z))%p !=\
                (A(o(y,x),o(z,x))*A(x,z)*B(u(x,z),u(y,z))+B(o(y,x),o(z,x))*A(x,z)*A(u(x,z),u(y,z))+\
                delta*B(o(y,x),o(z,x))*A(x,z)*B(u(x,z),u(y,z))+B(o(y,x),o(z,x))*A(x,z)*V(u(x,z),u(y,z))+\
                B(o(y,x),o(z,x))*B(x,z)*B(u(x,z),u(y,z))+B(o(y,x),o(z,x))*V(x,z)*B(u(x,z),u(y,z))+\
                V(o(y,x),o(z,x))*A(x,z)*B(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*V(u(x,y),o(z,y))*A(y,z))%p !=\
                (A(o(y,x),o(z,x))*A(x,z)*V(u(x,z),u(y,z))+V(o(y,x),o(z,x))*A(x,z)*A(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*A(u(x,y),o(z,y))*V(y,z)+V(x,y)*A(u(x,y),o(z,y))*A(y,z))%p!=\
                (A(o(y,x),o(z,x))*V(x,z)*A(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (B(x,y)*B(u(x,y),o(z,y))*A(y,z)+B(x,y)*V(u(x,y),o(z,y))*V(y,z))%p!=\
                (A(o(y,x),o(z,x))*B(x,z)*B(u(x,z),u(y,z))+V(o(y,x),o(z,x))*V(x,z)*B(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*B(u(x,y),o(z,y))*B(y,z)+V(x,y)*V(u(x,y),o(z,y))*B(y,z))%p!=\
                (B(o(y,x),o(z,x))*B(x,z)*A(u(x,z),u(y,z))+B(o(y,x),o(z,x))*V(x,z)*V(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (B(x,y)*B(u(x,y),o(z,y))*V(y,z)+B(x,y)*V(u(x,y),o(z,y))*A(y,z))%p!=\
                (A(o(y,x),o(z,x))*B(x,z)*V(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*B(u(x,y),o(z,y))*V(y,z))%p!=\
                (B(o(y,x),o(z,x))*B(x,z)*V(u(x,z),u(y,z))+B(o(y,x),o(z,x))*V(x,z)*A(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (V(x,y)*B(u(x,y),o(z,y))*A(y,z))%p!=\
                (A(o(y,x),o(z,x))*V(x,z)*B(u(x,z),u(y,z))+V(o(y,x),o(z,x))*B(x,z)*B(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*V(u(x,y),o(z,y))*B(y,z)+V(x,y)*B(u(x,y),o(z,y))*B(y,z))%p!=\
                (V(o(y,x),o(z,x))*B(x,z)*A(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (V(x,y)*V(u(x,y),o(z,y))*A(y,z))%p != (A(o(y,x),o(z,x))*V(x,z)*V(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (A(x,y)*V(u(x,y),o(z,y))*V(y,z))%p != (V(o(y,x),o(z,x))*V(x,z)*A(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (V(x,y)*B(u(x,y),o(z,y))*V(y,z))%p != (V(o(y,x),o(z,x))*B(x,z)*V(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
                if (V(x,y)*V(u(x,y),o(z,y))*V(y,z))%p != (V(o(y,x),o(z,x))*V(x,z)*V(u(x,z),u(y,z)))%p:
                    all_good = False
                    break
    if all_good:
        return(all_good, delta, omega)
    else:
        return(False, delta, omega)

def bvb_m(code_label, bvb):
    # For a given labeled peer code code_label of a virtual knotoid diagram, and a choosen biquandle virtual bracket contained in dataset, this function computes the biquandle virtual bracket matrix of code_label.
    u = symbols("u")
    for code in data.parsed_codesV:
        if code.label==code_label:
            writheV = writhe(code)
            if writheV<0:
                weight = (data.bvb[bvb].omega**(-writheV))%data.bvb[bvb].bvbp
            else:
                weight = (inv(data.bvb[bvb].omega,data.bvb[bvb].bvbp)**(writheV))%data.bvb[bvb].bvbp
            bvb_matrix = zeros(len(data.bvb[bvb].under))
            coloring_list, coloring_filt = color(code, data.bvb[bvb].bp ,data.bvb[bvb].under, data.bvb[bvb].over)
            free_ind = free_indices(code)
            for coloring in coloring_list:
                state_values = []
                for state_label in all_arrays(len(free_ind),3):
                    coef = state_coefficient(code, state_label, free_ind)
                    coef_val = coef_eval(coef,coloring,data.bvb[bvb].A,data.bvb[bvb].B,data.bvb[bvb].V,\
                                            data.bvb[bvb].C,data.bvb[bvb].D,data.bvb[bvb].U)
                    comp_no = state_comp_no(code,state_label,free_ind)
                    state_value = coef_val*(data.bvb[bvb].delta**comp_no)
                    state_values.append((weight*state_value)%data.bvb[bvb].bvbp)
                total_sum = 0
                for i in state_values:
                    total_sum=total_sum+i
                bvb_matrix[coloring[0],coloring[-1]] += u**(total_sum%data.bvb[bvb].bvbp)
    return(bvb_matrix)

def state_coefficient(peerCode, state_label, free_indices):
    # Computes the coefficient of a given state.
    pCode = peerCode.pCode
    signs = peerCode.signs
    index = peerCode.index
    c = len(pCode)
    n = 2*c
    coefficient_array = []
    for i in range(len(free_indices)):
        k = free_indices[i]
        y_k = abs(pCode[k])
        if pCode[k]>0 and signs[k] == "-":
            if state_label[i] == 0:
                coefficient_array.append(["D", 2*k+1, y_k])
            if state_label[i] == 1:
                coefficient_array.append(["C", 2*k+1, y_k])
            if state_label[i] == 2:
                coefficient_array.append(["U", 2*k+1, y_k])    
        if pCode[k]>0 and signs[k] == "+":
            if state_label[i] == 0:
                coefficient_array.append(["A", y_k, 2*k+1])
            if state_label[i] == 1:
                coefficient_array.append(["B", y_k, 2*k+1])
            if state_label[i] == 2:
                coefficient_array.append(["V", y_k, 2*k+1])
        if pCode[k]<0 and signs[k] == "+":
            if state_label[i] == 0:
                coefficient_array.append(["D", y_k+1, 2*k])
            if state_label[i] == 1:
                coefficient_array.append(["C", y_k+1, 2*k])
            if state_label[i] == 2:
                coefficient_array.append(["U", y_k+1, 2*k])
        if pCode[k]<0 and signs[k] == "-":
            if state_label[i] == 0:
                coefficient_array.append(["A", 2*k, y_k+1])
            if state_label[i] == 1:
                coefficient_array.append(["B", 2*k, y_k+1])
            if state_label[i] == 2:
                coefficient_array.append(["V", 2*k, y_k+1])
    return coefficient_array

def label_switcher(old,new,labels):
    # When a smoothing is utilized at a crossing, relabels the joining semi-arcs so that they have the same label.
    old_value = labels[old]
    new_value = labels[new]
    for i in range(len(labels)):
        if labels[i] == old_value:
            labels[i] = new_value
    return labels
    
def state_comp_no(peerCode, state_label, free_indices):
    # Computes the number of components in a state of a virtual knotoid diagram.
    pCode = peerCode.pCode
    signs = peerCode.signs
    index = peerCode.index
    c = len(pCode)
    n = 2*c
    labels = np.array(range(n), dtype = "int")
    if index >= 0:
        if abs(pCode[index])+1 == 2*len(pCode):
            labels = label_switcher((2*index)%n,(2*index+1)%n,labels)
        if abs(pCode[index])+1 != 2*len(pCode):
            labels = label_switcher(abs(pCode[index]),(abs(pCode[index])+1)%n,labels)
            for j in range(len(pCode)):
                if abs(pCode[j])+1 == 2*len(pCode):
                    labels = label_switcher((2*j)%n,(2*j+1)%n,labels)
    for k in range(c):
        y_k = abs(pCode[k])
        if signs[k] == "*":
            labels = label_switcher(y_k,(y_k+1)%n,labels)
            labels = label_switcher((2*k)%n,(2*k+1)%n,labels)
    for i in range(len(free_indices)):
        k = free_indices[i]
        #print("Şu anda ",k, "indisi için işlem yapıyorum.")
        y_k = abs(pCode[k])
        if pCode[k]>0 and signs[k] == "-":
            if state_label[i] == 0:
                labels = label_switcher((2*k)%n,y_k,labels)
                labels = label_switcher((y_k+1)%n,(2*k+1)%n,labels)
            if state_label[i] == 1:
                labels = label_switcher((2*k)%n,(y_k+1)%n,labels)
                labels = label_switcher(y_k,(2*k+1)%n,labels)
            if state_label[i] == 2:
                labels = label_switcher((2*k)%n,(2*k+1)%n,labels)
                labels = label_switcher(y_k,(y_k+1)%n,labels)   
        if pCode[k]>0 and signs[k] == "+":
            if state_label[i] == 0:
                labels = label_switcher((2*k)%n,(y_k+1)%n,labels)
                labels = label_switcher(y_k,(2*k+1)%n,labels)
            if state_label[i] == 1:
                labels = label_switcher((2*k)%n,y_k,labels)
                labels = label_switcher((y_k+1)%n,(2*k+1)%n,labels)
            if state_label[i] == 2:
                labels = label_switcher((2*k)%n,(2*k+1)%n,labels)
                labels = label_switcher(y_k,(y_k+1)%n,labels) 
        if pCode[k]<0 and signs[k] == "+":
            if state_label[i] == 0:
                labels = label_switcher((2*k)%n,y_k,labels)
                labels = label_switcher((y_k+1)%n,(2*k+1)%n,labels)
            if state_label[i] == 1:
                labels = label_switcher((2*k)%n,(y_k+1)%n,labels)
                labels = label_switcher(y_k,(2*k+1)%n,labels)
            if state_label[i] == 2:
                labels = label_switcher((2*k)%n,(2*k+1)%n,labels)
                labels = label_switcher(y_k,(y_k+1)%n,labels)
        if pCode[k]<0 and signs[k] == "-":
            if state_label[i] == 0:
                labels = label_switcher((2*k)%n,(y_k+1)%n,labels)
                labels = label_switcher(y_k,(2*k+1)%n,labels)
            if state_label[i] == 1:
                labels = label_switcher((2*k)%n,y_k,labels)
                labels = label_switcher((y_k+1)%n,(2*k+1)%n,labels)
            if state_label[i] == 2:
                labels = label_switcher((2*k)%n,(2*k+1)%n,labels)
                labels = label_switcher(y_k,(y_k+1)%n,labels)
    cn = len(np.unique(labels))
    return cn

def free_indices(peerCode):
    free_indices = [] 
    pCode = peerCode.pCode
    signs = peerCode.signs
    index = peerCode.index
    redundant_indices=[index]
    if abs(pCode[index])+1 != 2*len(pCode):
            for j in range(len(pCode)):
                if abs(pCode[j])+1 == 2*len(pCode):
                    redundant_indices.append(j)
    for i in range(len(pCode)):
        if signs[i] == "*":
            redundant_indices.append(i)
    for i in range(len(pCode)):
        if i not in redundant_indices:
            free_indices.append(i)
    return free_indices        
    
def equationDeriver (pCode,signs,index):
    # Derives the equations of labels obtained from each crossing in a given virtual knotoid diagram.
    c=len(pCode)
    y=list(map(abs,pCode))
    n=2*c
    eqs=[]
    redundantSemiarcs=[]
    if index < 0:
        for k in range(c):
            if pCode[k]>0 and signs[k]=="+":
                eqs.append(("over",2*k+1,y[k],2*k))
                eqs.append(("under",y[k],2*k+1,(y[k]+1)%n))
            elif pCode[k]>0 and signs[k]=="-":
                eqs.append(("over",y[k],2*k+1,(y[k]+1)%n))
                eqs.append(("under",2*k+1,y[k],2*k))
            elif pCode[k]<0 and signs[k]=="+":
                eqs.append(("over",2*k,y[k]+1,(2*k+1)%n))
                eqs.append(("under",y[k]+1,2*k,y[k]))
            elif pCode[k]<0 and signs[k]=="-":
                eqs.append(("over",y[k]+1,2*k,y[k]))
                eqs.append(("under",2*k,y[k]+1,(2*k+1)%n))
            elif signs[k]=="*":
                eqs.append(("virtual",2*k,(2*k+1)%(2*c),1))
                eqs.append(("virtual",y[k],(y[k]+1)%(2*c),1))
                redundantSemiarcs.append(y[k]+1)
                redundantSemiarcs.append(2*k+1)
    if index>=0 and y[index]==2*c-1:
        for k in range(c):
            if pCode[k]>0 and signs[k]=="+" and index!=k:
                eqs.append(("over",2*k+1,y[k],2*k))
                eqs.append(("under",y[k],2*k+1,y[k]+1))
            elif pCode[k]>0 and signs[k]=="-" and index!=k:
                eqs.append(("over",y[k],2*k+1,y[k]+1))
                eqs.append(("under",2*k+1,y[k],2*k))
            elif pCode[k]<0 and signs[k]=="+" and index!=k:
                eqs.append(("over",2*k,y[k]+1,2*k+1))
                eqs.append(("under",y[k]+1,2*k,y[k]))
            elif pCode[k]<0 and signs[k]=="-" and index!=k:
                eqs.append(("over",y[k]+1,2*k,y[k]))
                eqs.append(("under",2*k,y[k]+1,2*k+1))
            elif signs[k]=="*":
                eqs.append(("virtual",2*k,2*k+1,1))
                eqs.append(("virtual",y[k],y[k]+1,1))
                redundantSemiarcs.append(y[k]+1)
                redundantSemiarcs.append(2*k+1)
            eqs.append(("equal",2*index,2*index+1,1))
            redundantSemiarcs.append(2*index+1)
    elif index>=0 and y[index]!=2*c-1:
        for i in range(c):
            if y[i]==2*c-1:
                index2=i
        for k in range(c):
            if pCode[k]>0 and signs[k]=="+" and index!=k and index2!=k:
                eqs.append(("over",2*k+1,y[k],2*k))
                eqs.append(("under",y[k],2*k+1,y[k]+1))
            elif pCode[k]>0 and signs[k]=="-" and index!=k and index2!=k:
                eqs.append(("over",y[k],2*k+1,y[k]+1))
                eqs.append(("under",2*k+1,y[k],2*k))
            elif pCode[k]<0 and signs[k]=="+" and index!=k and index2!=k:
                eqs.append(("over",2*k,y[k]+1,2*k+1))
                eqs.append(("under",y[k]+1,2*k,y[k]))
            elif pCode[k]<0 and signs[k]=="-" and index!=k and index2!=k:
                eqs.append(("over",y[k]+1,2*k,y[k]))
                eqs.append(("under",2*k,y[k]+1,2*k+1))
            elif signs[k]=="*":
                eqs.append(("virtual",2*k,2*k+1,1))
                eqs.append(("virtual",y[k],y[k]+1,1))
                redundantSemiarcs.append(y[k]+1)
                redundantSemiarcs.append(2*k+1)
        eqs.append(("equal",y[index],y[index]+1,1))
        redundantSemiarcs.append(y[index]+1)
        eqs.append(("equal",2*index2,2*index2+1,1))
        redundantSemiarcs.append(2*index2+1)
        redundantSemiarcs.append(2*c-1)
        eqs.append(("equal",2*c-1,2*c-2,1))
    return eqs, redundantSemiarcs

def color (code, p, u_matrix, o_matrix):
    # Computes the counting number of a virtual knotoid diagram "code" with the given biquandle operation matrices over Z_p. 
    pCode = code.pCode
    signs = code.signs
    index = code.index

    def over (x,y):
        return o_matrix[x][y]
    def under (x,y):
        return u_matrix[x][y]
    
    equations, redundantSemiarcs = equationDeriver(pCode,signs,index)
    equations = list(dict.fromkeys(equations))

    coloring_number = 0
    coloring_list = []
    coloring_list_filt = []
    for assignment in itertools.product(range(p), repeat=2*len(pCode)):
        x = assignment
        all_good = True
        for op, i, j, k in equations:
            if op == "over":
                if over(x[i], x[j]) != x[k]:
                    all_good = False
                    break
            elif op == "under":
                if under(x[i], x[j]) != x[k]:
                    all_good = False
                    break
            elif op== "virtual":
                if x[i]!=x[j]:
                    all_good = False
                    break
            elif op== "equal":
                if x[i]!=x[j]:
                    all_good = False
                    break
        if all_good:
            coloring_list.append(x)
            redundantSemiarcs.sort(reverse = True)
            x_filtered = tuple(val for i, val in enumerate(x) if i not in redundantSemiarcs)
            coloring_list_filt.append(x_filtered)
            coloring_number=coloring_number+1
    return coloring_list, coloring_list_filt

def coef_eval(coefficient_array, coloring, A, B, V, C, D, U):
    # Evaluates the value of the coefficient of a state for the given coefficient matrix of a biquandle virtual bracket.
    coefficient = 1
    for label in coefficient_array:
        if label[0] == "A":
            coefficient = coefficient* A[coloring[label[1]]][coloring[label[2]]]
        if label[0] == "B":
            coefficient = coefficient* B[coloring[label[1]]][coloring[label[2]]]
        if label[0] == "V":
            coefficient = coefficient* V[coloring[label[1]]][coloring[label[2]]]
        if label[0] == "C":
            coefficient = coefficient* C[coloring[label[1]]][coloring[label[2]]]
        if label[0] == "D":
            coefficient = coefficient* D[coloring[label[1]]][coloring[label[2]]]
        if label[0] == "U":
            coefficient = coefficient* U[coloring[label[1]]][coloring[label[2]]]
    return coefficient

def writhe(code):
    # Computes the writhe of a given virtual knotoid diagram "code".
    pCode = code.pCode
    signs = code.signs
    index = code.index
    n = 0
    p = 0
    for i in range(len(pCode)):
        if i!=index:
            if (pCode[i]<0 and signs[i]=="-") or (pCode[i]>0 and signs[i]=="+"):
                p = p+1
            if (pCode[i]<0 and signs[i]=="+") or (pCode[i]>0 and signs[i]=="-"):
                n = n+1
    return(p-n)


