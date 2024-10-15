#!/usr/bin/env python
# coding: utf-8

from io import FileIO
from sympy import *
from IPython.display import display, Math, Markdown, Latex
init_printing(use_latex='mathjax')
# from sympy import diff as D

ReturnMatrices = true

a,b,c,d,e,f,g,h,x,y,x,z,u,v,w,k,t,r,l,ℓ = symbols('a b c d e f g h x y x z u v w k t r \\ell \\ell')
α,β,γ,δ,θ,φ,φ1,φ2 = symbols('alpha beta gamma delta theta varphi varphi_1 varphi_2')
θ1,θ2,θ3 = symbols('theta_1 theta_2 theta_3')
kα,kβ,kγ,kδ,kθ,kφ = symbols('k_{\\alpha} k_{\\beta} k_{\\gamma} k_{\\delta} k_{\\theta} k_{\\varphi}')
kθ1,kθ2,kθ3 = symbols('k_{\\theta_1} k_{\\theta_2} k_{\\theta_3}')
kx,ky,kz,kφ1,kφ2 = symbols('k_x k_y k_z k_{\\varphi_1} k_{\\varphi_2}')
Ck = {α:kα,β:kβ,γ:kγ,δ:kδ,θ:kθ,φ:kφ,x:kx,y:ky,z:kz,φ1:kφ1,φ2:kφ2,θ1:kθ1,θ2:kθ2,θ3:kθ3}
ℓα,ℓβ,ℓγ,ℓδ,ℓθ,ℓφ = symbols('\ell_{\\alpha} \ell_{\\beta} \ell_{\\gamma} \ell_{\\delta} \ell_{\\theta} \ell_{\\varphi}')
ℓx,ℓy,ℓz,ℓφ1,ℓφ2 = symbols('\ell_x \ell_y \ell_z \ell_{\\varphi_1} \ell_{\\varphi_2}')
ℓθ1,ℓθ2,ℓθ3 = symbols('\ell_{\\theta_1} \ell_{\\theta_2} \ell_{\\theta_3}')
Ls = {α:ℓα,β:ℓβ,γ:ℓγ,δ:ℓδ,θ:ℓθ,φ:ℓφ,x:ℓx,y:ℓy,z:ℓz,φ1:ℓφ1,φ2:ℓφ2,θ1:ℓθ1,θ2:ℓθ2,θ3:ℓθ3}

x1,x2,y1,y2,z1,z2 = symbols('x_1 x_2 y_1 y_2 z_1 z_2')
αt,βt,γt,δt,θt,φt = symbols('\\dot{\\alpha} \\dot{\\beta} \\dot{\\gamma} \\dot{\\delta} \\dot{\\theta} \\dot{\\varphi}')
θ1t,θ2t,θ3t = symbols('\\dot{\\theta_1} \\dot{\\theta_2} \\dot{\\theta_3}')
x1t,x2t,y1t,y2t,z1t,z2t = symbols('\\dot{x_1} \\dot{x_2} \\dot{y_1} \\dot{y_2} \\dot{z_1} \\dot{z_2}')
xt,yt,zt = symbols('\\dot{x} \\dot{y} \\dot{z}')
Ve = {α:αt,β:βt,γ:γt,δ:δt,θ:θt,θ1:θ1t,θ2:θ2t,θ3:θ3t,φ:φt,x:xt,y:yt,z:zt,x1:x1t,x2:x2t,y1:y1t,y2:y2t,z1:z1t,z2:z2t}

def subscrito2(primeiro: None, segundo: None):
    return symbols( latex(k)+ '_{' + latex(primeiro) + '.' + latex(segundo) +'}')

def MD_Matriz(st, mt, row, col):
    V = st + "\\begin{bmatrix}"
    for i in range(row):
        for j in range(col):
            V += latex(mt[i,j])
            if j != col-1:
                V += ' & '
            else:
                V += '\\\ '
    V += "\end{bmatrix}"
    return V

def PrintMD_MultDOF(Cg,f,r,F,J,K,L):
    display( Markdown("### Matrizes F e J") )
    FJ = MD_Matriz("F = ", F, r, f) + MD_Matriz("\qquad J = ", J, r, r)
    display(Markdown('$ %s  $' %FJ))
    
    printaK( J, F, K, f, Cg )

    # Impressão dos Coeficientes de Velocidade
    display( Markdown("### Coeficientes de Velocidade") )
    MK = [latex(i) for i in K]
    KS = []
    for i in range(r):
        for j in range(f):
            KS.append( latex(Ck[ Cg[j] ])  + latex( Cg[f+i] ) )
    KS = [i.replace('}',' ')+'}' for i in KS]
    KS = [i.replace('}',' ')+'}' for i in KS]
    LN = ""
    for i in range(len(KS)):
        LN += KS[i]+" = "+MK[i]+" \qquad "
    display(Markdown('$ %s $' %LN))

    try:
        printaL( J, F, K, f, Cg )
    except:
        display(Markdown('### 🔴🦇 O cálculo das matrizes de resolução deram errados, continuando o código... 🦇🔴'))
        # display(Markdown('Tente dexar a matriz no mesmo formato dessa:'))
        # display(Markdown('$$\\begin{bmatrix} 🦇 & 🦇 & 0 & 0 \\\\ 🦇 & 🦇 & 0 & 0 \\\\ 🦇 & 0 & 🦇 & 🦇 \\\\ 🦇 & 0 & 🦇 & 🦇 \\end{bmatrix}$$'))
        # display(Markdown('Ou dessa:'))
        # display(Markdown('$$\\begin{bmatrix} 🦇 & 🦇 & 0 & 0 \\\\ 🦇 & 🦇 & 0 & 0 \\\\ 0 & 🦇 & 🦇 & 🦇 \\\\ 0 & 🦇 & 🦇 & 🦇 \\end{bmatrix}$$'))

    # Impressão dos Coeficientes da Aceleração
    display( Markdown("### Coeficientes da Aceleração") )
    ML = [latex(i) for i in L]
    MS = []
    for i in range(r):
        for j in range(f):
            MS.append( latex(Ls[ Cg[j] ])  + latex( Cg[f+i] ) )
    MS = [i.replace('}',' ')+'}' for i in MS]
    LN = ""
    for i in range(len(KS)):
        LN += MS[i]+" = "+ML[i]+" \qquad "
    display(Markdown('$ %s $' %LN))
    return

# Fora zé Vampiro
def printaK( J, F, K, f, Cg ):
    if len(J) == 4:
        inversa = latex(1/simplify(det(J))) + latex(simplify((J**-1)*det(J)))
        invesa_negativa = latex(-1/simplify(det(J))) + latex(simplify((J**-1)*det(J)))
        display( Markdown('### Inversa de J') )
        display(Math(inversa))

    if (f==1):
        if len(J) == 4:
            display( Markdown('### Matrizes de resolução velocidades (Não somos robô  🦇) 😈') )
            display(Math(latex(Matrix([[Ck[Cg[f + 0]]], [Ck[Cg[ f + 1]]]])) + '=' + (invesa_negativa) + latex(F)))
        elif len(J) == 16:
            J_provisorio_1 = J[:2,:2]
            J_provisorio_2 = J[2:4,2:4]

            J_intermediario = latex(Ck[Cg[f]]) + latex(J[2:4, 0]) + '+' + latex(Ck[Cg[f + 1]]) + latex(J[2:4, 1])
            
            inversa_1 = latex(Matrix([[Ck[Cg[f + 0]]], [Ck[Cg[ f + 1]]]])) + '=' + latex(-1/simplify(det(J_provisorio_1))) + latex(simplify((J_provisorio_1**-1)*det(J_provisorio_1))) + latex(F[0:2, 0])
            inversa_2 = latex(Matrix([[Ck[Cg[f + 2]]], [Ck[Cg[ f + 3]]]])) + '=' + latex(-1/simplify(det(J_provisorio_2))) + latex(simplify((J_provisorio_2**-1)*det(J_provisorio_2))) + '(' + latex(F[2:4, 0]) + '+' + J_intermediario + ')'

            display( Markdown('### Matrizes de resolução velocidades (Não somos robô  🦇) 😈') )
            display(Math(inversa_1))
            display(Math(inversa_2))

    elif f == 2:

        if len(J) == 4:
            display( Markdown('### Matrizes de resolução velocidades (Não somos robô  🦇) 😈') )
            
            subscrito_1 = latex(Matrix([ [ (subscrito2(Cg[0], Cg[f])) ], [ (subscrito2(Cg[0], Cg[f + 1])) ] ]))
            display(Math( subscrito_1 + '=' + (invesa_negativa)+latex(F[0:,0])))
            
            subscrito_2 = latex(Matrix([ [ subscrito2(Cg[1], Cg[f]) ], [ subscrito2(Cg[1], Cg[f + 1]) ] ]))
            display(Math( subscrito_2 + '=' + (invesa_negativa)+latex(F[0:,1])))

        elif len(J) == 16:
            J_provisorio_1 = J[:2,:2]
            J_provisorio_2 = J[2:4,2:4]
            
            J_intermediario1 = latex(subscrito2(Cg[0], Cg[f + 0])) + latex(J[2:4, 0]) + '+' + latex(subscrito2(Cg[0], Cg[f + 1])) + latex(J[2:4, 1])
            J_intermediario2 = latex(subscrito2(Cg[1], Cg[f + 0])) + latex(J[2:4, 0]) + '+' + latex(subscrito2(Cg[1], Cg[f + 1]))  + latex(J[2:4, 1])
            
            subscrito_1 = latex(Matrix([[subscrito2(Cg[0], Cg[f])], [subscrito2(Cg[0], Cg[f + 1])]]))
            subscrito_2 = latex(Matrix([[subscrito2(Cg[1], Cg[f])], [subscrito2(Cg[1], Cg[f + 1])]]))
            inversa_1_1 = subscrito_1 + '=' + latex(-1/simplify(det(J_provisorio_1))) + latex(simplify((J_provisorio_1**-1)*det(J_provisorio_1))) + latex(F[0:2, 0])
            inversa_1_2 = subscrito_2 + '=' + latex(-1/simplify(det(J_provisorio_1))) + latex(simplify((J_provisorio_1**-1)*det(J_provisorio_1))) + latex(F[0:2, 1])

            subscrito_1 = latex(Matrix([[subscrito2(Cg[0], Cg[f + 2])], [subscrito2(Cg[0], Cg[f + 3])]]))
            subscrito_2 = latex(Matrix([[subscrito2(Cg[1], Cg[f + 2])], [subscrito2(Cg[1], Cg[f + 3])]]))
            inversa_2_1 = subscrito_1 + '=' + latex(-1/simplify(det(J_provisorio_2))) + latex(simplify((J_provisorio_2**-1)*det(J_provisorio_2))) + '(' + latex(F[2:4, 0]) + '+' + J_intermediario1 + ')'
            inversa_2_2 = subscrito_2 + '=' + latex(-1/simplify(det(J_provisorio_2))) + latex(simplify((J_provisorio_2**-1)*det(J_provisorio_2))) + '(' + latex(F[2:4, 1]) + '+' + J_intermediario2 + ')'

            display( Markdown('### Matrizes de resolução velocidades (Não somos robô  🦇) 😈') )
            display(Math(inversa_1_1))
            display(Math(inversa_1_2))
            display(Math(inversa_2_1))
            display(Math(inversa_2_2))

    elif f == 3:
        if len(J) == 4:
            display( Markdown('### Matrizes de resolução velocidades (Não somos robô  🦇) 😈') )
            
            subscrito_1 = latex(Matrix([ [ (subscrito2(Cg[0], Cg[f])) ], [ (subscrito2(Cg[0], Cg[f + 1])) ] ]))
            display(Math( subscrito_1 + '=' + (invesa_negativa)+latex(F[0:,0])))
            
            subscrito_2 = latex(Matrix([ [ subscrito2(Cg[1], Cg[f]) ], [ subscrito2(Cg[1], Cg[f + 1]) ] ]))
            display(Math( subscrito_2 + '=' + (invesa_negativa)+latex(F[0:,1])))

            subscrito_3 = latex(Matrix([ [ subscrito2(Cg[2], Cg[f]) ], [ subscrito2(Cg[2], Cg[f + 1]) ] ]))
            display(Math( subscrito_3 + '=' + (invesa_negativa)+latex(F[0:,2])))

def printaL( J, F, K, f, Cg ):
    def Derivative2( cord1, cord2, cord3):
        return Derivative( symbols(latex(k) + '_{' + latex(cord1) + '.' + latex(cord2) + '}'), cord3 )

    matriz_1 = [ [ [ Derivative2(cord_principal, cord_secund, cord_principal2) for cord_principal2 in Cg[:f]] for cord_secund in Cg[f:] ] for cord_principal in Cg[:f] ]
    matriz_2 = [ [ [ Derivative2(cord_principal, cord_secund, cord_secund2) for cord_secund2 in Cg[f:] ] for cord_secund in Cg[f:] ] for cord_principal in Cg[:f] ]

    display( Markdown('### Matrizes de resolução acelerações (Não somos robô  🦇) 😈') )

    for coord_principal in Cg[:f]:
        
        matriz_principal = Matrix( Cg[:f] )
        matriz_secundaria = Matrix( Cg[f:] )
        eq = latex( Matrix([ [symbols( latex(l)+'_{' + latex(coord_principal) + '.' + latex(coord_secund) + '}' )]  for coord_secund in Cg[f:] ] ) )
        lMatrix = Math( eq + '=' + latex( 1/Ve[coord_principal] ) + '(' + latex(Matrix(matriz_1[ Cg.index(coord_principal) ])) + latex( Matrix([[Ve[item]] for item in Cg[:f]]) ) + '+' + latex( Matrix(matriz_2[ Cg.index(coord_principal) ]) ) + latex( Matrix([[Ve[item]] for item in Cg[f:]]) )  + ')' ) 
        
        display( Markdown('##### Matriz L') )
        display( lMatrix )

        printMatrix = []

        for coordenada in Cg:
            indice_matriz_principal = Cg.index(coord_principal)

            for coord_secund in Cg[f:]:
                indice_matriz_secund = Cg.index(coord_secund ) - f

                ka = K[indice_matriz_secund:indice_matriz_secund+1 , indice_matriz_principal:indice_matriz_principal+1]
                ka= ka[0]

                derivada = latex( diff( ka, coordenada ) )
                
                eq1 = latex(Derivative2( coord_principal, coord_secund, coordenada ))

                # display( Math(eq1 + '=' + derivada ) )

                printMatrix.append(eq1 + '=' + derivada +" \qquad ")

        display( Markdown('##### Derivadas') )
        display( Markdown( '$ %s $' %' '.join(printMatrix)) )
        display()

        

def PrintMD_OneDOF(Cg,r,F,J,K,L,f):
    display( Markdown("### Matrizes F e J") )

    MF = [latex(i) for i in F]
    VF = "F = \\begin{Bmatrix}"

    for i in range(len(MF)):
        VF += MF[i]
        if i != len(MF)-1:
            VF += '\\\ '
    VF += "\end{Bmatrix}"

    LN = VF + MD_Matriz("\qquad J = ", J, r, r)
    display(Markdown('$ %s $' %LN))

    try:
        printaK( J, F, K, f, Cg )
    except:
        display(Markdown('### 🔴🦇 O cálculo das matrizes de resolução deram errados, continuando o código... 🦇🔴'))
        # display(Markdown('Tente dexar a matriz no mesmo formato dessa:'))
        # display(Markdown('$$\\begin{bmatrix} 🦇 & 🦇 & 0 & 0 \\\\ 🦇 & 🦇 & 0 & 0 \\\\ 🦇 & 0 & 🦇 & 🦇 \\\\ 🦇 & 0 & 🦇 & 🦇 \\end{bmatrix}$$'))
        # display(Markdown('Ou dessa:'))
        # display(Markdown('$$\\begin{bmatrix} 🦇 & 🦇 & 0 & 0 \\\\ 🦇 & 🦇 & 0 & 0 \\\\ 0 & 🦇 & 🦇 & 🦇 \\\\ 0 & 🦇 & 🦇 & 🦇 \\end{bmatrix}$$'))


    # Impressão dos Coeficientes de Velocidade
    display( Markdown("### Coeficientes de Velocidade") )
    MK = [latex(i) for i in K]
    KS = [latex(Ck[i]) for i in Cg[1:]]
    LN = ""
    for i in range(len(KS)):
        LN += KS[i]+" = "+MK[i]+" \qquad "
    display(Markdown('$ %s $' %LN))

    # Impressão dos Coeficientes da Aceleração
    display( Markdown("### Coeficientes da Aceleração") )
    Ms = Matrix([ Ls[i] for i in Cg[1:] ])
    ML = [latex(i) for i in L]
    MS = [latex(i) for i in Ms]
    LN = ""
    for i in range(len(KS)):
        LN += MS[i]+" = "+ML[i]+" \qquad "
    display(Markdown('$ %s $' %LN))
    return

def All_Matrix(Cg,Eq,f):
    Eq = Matrix(Eq)
    J = Eq.jacobian(Cg[f:])     # Jacobiano do sistema
    F = Eq.jacobian([Cg[:f]])   # Obtenção da matriz F
    K = simplify(-(J**-1)*F)    # Obtenção da matriz K
    if f == 1:
        Ks = Matrix([ Ck[i] for i in Cg[1:] ])
        L = simplify( K.jacobian([Cg[0]]) + K.jacobian([Cg[f:]])*Ks )
    else:
        L = Matrix([])
        Vv = [ Ve[i] for i in Cg ]
        for i in range(f):          # Obtenção da matriz L
            Lprov = K.col(i).jacobian(Cg[:f])*((1/Vv[i])*Matrix(Vv[:f])) + K.col(i).jacobian(Cg[f:])*((1/Vv[i])*Matrix(Vv[f:]))
            L = L.row_join(Lprov)
    return F, K, J, L

def MecSolve(Cg,Eq,flag=0):     # Cg,Vg Coords e Veloc Generalizadas
    f = len(Cg)-len(Eq)         # Graus de liberdade
    r = len(Eq)                 # Qtd. de Eq. de restrição

    F,K,J,L = All_Matrix(Cg,Eq,f)    

    if flag == true:
        return F, J, K, L

    if f == 1:
        PrintMD_OneDOF(Cg,r,F,J,K,L,f)
    else:
        PrintMD_MultDOF(Cg,f,r,F,J,K,L)