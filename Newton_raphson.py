from Modulo_elasticidade_var import modE

def nao_incremental(AL: float, F: float, u: float, n: int
                   , tol: float = 1.e-11,max_it: int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 00/00/0000
    ---------------------------------------------------------------------------
    nao_incremental : Alg Newton-Raphson normal
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    AL     - A/L, onde A é a area e L é o comprimento
    F      - vetor F no passo de tempo n+1
    u      - deslocamento(atual)
    n      - passo de tempo
    tol    - tolerancia do metodo nao linear
    max_it - numero maximo de iteracoes
    ---------------------------------------------------------------------------
    Saida:
    ---------------------------------------------------------------------------
    u - deslocamentos atualizados
    E - modulo de elasticidade atualizado
    ---------------------------------------------------------------------------
    OBS:
    n - tempo
    i - iteracoes
    **************************************************************************
    '''
    for i in range(max_it):
        # modulo de elesticidade atuaizada em funcao de u
        E = modE(n, u)
        # matriz de rigidez atualizado (K(n+1,i))
        K = E*AL
        # residuo ( R(n+1,i) = F(n+1,i) - K(n+1,i) * u(n+1,i) )
        R = F - K*u
        # ---- convergencia
        if abs(R) < tol:
            break

        # du(n+1,i+1) = K^-1(n+1,i) * R
        du = R/K
        # u(n+1,i+1) = u(n+1,i) + du(n+1,i+1)
        u += du

    print(f"NR : passo de carga {n}: |F- Ku| = {abs(R):.6e}")
    return u, E

def incremental(AL: float, F: float, u: float, n: int
                , tol: float = 1.e-11,max_it: int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 00/00/0000
    ---------------------------------------------------------------------------
    incremental : Alg Newton-Raphson incremental
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    AL     - A/L, onde A é a area e L é o comprimento
    F      - vetor F no passo de tempo n+1
    u      - deslocamento(atual)
    n      - passo de tempo
    tol    - tolerancia do metodo nao linear
    max_it - numero maximo de iteracoes
    ---------------------------------------------------------------------------
    Saida:
    ---------------------------------------------------------------------------
    u - deslocamentos atualizados
    E - modulo de elasticidade atualizado
    ---------------------------------------------------------------------------
    OBS:
    n - tempo
    i - iteracoes
    **************************************************************************
    '''

    du: float = 0.0 # du(n+1, 0)
    u0: float = u   # u(n)
    for i in range(max_it):
        # modulo de elesticidade atuaizada em funcao de u
        E = modE(n, u0 + du)
        # matriz de rigidez atualizado (K(n+1,i))
        K = E*AL
        # residuo (dR(n+1,i) = F(n+1,i) - K(n+1,i)*u(n) - K(n+1,i)*du(n+1,i))
        dR = F - K*u0 - K*du
        # convergencia
        if abs(dR) < tol:
            break
        # ddu(n+1,i+1) = K^-1(n+1,i) * dR
        ddu = dR/K

        # du(n+1,i+1) = du(n+1,i) + ddu(n+1,i+1)
        du += ddu

    print(f"NRI : passo de carga {n}: |F- Ku| = {abs(dR):.6e}")
    return u0 + du, E

def incremental_cons(AL: float, dF: float, u: float, n: int
                , tol: float = 1.e-11,max_it: int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 00/00/0000
    ---------------------------------------------------------------------------
    incremental_cons : Alg Newton-Raphson com equacao constitutiva incremental
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    AL     - A/L, onde A é a area e L é o comprimento
    dF     - vetor F(n+1) - F(n)
    u      - deslocamento(atual)
    n      - passo de tempo
    tol    - tolerancia do metodo nao linear
    max_it - numero maximo de iteracoes
    ---------------------------------------------------------------------------
    Saida:
    ---------------------------------------------------------------------------
    u - deslocamentos atualizados
    E - modulo de elasticidade atualizado
    ---------------------------------------------------------------------------
    OBS:
    n - tempo
    i - iteracoes
    **************************************************************************
    '''

    du: float = 0.0 # du(n+1, 0)
    u0: float = u   # u(n)
    for i in range(max_it):
        # modulo de elesticidade atuaizada em funcao de u
        E = modE(n, u0 + du)
        # matriz de rigidez atualizado (K(n+1,i))
        K = E*AL
        # residuo (dR(n+1,i) = F(n+1,i) - F(n) - K(n+1,i)*du(n+1,i))
        dR = dF - K*du
        # convergencia
        if abs(dR) < tol:
            break

        # ddu(n+1,i+1) = K^-1(n+1,i) * dR
        ddu = dR/K

        # du(n+1,i+1) = du(n+1,i) + ddu(n+1,i+1)
        du += ddu

    print(f"NRCI : passo de carga {n}: |F- Ku| = {abs(dR):.6e}")
    return u0 + du, E
