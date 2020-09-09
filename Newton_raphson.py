def nao_incremental(F: float, u: float, n: int, K
                   , tol: float = 1.e-11,max_it: int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 05/09/2020
    ---------------------------------------------------------------------------
    nao_incremental : Alg Newton-Raphson normal
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    F      - vetor F no passo de tempo n+1
    u      - deslocamento(atual)
    n      - passo de tempo
    K      - matriz (Objeto)
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
        # matriz de rigidez atualizado (K(n+1,i))
        Ki = K.update(n, u)
        # residuo ( R(n+1,i) = F(n+1,i) - K(n+1,i) * u(n+1,i) )
        R = F - Ki*u
        # ---- convergencia
        if abs(R) < tol:
            break

        # du(n+1,i+1) = K^-1(n+1,i) * R
        du = R/Ki
        # u(n+1,i+1) = u(n+1,i) + du(n+1,i+1)
        u += du

    print(f"NR : passo de carga {n}: |F- Ku| = {abs(R):.6e}")
    return u, K.E

def incremental(F: float, u: float, n: int, K
                , tol: float = 1.e-11,max_it: int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 05/09/2020
    ---------------------------------------------------------------------------
    incremental : Alg Newton-Raphson incremental
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    F      - vetor F no passo de tempo n+1
    u      - deslocamento(atual)
    n      - passo de tempo
    K      - matriz (Objeto)
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
        # matriz de rigidez atualizado (K(n+1,i))
        Ki = K.update(n, u0 + du)
        # residuo (dR(n+1,i) = F(n+1,i) - K(n+1,i)*u(n) - K(n+1,i)*du(n+1,i))
        dR = F - Ki*u0 - Ki*du
        # convergencia
        if abs(dR) < tol:
            break
        # ddu(n+1,i+1) = K^-1(n+1,i) * dR
        ddu = dR/Ki

        # du(n+1,i+1) = du(n+1,i) + ddu(n+1,i+1)
        du += ddu

    print(f"NRI : passo de carga {n}: |F- Ku| = {abs(dR):.6e}")
    return u0 + du, K.E

def incremental_cons(dF: float, u: float, n: int, K
                , tol: float = 1.e-11,max_it: int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 05/09/2020
    ---------------------------------------------------------------------------
    incremental_cons : Alg Newton-Raphson com equacao constitutiva incremental
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    dF     - vetor F(n+1) - F(n)
    u      - deslocamento(atual)
    n      - passo de tempo
    K      - matriz (Objeto)
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
        # matriz de rigidez atualizado (K(n+1,i))
        Ki = K.update(n, u0 + du)
        # residuo (dR(n+1,i) = F(n+1,i) - F(n) - K(n+1,i)*du(n+1,i))
        dR = dF - Ki*du
        # convergencia
        if abs(dR) < tol:
            break

        # ddu(n+1,i+1) = K^-1(n+1,i) * dR
        ddu = dR/Ki

        # du(n+1,i+1) = du(n+1,i) + ddu(n+1,i+1)
        du += ddu

    print(f"NRCI : passo de carga {n}: |F- Ku| = {abs(dR):.6e}")
    return u0 + du, K.E
