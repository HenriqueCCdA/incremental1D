from Plot import plot_res, plot_f, plot_E
from Modulo_elasticidade_var import modE
from Newton_raphson import incremental, nao_incremental, incremental_cons

def loads(q: float, n:int = 100):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 00/00/0000
    ---------------------------------------------------------------------------
    loads : caregamento
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    q  - carga maxima
    df - incremento de carga
    n  - numero de cargas
    ---------------------------------------------------------------------------
    Saida:
    ---------------------------------------------------------------------------
    retorna os passo de carga
    ---------------------------------------------------------------------------
    OBS:
    **************************************************************************
    '''
    f = []
    for i in range(0, n + 1):
        if i <= 10:
            f.append(q/10.0*i)
        elif 10 < i <= 20:
            f.append(q)
        elif 20 < i <= 29:
            f.append(-(q/10.0)*(i-20) + q)
        else:
            f.append(0.0)
    return f

def linear(AL: float, F: float, u:float, n: int):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 00/00/0000
    ---------------------------------------------------------------------------
    linear : solucao linear
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    AL     - A/L, onde A é a area e L é o comprimento
    F      - vetor F no passo de tempo n+1
    u      - deslocamento(atual)
    n      - passo de tempo
    ---------------------------------------------------------------------------
    Saida:
    ---------------------------------------------------------------------------
    u - deslocamentos atualizados
    E - modulo de elasticidade atualizado
    ---------------------------------------------------------------------------
    OBS:
    **************************************************************************
    '''
    E = modE(n , u)
    K = E*AL
    u = F/K

    # ... checa o equilibrio
    E = modE(n,u)*AL
    K = E*AL
    R = F - K*u
    print(f"Linear : passo de carga {n}: |F- Ku| = {abs(R):.6e}")

    return u, E


def main():

    # area transversal
    A: float = 1.0
    # comprimento
    L: float = 1.0
    AL = A/L

    # forcao aplicada ao longo do tempo
    step_loads = loads(-1.1, n = 40)

    # --- Nao incremental

    # deslocamento inicial
    u = 0.0

    res_n_inc = [u]
    E_n_inc = [modE(0, 0)]
    for n in range(1, len(step_loads)):
        F = step_loads[n]
        # newton-rapshon
        u, E = nao_incremental(AL, F, u, n)
        # armazenando o modulo de elasticidade atualizado
        E_n_inc.append(E)
        # armazenando a resposta
        res_n_inc.append(u)
    # ------------------------------------------------------------------------

    # --- Incremental

    # deslocamento inicial
    u = 0.0

    res_inc = [u]
    E_inc = [modE(0, 0)]
    for n in range(1, len(step_loads)):
        F = step_loads[n]
        # newton-rapshon
        u, E = incremental(AL, F, u, n)
        # armazenando o modulo de elasticidade atualizado
        E_inc.append(E)
        # armazenando a resposta
        res_inc.append(u)
    # ------------------------------------------------------------------------

    # --- Linear

    # deslocamento inicial
    u = 0.0

    res_linear = [u]
    E_linear = [modE(0, 0)]
    for n in range(1, len(step_loads)):
        F = step_loads[n]
        # atualizando o modulo de elasticidade
        # linear
        u, E = linear(AL, F, u, n)
        # armazenando o modulo de elasticidade atualizado
        E_linear.append(E)
        # armazenando a resposta
        res_linear.append(u)
    # ------------------------------------------------------------------------

    # --- incrementa equacao constitutiva

    # deslocamento inicial
    u = 0.0

    res_inc_cons = [u]
    E_inc_cons = [modE(0, 0)]
    for n in range(1, len(step_loads)):
        # df
        if n == 0:
            dF = step_loads[n]
        else:
            dF = step_loads[n] - step_loads[n-1]
        # modulo de elasticidade atuaizada em funcao do tempo
        E = modE(n, u)
        # linear
        u, E = incremental_cons(AL, dF, u, n)
        # armazenando o modulo de elasticidade atualizado
        E_inc_cons.append(E)
        # armazenando a resposta
        res_inc_cons.append(u)
    # ------------------------------------------------------------------------

    plot_res(res_inc, res_n_inc, res_linear, res_inc_cons, L)
    plot_f(step_loads)
    plot_E((E_inc,E_n_inc, E_inc_cons, E_linear))

if __name__ == "__main__":
    main()