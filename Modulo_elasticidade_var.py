from math import log10

global tVar
global uVar

def modE(n: float, u: float):
    '''
    ***************************************************************************
    data criacao:     30/08/2020
    data modificacao: 00/00/0000
    ---------------------------------------------------------------------------
    modE_u : variacao do modulo de elestaticdade com o deslocamento
    ---------------------------------------------------------------------------
    Entrada:
    ---------------------------------------------------------------------------
    n     - passo de tempo
    u     - deslocamento
    ---------------------------------------------------------------------------
    Saida:
    ---------------------------------------------------------------------------
    E - modulo de elasticidade atualizado
    ---------------------------------------------------------------------------
    OBS:
    **************************************************************************
    '''

    uVar, tVar = False, True
    if tVar:
        gt = 4*log10(10*n+10)
    else:
        gt = 4

    if uVar:
        ui = abs(u)
        fx = (3.176*ui**4 - 8.3625*ui**3 + 7.9172*ui**2 - 3.2258*ui + 0.9994)
    else:
        fx = 1.0

    return gt*fx
