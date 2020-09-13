from math import log10


class ModElasticidade():

    def __init__(self, u_var: bool, t_var: bool, modE: float = 1.0):
        '''
        ***********************************************************************
        data criacao:     04/08/2020
        data modificacao: 00/00/0000
        -----------------------------------------------------------------------
        modE_u : variacao do modulo de elestaticdade com o deslocamento
        -----------------------------------------------------------------------
        Entrada:
        -----------------------------------------------------------------------
        u_var - variacao da modulo de elesticidade em funcao de x
        t_var - variacao da modulo de elesticidade em funcao de t(n)
        modE  - modulo de elasticidade inicial
        -----------------------------------------------------------------------
        Saida:
        -----------------------------------------------------------------------
        E - modulo de elasticidade atualizado
        -----------------------------------------------------------------------
        OBS:
        ***********************************************************************
        '''
        self.__u_var = u_var
        self.__t_var = t_var
        self.__E = modE
        self.__E0 = modE


    def update(self, n: float, u: float):
        '''
        **********************************************************************
        data criacao:     05/09/2020
        data modificacao: 00/00/0000
        ----------------------------------------------------------------------
        update : atualiza o modulo de elasticidade
        ----------------------------------------------------------------------
        Entrada:
        ----------------------------------------------------------------------
        n     - passo de tempo
        u     - deslocamento
         ----------------------------------------------------------------------
        Saida:
        -----------------------------------------------------------------------
        E - modulo de elasticidade atualizado
        -----------------------------------------------------------------------
        OBS:
        ***********************************************************************
        '''

        if self.__t_var:
            gt = 4*log10(10*n+10)
        else:
            gt = self.__E0

        if self.__u_var:
            ui = abs(u)
            fx = (3.176*ui**4 - 8.3625*ui**3 + 7.9172*ui**2 - 3.2258*ui + 0.9994)
        else:
            fx = 1.0

        self.__E = gt * fx
        return self.__E

    @property
    def E(self):
        '''
        retorna o modulo de elasticidade atualizado
        '''
        return self.__E