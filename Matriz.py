class Matriz():

    def __init__(self, modE, L: float = 1.0, A: float = 1.0):

        self.__L = L
        self.__A = A
        self.__E = modE
        self.__AL = L*A


    def update(self, n: float, u: float):
        '''
        **********************************************************************
        data criacao:     05/09/2020
        data modificacao: 00/00/0000
        ----------------------------------------------------------------------
        update : atualiza a matriz de coeficientes
        ----------------------------------------------------------------------
        Entrada:
        ----------------------------------------------------------------------
        n     - passo de tempo
        u     - deslocamento
         ----------------------------------------------------------------------
        Saida:
        -----------------------------------------------------------------------
        K - atualizado
        -----------------------------------------------------------------------
        OBS:
        ***********************************************************************
        '''

        E = self.__E.update(n, u)
        return E*self.__AL

    def get_E(self):
        '''
        retorna o modulo de elasticidade atualizado
        '''
        return self.__E.E

    @property
    def L(self):
        '''
        retorna a comprimento
        '''
        return self.__L

    @property
    def A(self):
        '''
        retorna a area
        '''
        return self.__A

    @property
    def E(self):
        '''
        retorna o modulo de elasticidade atualizado
        '''
        return self.get_E()
