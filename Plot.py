import matplotlib.pyplot as plt

def plot_res(res_inc: list, res_n_inc: list
             , res_linear: float, res_inc_cons: float
             , L: float):

    x = [i for i in range(len(res_inc))]
    y_inc = [L  + x for x in res_inc]
    y_n_inc = [L  + x for x in res_n_inc]
    y_linear= [L  + x for x in res_linear]
    y_inc_cons= [L  + x for x in res_inc_cons]

    fig, ax = plt.subplots()

    # plot res
    ax.plot(x, y_linear, label = 'Lin', ls='--', marker = '^')
    ax.plot(x, y_n_inc, label = 'NR', ls='--', marker='*')
    ax.plot(x, y_inc, label = 'NRI', ls = '--', marker='s')
    ax.plot(x, y_inc_cons, label = 'NRCI', ls='--', marker = 'v')
    ax.hlines(1,0, 40, ls='--', colors='black')
    ax.vlines(10,-10, 10, ls='--', colors='red')
    ax.vlines(20,-10, 10, ls='--', colors='red')
    ax.legend(loc = 'lower right')
    ax.set_ylabel('Comprimento')
    ax.set_xlabel('Tempo')
    ax.set_ylim((0.4, L*1.2))
    ax.set_xlim((0, len(x)))
#   plt.show()
    plt.savefig('desloc.png', dpi = 300)

def plot_f(f: list):

    x = [i for i in range(len(f))]

    fig, ax = plt.subplots()

    ax.plot(x, f, ls = '-')
    ax.set_ylabel('Carga')
    ax.set_xlabel('Tempo')
    ax.set_ylim((min(f)-0.1,max(f)+0.1))
    ax.set_xlim((0, 40))
#   plt.show()
    plt.savefig('forces.png', dpi = 300)

def plot_E(E: list):

    E_inc, E_n_inc, E_n_inc_cons, E_linear = E

    x = [i for i in range(len(E[0]))]

    fig, ax = plt.subplots()

    ax.plot(x, E_linear, label = 'Lin', ls = '--', marker='^')
    ax.plot(x, E_n_inc, label = 'NR', ls = '--', marker='*')
    ax.plot(x, E_inc, label = 'NRI', ls = '--', marker='s')
    ax.plot(x, E_n_inc_cons, label = 'NRCI', ls = '--', marker='v')
    ax.vlines(10,0, 50, ls='--', colors='red')
    ax.vlines(20,0, 50, ls='--', colors='red')
    ax.set_ylabel('Modulo de Elasticidade')
    ax.set_xlabel('Tempo')
    ax.legend(loc = 'lower right')
    ax.set_ylim((min(E_inc)-0.1,max(E_inc)+0.1))
    ax.set_xlim((0, len(x)))
#   plt.show()
    plt.savefig('modE.png', dpi = 300)
