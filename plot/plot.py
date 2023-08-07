import matplotlib.pyplot as plt
import numpy as np


def main():
    """"""

    # ...
    solution = np.loadtxt("solution.csv", delimiter=",", dtype=float)
    assert len(solution) == len(solution[0])
    L = len(solution)

    # ...
    u = lambda x, y: solution[L - 1 - y, x]

    # ...
    xs = np.linspace(0, L - 1, L, dtype=int)
    ys = np.linspace(0, L - 1, L, dtype=int)
    X, Y = np.meshgrid(xs, ys)
    Z = u(X, Y)

    # ...
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u(x, y)')
    plt.show()


if __name__ == '__main__':
    main()
