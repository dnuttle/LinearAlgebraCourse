import matplotlib.pyplot as plt
def plotv2d(v, axlims=[-4, 4, -4, 4]):
    plt.plot([0, v[0]], [0, v[1]])
    plt.axis('equal')
    plt.plot([axlims[0], axlims[1]],[0, 0],'k--')
    plt.plot([0,0], [axlims[0], axlims[1]],'k--')
    plt.grid()
    plt.axis((axlims))
    plt.show()