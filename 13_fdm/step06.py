import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

nx = 41
ny = 41
dx = 2./(nx-1)
dy = 2./(ny-1)
nt = 50
dt = .01
x = numpy.linspace(0,2,nx)
y = numpy.linspace(0,2,ny)
X, Y = numpy.meshgrid(x,y) 
u = numpy.ones((ny,nx))
v = numpy.ones((ny,nx))
u[10:20, 10:20] = 2
v[10:20, 10:20] = 2

fig = pyplot.figure(figsize=(11,7), dpi=100)

for n in range(nt):
    un = u.copy()
    vn = v.copy()
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            u[j,i] = un[j,i]-un[j,i]*dt/dx*(un[j,i]-un[j,i-1])-vn[j,i]*dt/dy*(un[j,i]-un[j-1,i])
            v[j,i] = vn[j,i]-un[j,i]*dt/dx*(vn[j,i]-vn[j,i-1])-vn[j,i]*dt/dy*(vn[j,i]-vn[j-1,i])
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.coolwarm)
    ax.set_zlim3d(1, 2) 
    pyplot.pause(0.05)
    pyplot.clf()
pyplot.show()
