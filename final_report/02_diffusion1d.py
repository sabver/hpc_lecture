import numpy
from matplotlib import pyplot
from matplotlib import animation, rc
from IPython.display import HTML
%matplotlib inline

import json
json_data = {}
def readjson(file_path):
    f = open(file_path)
    line = f.readline()
    # print(line)
    f.close()  
    array = json.loads(line) 
    return array

def readlist(index):
    return json_data[str(index)]

json_data = readjson("json.txt")

nx = 41
dx = 2 / (nx - 1)
nt = 20    #the number of timesteps we want to calculate
nu = 0.3   #the value of viscosity
sigma = .2 #sigma is a parameter, we'll learn more about it later
dt = sigma * dx**2 / nu #dt is defined using sigma ... more later!

x = numpy.linspace(0, 2, nx)
u = numpy.ones(nx)      #a numpy array with nx elements all equal to 1.
u[int(.5 / dx):int(1 / dx + 1)] = 2  #setting u = 2 between 0.5 and 1 as per our I.C.s
fig, ax = pyplot.subplots()
line, = ax.plot(x, u)

def my_diffusion(n,u,u_old):
    u = readlist(n)
    line.set_data(x, u)

def diffusion(n, u, u_old):
    u_old = u.copy()
    for i in range(1, nx-1):
        u[i] = u_old[i] + nu * dt / dx**2 * (u_old[i+1] - 2 * u_old[i] + u_old[i-1])
        line.set_data(x, u)

u_old = u
anim = animation.FuncAnimation(fig, my_diffusion, fargs=(u,u_old), frames=nt)
HTML(anim.to_html5_video())