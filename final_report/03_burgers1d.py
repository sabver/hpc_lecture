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

def burgers(n, u, u_old):
    u_old = u.copy()
    for i in range(1, nx-1):
        u[i] = u_old[i] - u_old[i] * dt / dx *(u_old[i] - u_old[i-1]) + nu * dt / dx**2 *\
                (u_old[i+1] - 2 * u_old[i] + u_old[i-1])
    u[0] = u_old[0] - u_old[0] * dt / dx * (u_old[0] - u_old[-2]) + nu * dt / dx**2 *\
                (u_old[1] - 2 * u_old[0] + u_old[-2])
    u[-1] = u[0]
    u_analytical = numpy.asarray([ufunc(n * dt, xi, nu) for xi in x])
    line.set_data(x, u)
    line2.set_data(x, u_analytical)

def my_burgers(n,u,u_old):
    u = readlist(n)
    u_analytical = numpy.asarray([ufunc(n * dt, xi, nu) for xi in x])
    line.set_data(x, u)
    line2.set_data(x, u_analytical)

u_old = u
anim = animation.FuncAnimation(fig, my_burgers, fargs=(u,u_old), frames=nt, interval=50)
HTML(anim.to_html5_video())