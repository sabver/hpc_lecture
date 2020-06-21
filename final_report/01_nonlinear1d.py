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

# for i in range(0, nt):
#     print(readlist(i))

def convection(n, u, u_old):
    u_old = u.copy()
    for i in range(1, nx):
        u[i] = u_old[i] - u_old[i] * dt / dx * (u_old[i] - u_old[i-1])
        line.set_data(x, u)

def my_convection(n,u,u_old):
    u = readlist(n)
    line.set_data(x, u)

u_old = u
anim = animation.FuncAnimation(fig, my_convection, fargs=(u,u_old), frames=nt)
HTML(anim.to_html5_video())