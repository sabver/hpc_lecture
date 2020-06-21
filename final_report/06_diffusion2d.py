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
###Ru_old through nt timesteps
def diffusion(n, u, u_old, surf):
    u_old = u.copy()
    u[1:-1, 1:-1] = (u_old[1:-1,1:-1] + 
                    nu * dt / dx**2 * (u_old[1:-1, 2:] - 2 * u_old[1:-1, 1:-1] + u_old[1:-1, 0:-2]) +
                    nu * dt / dy**2 * (u_old[2:,1: -1] - 2 * u_old[1:-1, 1:-1] + u_old[0:-2, 1:-1]))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    
    surf[0].remove()
    surf[0] = ax.plot_surface(X, Y, u[:], cmap=cm.seismic)

def my_diffusion(n, u, u_old, surf):
    u_old = u.copy()
    row, col = u.shape
    for j in range(1, row-1):
        for i in range(1, col-1):
            u[j, i] = readlist(n)[j][i] 

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    
    surf[0].remove()
    surf[0] = ax.plot_surface(X, Y, u[:], cmap=cm.seismic)

u_old = u
anim = animation.FuncAnimation(fig, diffusion, fargs=(u,u_old,surf), frames=nt, interval=50)
HTML(anim.to_html5_video())