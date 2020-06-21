import json
json_data_u = {}
json_data_v = {}
def readjson(file_path):
    f = open(file_path)
    line = f.readline()
    # print(line)
    f.close()  
    array = json.loads(line) 
    return array

def readlistu(index):
    return json_data_u[str(index)]

def readlistv(index):
    return json_data_v[str(index)]

json_data_u = readjson("u_json.txt")
json_data_v = readjson("v_json.txt")

def burgers(n, u, u_old, v, v_old, surf):
    u_old = u.copy()
    v_old = v.copy()

    u[1:-1, 1:-1] = (u_old[1:-1, 1:-1] -
                     dt / dx * u_old[1:-1, 1:-1] * 
                     (u_old[1:-1, 1:-1] - u_old[1:-1, 0:-2]) - 
                     dt / dy * v_old[1:-1, 1:-1] * 
                     (u_old[1:-1, 1:-1] - u_old[0:-2, 1:-1]) + 
                     nu * dt / dx**2 * 
                     (u_old[1:-1,2:] - 2 * u_old[1:-1, 1:-1] + u_old[1:-1, 0:-2]) + 
                     nu * dt / dy**2 * 
                     (u_old[2:, 1:-1] - 2 * u_old[1:-1, 1:-1] + u_old[0:-2, 1:-1]))
    
    v[1:-1, 1:-1] = (v_old[1:-1, 1:-1] - 
                     dt / dx * u_old[1:-1, 1:-1] *
                     (v_old[1:-1, 1:-1] - v_old[1:-1, 0:-2]) -
                     dt / dy * v_old[1:-1, 1:-1] * 
                    (v_old[1:-1, 1:-1] - v_old[0:-2, 1:-1]) + 
                     nu * dt / dx**2 * 
                     (v_old[1:-1, 2:] - 2 * v_old[1:-1, 1:-1] + v_old[1:-1, 0:-2]) +
                     nu * dt / dy**2 *
                     (v_old[2:, 1:-1] - 2 * v_old[1:-1, 1:-1] + v_old[0:-2, 1:-1]))
     
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    
    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1
    
    surf[0].remove()
    surf[0] = ax.plot_surface(X, Y, u[:], cmap=cm.seismic)

def my_burgers(n, u, u_old, v, v_old, surf):
    u_old = u.copy()
    v_old = v.copy()
    row, col = u.shape
    for j in range(1, row-1):
        for i in range(1, col-1):
            u[j, i] = readlistu(n)[j][i]
            v[j, i] = readlistv(n)[j][i] 
     
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1
    
    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1
    
    surf[0].remove()
    surf[0] = ax.plot_surface(X, Y, u[:], cmap=cm.seismic)

u_old = u
v_old = v
anim = animation.FuncAnimation(fig, my_burgers, fargs=(u,u_old,v,v_old,surf), frames=nt, interval=50)
HTML(anim.to_html5_video())