u = numpy.ones((ny, nx))
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2
# print(u)
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

# print(readlist(0)[0][0])

# def my_convection(n, u, u_old, surf):
#     u_old = u.copy()
#     row, col = u.shape
#     for j in range(1, row):
#         for i in range(1, col):
#             u[j, i] = readlist(n)[j][i]
#             u[0, :] = 1
#             u[-1, :] = 1
#             u[:, 0] = 1
#             u[:, -1] = 1    
#     # u = readlist(n)
#     surf[0].remove()
#     surf[0] = ax.plot_surface(X, Y, u[:], cmap=cm.seismic)

def convection(n, u, u_old, surf):
    u_old = u.copy()
    row, col = u.shape
    for j in range(1, row):
        for i in range(1, col):
            u[j, i] = (u_old[j, i] - (c * dt / dx * (u_old[j, i] - u_old[j, i - 1])) -
                                  (c * dt / dy * (u_old[j, i] - u_old[j - 1, i])))
            u[0, :] = 1
            u[-1, :] = 1
            u[:, 0] = 1
            u[:, -1] = 1
    surf[0].remove()
    surf[0] = ax.plot_surface(X, Y, u[:], cmap=cm.seismic)

u_old = u
anim = animation.FuncAnimation(fig, convection, fargs=(u,u_old,surf), frames=nt, interval=50)
HTML(anim.to_html5_video())