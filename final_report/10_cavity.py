import json
json_data_p = {}
json_data_u = {}
json_data_v = {}
def readjson(file_path):
    f = open(file_path)
    line = f.readline()
    # print(line)
    f.close()  
    array = json.loads(line) 
    return array

json_data_u = numpy.asarray(readjson("u_json.txt")) 
json_data_v = numpy.asarray(readjson("v_json.txt")) 
json_data_p = numpy.asarray(readjson("p_json.txt")) 

u = numpy.zeros((ny, nx))
v = numpy.zeros((ny, nx))
p = numpy.zeros((ny, nx))
b = numpy.zeros((ny, nx))
nt = 100
# u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu)
u = json_data_u
v = json_data_v
p = json_data_p