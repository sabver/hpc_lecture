import json
json_data = {}
def readjson(file_path):
    f = open(file_path)
    line = f.readline()
    print(line)
    f.close()  
    array = json.loads(line) 
    return array

def readlist(index):
    return json_data[str(index)]


json_data = readjson("json.txt")

for i in range(0, nt):
    print(readlist(i))



