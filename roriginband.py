#!/usr/bin/python
# Read the input data from a file
with open("C:/Users/dell/Desktop/fourth grade/Cmca_12/strain/tensile/010/BAND.dat", "r") as file:
    lines = file.readlines()

current_band = None
path = {}
energy = {}
output = []

for line in lines:
    if line.startswith("# Band-Index"):
        current_band = int(line.split()[-1])
        energy[current_band] = []
        path[current_band] = []
    elif line.strip() and current_band is not None:
        values = line.split()
        energy[current_band].append(float(values[-1]))
        path[current_band].append(float(values[0]))

out_path_list = path[1][:]


for key, value in path.items():
    if value[0] > value[1]:
        energy[key].reverse()

for i in range(len(out_path_list)):
    temp = []
    temp.append(out_path_list[i])
    for key, value in energy.items():
        temp.append(value[i])
    
    output.append(" ".join([f"{values}" for values in temp]))   

with open("C:/Users/dell/Desktop/fourth grade/Cmca_12/strain/tensile/010/m.dat", "w") as file:
    file.write("\n".join(output))