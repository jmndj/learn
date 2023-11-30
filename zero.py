zero = []
with open("zero", "r") as f:
    line = f.readline()
    while line:
        zero.append(line.rstrip('\n'))
        line = f.readline()
for i in zero:
    if i == '':
        break
    new = float(i)*1000/1.602176634*10**19/(6.02214076*10**23)
    print(new)

