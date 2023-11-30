kn = 0  #横坐标k的数目
kpath = []
p = []
hang =[]
a = [1,2,3,4]
with open("C:/Users/dell/Desktop/paper/data/100010/5/freq.plot","r") as file:
    lines = file.readlines()
    with open("C:/Users/dell/Desktop/paper/data/100010/5/freqs.plot","w") as f:
        for line in lines:
            k = line.split()
            if len(k) == 2:
                kn = kn+1      
            else:
                print(kn)
                break    #确定横坐标K的数目
        for line in lines:
            k = line.split()
            if len(k) == 2:
                kpath.append(k[0])
                p.append(k[1])   # 将横坐标K以及声子频率全部获取
        for i in range(0,kn):
            hang.append(kpath[i])
            for j in range(0,len(kpath)):
                if kpath[j] == kpath[i]:
                    hang.append(float(p[j])/33.3)
            for i in hang:
                print(float(i),'',end='',file=f)
            print('',file= f)
            hang =[]          


