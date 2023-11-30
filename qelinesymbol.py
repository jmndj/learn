def reverse_lines_in_range(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # 注意Python的索引是从0开始的，所以我们需要减1
    for i in range(1,73,2):
        lines[i*121:i*121+121] = reversed(lines[i*121:i*121+121])

    with open(output_file, 'w') as file:
        file.writelines(lines)

# 使用函数
def count_lines(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return len(lines)

# 使用函数
def remove_blank_lines(input_file, output_file):
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # 去除空行
    lines = [line for line in lines if line.strip()]

    with open(output_file, 'w') as file:
        file.writelines(lines)

# 使用函数
remove_blank_lines("C:/Users/dell/Desktop/paper/data/100010/5/bands.dat","C:/Users/dell/Desktop/paper/data/100010/5/bands.dat" )
#num_lines = count_lines('input.txt')
#print(f'The file has {num_lines} lines.')
reverse_lines_in_range("C:/Users/dell/Desktop/paper/data/100010/5/bands.dat", "C:/Users/dell/Desktop/paper/data/100010/5/bands.dat")
#print(count_lines("C:/Users/dell/Desktop/paper/data/none/bands.dat")/121)
