# 假设数据文件名为data.txt，第一列是名称，第二列是数据
filename = "info_big_bac_list"

# 用字典来保存已经打开的文件对象，避免频繁打开和关闭文件
file_dict = {}

with open(filename, "r") as f:
    for line in f:
        # 使用split方法分割每一行
        parts = line.strip().split(" ")
        name = parts[1]
        data = "\t".join(parts)

        # 获取前四个字符
        prefix = name[:4]

        # 如果该前缀的文件还没有打开，则打开一个新文件并保存到字典中
        if prefix not in file_dict:
            new_file = open(f"each_seq/{prefix}.txt", "w")
            file_dict[prefix] = new_file
        # 否则从字典中获取已经打开的文件对象
        else:
            new_file = file_dict[prefix]

        # 写入数据到文件中
        new_file.write(f"{data}\n")

# 关闭所有文件
for prefix in file_dict:
    file_dict[prefix].close()
