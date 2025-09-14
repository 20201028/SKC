import csv
def read_brightkite_ids(file_path):
    user_ids = list()
    with open(file_path, 'r', encoding='utf-8') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                user_id = parts[0]
                user_ids.append(user_id)
    return user_ids

def extract_dblp_data(dblp_file_path, user_ids):
    extracted_data = {}
    att = set()
    avg = 0
    with open(dblp_file_path, 'r', encoding='utf-8') as file:
        user_num = 0
        for line in file:
            if user_num >= len(user_ids):
                break
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                user_id = user_ids[user_num]
                extracted_data[user_id] = parts[1]
                att_list = parts[1].split(',')
                avg += len(att_list)
                att = att.union(set(att_list))
                user_num += 1
        print("att num:" + str(len(att)) + ", avg: " + str(avg/len(extracted_data)))
    return extracted_data

def write_extracted_data(extracted_data, output_file_path):
    with open(output_file_path, 'w', encoding='utf-8', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for user_id, attributes in extracted_data.items():
            writer.writerow([user_id, attributes])

def main(brightkite_file_path, dblp_file_path, output_file_path):
    user_ids = read_brightkite_ids(brightkite_file_path)
    extracted_data = extract_dblp_data(dblp_file_path, user_ids)
    write_extracted_data(extracted_data, output_file_path)

if __name__ == "__main__":
    # brightkite_file_path = 'Brightkite_processed_data.txt'  # 替换为你的Brightkite_processed_data.txt路径
    # dblp_file_path = 'pokec_attribute.txt'  # 替换为你的dblpAtt.txt路径
    # output_file_path = 'bk_extract_att.txt'  # 替换为你想要的输出文件路径
    brightkite_file_path = 'Gowalla_processed_data.txt'
    dblp_file_path = 'dblpAtt.txt'  # 替换为你的dblpAtt.txt路径
    output_file_path = 'go_extract_att.txt'  # 替换为你想要的输出文件路径
    main(brightkite_file_path, dblp_file_path, output_file_path)