import random
import time
import os

def generate_scaled_graph_simple(original_file: str, num_copies: int, output_file: str, connection_prob: float = 0.001, node_num: int = 1000):
    """
    简化版的图扩展生成器
    """
    print(f"生成 {num_copies} 倍扩展图...")
    
    # # 分析原始图
    # node_set = set()
    # with open(original_file, 'r') as f:
    #     for line in f:
    #         line = line.strip()
    #         if line and not line.startswith('#'):
    #             parts = line.split()
    #             if len(parts) >= 2:
    #                 node_set.add(int(parts[0]))
    #                 node_set.add(int(parts[1]))
    
    # original_nodes = len(node_set)
    # print(f"原始图节点数: {original_nodes}")
    
    start_time = time.time()
    connections_made = 0
    
    with open(output_file, 'w') as out_f:
        # 写入元数据
        # out_f.write(f"#vNum: {original_nodes * num_copies}\n")
        # out_f.write(f"#eNum: estimated\n")  # 边数难以精确预估
        # out_f.write(f"#copies: {num_copies}\n")
        
        # 复制所有副本的内部边
        for copy_id in range(num_copies):
            offset = copy_id * node_num
            print(f"处理副本 {copy_id + 1}/{num_copies}")
            
            with open(original_file, 'r') as in_f:
                for line in in_f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            u = int(parts[0]) + offset
                            v = int(parts[1]) + offset
                            out_f.write(f"{u}\t{v}\n")
        
        # 添加连接边
        print("添加连接边...")
        for i in range(num_copies):
            for j in range(i + 1, num_copies):
                offset_i = i * node_num
                offset_j = j * node_num
                
                for _ in range(int(node_num * connection_prob)):
                    if random.random() < connection_prob:
                        node_i = random.randint(0, node_num-1) + offset_i
                        node_j = random.randint(0, node_num-1) + offset_j
                        out_f.write(f"{node_i}\t{node_j}\n")
                        connections_made += 1
    
    generation_time = time.time() - start_time
    file_size = os.path.getsize(output_file) / (1024 * 1024 * 1024)
    
    print(f"生成完成!")
    print(f"连接边数量: {connections_made}")
    print(f"文件大小: {file_size:.2f} GB")
    print(f"生成时间: {generation_time:.2f} 秒")

# 使用示例
if __name__ == "__main__":
    # 生成多个规模的图
    scales = [20, 40, 60, 80, 100]
    # scales = [2]
    
    for scale in scales:
        output_filename = f"graph_scale_{scale}x.txt"
        generate_scaled_graph_simple(
            original_file="dblp_graph.txt",
            num_copies=scale,
            output_file=output_filename,
            connection_prob=0.001,
            node_num = 1278674
        )
    #     print()
    # scale = 20
    # output_filename = f"graph_scale_{scale}x.txt"
    # generate_scaled_graph_simple(
    #     original_file="dblp_graph.txt",
    #     num_copies=scale,
    #     output_file=output_filename,
    #     connection_prob=0.001,
    #     node_num = 1278673
    # )