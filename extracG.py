def read_edges_from_file(file_path):
    edges = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                node1, node2 = map(int, line.strip().split())
                edges.append((node1, node2))
    return edges

def get_top_thrd_percent_nodes(edges, thrd):
    all_nodes = set()
    for node1, node2 in edges:
        all_nodes.add(node1)
        all_nodes.add(node2)
    sorted_nodes = sorted(all_nodes)
    top_thrd_percent_count = int(len(sorted_nodes) * thrd)
    top_thrd_percent_nodes = set(sorted_nodes[:top_thrd_percent_count])
    return top_thrd_percent_nodes

def filter_edges_with_top_thrd_percent_nodes(edges, top_thrd_percent_nodes):
    filtered_edges = []
    for node1, node2 in edges:
        if node1 in top_thrd_percent_nodes and node2 in top_thrd_percent_nodes:
            filtered_edges.append((node1, node2))
    return filtered_edges

# 假设文件名为 'data.txt'
dataName = "Brightkite_edges"
file_path = dataName + '.txt'
edges = read_edges_from_file(file_path)

thresholds = [0.2, 0.4, 0.6, 0.8]

for thrd in thresholds:
    top_thrd_percent_nodes = get_top_thrd_percent_nodes(edges, thrd)
    filtered_edges = filter_edges_with_top_thrd_percent_nodes(edges, top_thrd_percent_nodes)
    
    # 输出结果到文件
    output_file_path = f'{dataName}_{int(thrd*10)}.txt'
    with open(output_file_path, 'w') as output_file:
        for edge in filtered_edges:
            output_file.write(f"{edge[0]}\t{edge[1]}\n")