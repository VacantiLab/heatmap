node_csv = '/Users/nate/Desktop/temporary/nodes.csv'
edge_csv = '/Users/nate/Desktop/temporary/edges.csv'
gdf_file = '/Users/nate/Desktop/temporary/graph.gdf'


with open(gdf_file,'w') as file_write:
    with open(node_csv,'r') as file_read:
        node_line = 1
        for line in file_read:
            if node_line == 1:
                file_write.write('nodedef>name VARCHAR,color VARCHAR,z DOUBLE\n')
            node_line = node_line + 1
            if node_line != 1:
                file_write.write(line)
    with open(edge_csv,'r') as file_read:
        edge_line = 1
        for line in file_read:
            if edge_line == 1:
                file_write.write('edgedef>node1 VARCHAR,node2 VARCHAR\n')
            edge_line = edge_line + 1
            if edge_line != 1:
                file_write.write(line)
