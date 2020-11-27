def extract_edges(conn_mx):
    """
    Extracts edge strength values and their index in the connectivity matrix

    :param conn_max: connectivity matrix (symmetrical e.g. 268 * 268)
    :return edge_ix: list containg edge strength value of each edge in upper 
    triangle
    :return edge_values: list containing tuples with index of each edge in 
    upper triangle of matrix
    Author: Rory Boyle rorytboyle@gmail.com
    Date: 25/11/2020
    """
    import numpy as np
        
    # create empty lists to store variables
    edge_ix = []
    edge_values = []

    # Extract upper triangular part of connectivity matrix
    upper_triangle = np.triu(conn_mx, k=1)

    # Loop through upper triangle of connectivity matrix and extract edges and
    # their indices
    
    for i in range(len(conn_mx)):
        for j in range(len(conn_mx)):
            if upper_triangle[i, j] != 0:
                
                edge_ix.append((i, j))
                edge_values.append(upper_triangle[i, j])
                             
    return edge_ix, edge_values