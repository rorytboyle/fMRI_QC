def check_func_conn_motion(motion, conn_mx_dir, conn_mx_filename_str):
    """
    Creates a dataframe and dict summarising the relationship between
    motion (mean framewise displacement) and functional connectivity at
    the global level, regional/node level, and edge level.

    :param motion: pandas dataframe with columns 'subid' and 'mean_FWD'. Number
    of rows = number of participants.
    :param conn_mx_dir: Path name for folder containing ONLY connectivity
    matrix files.
    :param conn_mx_filename_str: string containing name of connectivity matrix
    files listed in conn_mx_dir. Files should be named using pattern 
    subid_conn_mx_filenames e.g. 0001_connMx_zscored.csv
        
    :return motion_check_df: p * 4 dataframe with columns = 'level', 'label',
    'r', 'p' which contains correlation between functional connectivity and 
    mean FWD at global, regional, and edge-level.
    p = 1 + number of regions + number of edges in upper triangle of 
    connectivity matrix. e.g. for 268 * 268 connectivity matrix, p = 1 + 268
    + ((268*268) - 268)/2 = 36047.
    'level' = whether value refers to global, regional/node-level, or
    edge-level functional connectiivty correlation.
    'label' = label/indices of region/node and edge.
    'r' = Pearson's r for correlation between functional connectivity and mean
    FWD.
    'p' = p-value for Pearson's r between functional connectivity and mean FWD.
    :return summary: dict containing p-value for correlation between global 
    functional connectivity and mean FWD (global_p), Pearson's r between global
    functional connectivity and mean FWD (global_r), indices of edges sig
    correlated with mean FWD (ix_correlated_edges), indices of regions/nodes
    sig correlated with mean FW (ix_correlated_regions), total number and % of
    edges and regions sig correlated with mean FWD.

    Author: Rory Boyle rorytboyle@gmail.com
    Date: 25/11/2020
    """
    import pandas as pd
    import numpy as np
    import scipy.stats as stats
    from os import path

    #%% 1) Quick check of input
    # Take first ppt and check whether a connectivity matrix exists
    conn_mx = path.join(conn_mx_dir, str(motion['subid'][0]) + 
                        conn_mx_filename_str)
    if not path.exists(conn_mx):
        print("Check conn_mx_dir and conn_mx_filename_str as connectivity \nmatrix does not exist for this motion['subid'][0]")
        return
    
    # check motion dataframe supplied correctly
    motion_cols = ['subid', 'mean_FWD']
    if motion.columns.tolist() != motion_cols:
        print("Supply motion dataframe with columns 'subid' and mean_FWD'")
        return

    #%% 2) Collate all edges and regional functional connectivity for all ppts
    all_edges = []
    all_regional_avgs = []
    
    # loop through ppts and get connectivity matrix files and extract info
    for ppt in motion['subid'].values.tolist():
        
        # get name of connectivity matrix .csv file
        conn_mx_name = conn_mx_dir + '\\' + str(ppt) + conn_mx_filename_str
        
        # load connectivity matrix
        conn_mx = pd.read_csv(conn_mx_name, header=None)
        
        # extract edge indices and edge values of upper triangle in conn mx
        edge_ix, edge_values = extract_edges(conn_mx)
        
        # add 'subid' label and actual subid to lists before zipping so they can be
        # converted to a dict and added to list of dicts before being converted to
        # dataframe
        edge_ix.insert(0, 'subid')
        edge_values.insert(0, str(ppt))
        
        edges = dict(zip(edge_ix, edge_values))
        all_edges.append(edges)
        
        # get regional connectivity values
        ## set diagonal to nan before calculating mean for each node
        conn_mx.values[[np.arange(conn_mx.shape[0])]*2] = np.nan
        regional_avg = conn_mx.mean().values.tolist()
    
        node_ix = list(range(268))
            
        # add 'subid' label and actual subid to lists before zipping so they can be
        # converted to a dict and added to list of dicts before being converted to
        # dataframe
        node_ix.insert(0, 'subid')
        regional_avg.insert(0, str(ppt))
        
        
        regional_values = dict(zip(node_ix, regional_avg))
        all_regional_avgs.append(regional_values)
        
    edge_df = pd.DataFrame.from_records(all_edges).set_index('subid')
    regional_df = pd.DataFrame.from_records(all_regional_avgs).set_index('subid')
    
    #%% 3) create df to store correlation between head motion and functional
    #      connectivity at all 3 levels (global, regional/node, and edge)
    cols = ['level', 'label', 'r', 'p']
    motion_check_df = pd.DataFrame(columns=cols)
    
    #%% 4) Get correlation of global functional connectivity with motion
    
    # get global functional connectivity values for all ppts - calculate global
    # functional connectivity (average of functional connectivity in upper triangle
    # of connectivity matrix)
    global_func_conn = edge_df.mean(axis=1)
    
    # Calculate correlation at global level
    global_r, global_p = stats.pearsonr(motion['mean_FWD'], global_func_conn)
    
    # add to motion_check_df
    to_add = pd.Series(['global', 'global', global_r, global_p], index=cols)
    motion_check_df = motion_check_df.append(to_add, ignore_index=True)
    
    #%% 5) Get correlation of regional functional connectivity values with motion
    
    # make list to have a quick count of number sig correlated regions/nodes
    correlated_regions = []
    
    # Calculate correlation at regional level
    for node in range(len(conn_mx)):
        r, p = stats.pearsonr(motion['mean_FWD'], regional_df.iloc[:, node])
        
        # make quick list of nodes that are correlated with FWD
        if p < .05:
            correlated_regions.append(node)
        
        # add to motion_check_df
        to_add = pd.Series(['regional', node, r, p], index=cols)
        motion_check_df = motion_check_df.append(to_add, ignore_index=True)
    
    #%% 6) Get correlation between functional connectivity in each edge and motion
    
    # make list to have a quick count of number sig correlated edges
    correlated_edges = []
    
    # Calculate correlation at edge level
    for edge in range(edge_df.shape[1]):
        r, p = stats.pearsonr(motion['mean_FWD'], edge_df.iloc[:, edge])
        
        # make quick list of nodes that are correlated with FWD
        if p < .05:
            correlated_edges.append(edge)
            
        # add to motion_check_df
        to_add = pd.Series(['edge', edge, r, p], index=cols)       
        motion_check_df = motion_check_df.append(to_add, ignore_index=True)
        
    #%% 7) Summarise results
    ## Summarise correlation with global functional connectivity
    # get global correlation values
    global_corr = motion_check_df.loc[motion_check_df['level'] == 'global'
                                      ].iloc[:,2:].values.tolist()
    # unnest list
    global_corr = sum(global_corr, [])
    
    # add to summary dict
    summary = {'global_r': global_corr[0], 'global_p': global_corr[1]}
    
    ## Summarise correlation with regional functional connectivity
    
    regional_corr = motion_check_df.loc[motion_check_df['level'] == 'regional']
    total_regions = len(regional_corr.index)
    regional_corr = regional_corr[regional_corr['p'] < .05]
    
    # get number, percent, and indices of correlated regions
    summary['num_correlated_regions'] = len(regional_corr.index)
    
    summary['percent_correlated_regions']  = len(regional_corr.index
           )/total_regions * 100
           
    summary['ix_correlated_regions'] = regional_corr.loc[:, 'label'].values.tolist()
       
    ## Summarise correlation with edge-level functional connectivity
    
    edge_corr = motion_check_df.loc[motion_check_df['level'] == 'edge']
    total_edges = len(edge_corr.index)
    edge_corr = edge_corr[edge_corr['p'] < .05]
    
    # get number, percent, and indices of correlated regions
    summary['num_correlated_edges'] = len(edge_corr.index)
    
    summary['percent_correlated_edges']  = len(edge_corr.index
           )/total_edges * 100
           
    summary['ix_correlated_edges'] = edge_corr.loc[:, 'label'].values.tolist()
    
    #%% 8) return values
    return motion_check_df, summary
