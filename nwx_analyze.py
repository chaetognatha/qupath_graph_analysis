#!/usr/bin/env python
# coding: utf-8
# author: Mattis
# email: mattisknulst@gmail.com
import argparse
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from libpysal import weights
from alive_progress import alive_bar
import multiprocessing as mp
import warnings

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--file', type=str, help='path to the csv file',
                    required=True)
parser.add_argument('-o', '--results_dir',
                    type=str,
                    help='path to the results directory'
                         ', default current directory',
                    default='.')
parser.add_argument('-c', '--critical_distance',
                    type=int, help='critical distance, default=30',
                    default=30)
parser.add_argument('-p', '--pair', type=str,
                    help='cell type pair', required=True)
parser.add_argument('-s', '--sep', help='csv separator, default tab',
                    default='\t')
parser.add_argument('-d', '--decimal', help='float decimal sign, default .',
                    default='.')
parser.add_argument('-t', '--tiff_dir', type=str,
                    help='path to the tiff directory', required=True)
parser.add_argument('-n', '--n_workers', type=int, default=1,
                    help="number of processes to use")




def cluster_cooccurrence(df, G, critical_distance, iterations=5):
    """
    This function calculates the cluster cooccurence, by creating a randomly
    distributed coordinate system with cell types using the max and min values
    from the corresponding columns in the data. Then it counts the average
    number of edges over a set of critical distances and divides that with
    with the same for the random case.
    :param df: dataframe
    returns: cluster cooccurence value
    """
    # get the max and min values for x and y coordinates
    x_min = df.loc[:,'X'].min()
    x_max = df.loc[:,'X'].max()
    y_min = df.loc[:,'Y'].min()
    y_max = df.loc[:,'Y'].max()
    # get number of rows
    n_rows = df.shape[0]
    results = {
        'between_all_ccr': [],
        'within_1_ccr': [],
        'within_2_ccr': [],
    }
    for i in range(iterations):
        # create a random coordinate system
        random_x = np.random.randint(x_min, x_max, size=(n_rows, 1))
        random_y = np.random.randint(y_min, y_max, size=(n_rows, 1))
        # create a dataframe with the random coordinates and name columns
        random_coordinates = pd.DataFrame(
            np.concatenate((random_x, random_y), axis=1),
            columns=["X", "Y"])
        # add the cell types to the random coordinates
        random_coordinates["Class"] = df["Class"]
        # create a graph from the random coordinates
        random_graph, random_weights = make_graph(random_coordinates,
                                                  critical_distance)

        # only count edges between cells of different types
        edges = [edge for edge in G.edges if G.nodes[edge[0]]['cell_type'] !=
                 G.nodes[edge[1]]['cell_type']]
        # same thing for random graph
        random_edges = [edge for edge in random_graph.edges if
                        random_graph.nodes[edge[0]]['cell_type'] !=
                        random_graph.nodes[edge[1]]['cell_type']]
        # calculate the cluster cooccurence
        if len(random_edges) > 0:
            cluster_cooccurence = len(edges) / len(random_edges)
        else:
            cluster_cooccurence = np.nan
        # pick out cell type 1
        cell_type_1 = df.Class.unique()[0]
        # pick out cell type 2
        cell_type_2 = df.Class.unique()[1]
        # count within group connections of cell type 1
        within_group_1 = [edge for edge in G.edges if
                            G.nodes[edge[0]]['cell_type'] == cell_type_1 and
                            G.nodes[edge[1]]['cell_type'] == cell_type_1]
        # count within group connections of cell type 2
        within_group_2 = [edge for edge in G.edges if
                            G.nodes[edge[0]]['cell_type'] == cell_type_2 and
                            G.nodes[edge[1]]['cell_type'] == cell_type_2]

        # count within group connections of cell type 1 in random graph
        random_within_group_1 = [edge for edge in random_graph.edges if
                                    random_graph.nodes[edge[0]]['cell_type'] ==
                                    cell_type_1 and
                                    random_graph.nodes[edge[1]]['cell_type'] ==
                                    cell_type_1]
        # count within group connections of cell type 2 in random graph
        random_within_group_2 = [edge for edge in random_graph.edges if
                                    random_graph.nodes[edge[0]]['cell_type'] ==
                                    cell_type_2 and
                                    random_graph.nodes[edge[1]]['cell_type'] ==
                                    cell_type_2]

        # calculate the cluster cooccurence for within group connections of
        # cell type 1
        if len(random_within_group_1) > 0:
            within_group_1_cooccurence = len(within_group_1) / \
                                        len(random_within_group_1)
        else:
            within_group_1_cooccurence = np.nan
        # calculate the cluster cooccurence for within group connections of
        # cell type 2
        if len(random_within_group_2) > 0:
            within_group_2_cooccurence = len(within_group_2) / \
                                        len(random_within_group_2)
        else:
            within_group_2_cooccurence = np.nan



        results['between_all_ccr'].append(cluster_cooccurence)
        results['within_1_ccr'].append(within_group_1_cooccurence)
        results['within_2_ccr'].append(within_group_2_cooccurence)

    # get average ccr
    results['between_all_ccr'] = np.nanmean(results['between_all_ccr'])
    results['within_1_ccr'] = np.nanmean(results['within_1_ccr'])
    results['within_2_ccr'] = np.nanmean(results['within_2_ccr'])

    return results


def make_graph(df, critical_distance):
    """
    Make a graph from a dataframe with x, y coordinates and cell type
    :param df: dataframe with x, y coordinates and cell type (qupath output)
    :param critical_distance: distance in pixels
    :return: graph and weights
    """
    # extract the spatial coordinates
    coordinates = df.loc[:,'X':'Y']
    # extract the cell types
    cell_types = df.Class
    # Creating a graph from coordinates
    positions = coordinates.to_numpy()
    # create a weights object
    # catch user warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = weights.DistanceBand.from_array(positions,
                                        threshold=critical_distance)

    # create a networkx graph from the weights object
    G = nx.from_numpy_matrix(w.full()[0])

    # add the coordinates to the graph
    for i, (x, y) in enumerate(positions):
        G.nodes[i]['pos'] = (x, y)

    # add the cell types to the graph
    for i, cell_type in enumerate(cell_types):
        G.nodes[i]['cell_type'] = cell_type

    return G, w


def plot_graph(G, pair, results_dir, cell_type_filter, image_file,
               critical_distance, prepend=''):
    """
    This function plots the graph
    :param G: graph
    :param pair: cell type pair
    :param results_dir: path to the results directory
    :param cell_type_filter: cell type filter
    :param image_file: path to the image file
    :param critical_distance: critical distance
    :param prepend: string to prepend to the filename
    :return: None
    """
    # get image
    my_img = plt.imread(image_file)
    # we don't want colons in file names
    nice = pair.replace(':', '_')
    out_file = f'{nice}_{prepend}image_{image}_distance_{critical_distance}.png'
    cell_types = [cell_type for node, cell_type in G.nodes(data='cell_type')]
    # plot the graph
    fig = plt.figure(figsize=(10, 10))
    y_lim, x_lim = my_img.shape[0], my_img.shape[1]
    extent = 0, x_lim, 0, y_lim
    # transpose image

    my_img = np.flipud(my_img)
    # draw the image
    plt.imshow(my_img, cmap='gray', extent=extent, interpolation='nearest')
    # cell type name and colors
    c1 = f"{cell_type_filter[0]} (lime)"
    c2 = f"{cell_type_filter[1]} (magenta)"
    color_map = ["lime" if cell_type == cell_type_filter[0] else "magenta" for
                 cell_type in
                 cell_types]
    # draw graph on image
    nx.draw_networkx(G, pos=nx.get_node_attributes(G, 'pos'),
                     node_color=color_map, node_size=10, alpha=0.5,
                     edge_color='yellow', width=0.5, with_labels=False)

    plt.title(f'{c1} {c2} image {image} distance {critical_distance}px')

    plt.savefig(f'{results_dir}/{out_file}')
    plt.close(fig)


def calculate_statistics(df, G, w, cell_type_filter, critical_distance):
    """
    This function calculates the statistics
    :param df: dataframe with x, y coordinates and cell type (qupath output)
    :param G: graph
    :param w: weights object
    :param cell_type_filter: cell type filter
    :param critical_distance: critical distance needed for ccr
    :return: tuple of statistics
    """
    # test that df contains class
    assert 'Class' in df.columns, "df does not contain Class column"

    aac = nx.attribute_assortativity_coefficient(G, 'cell_type')
    # count nr of cells in each class
    n_cells_class_1 = df.loc[
        df['Class'] == cell_type_filter[0], 'Class'].count()
    n_cells_class_2 = df.loc[
        df['Class'] == cell_type_filter[1], 'Class'].count()
    # count number of islands
    n_islands = len(w.islands)
    # get all nodes belonging to class 1
    class_1_nodes = [node for node, cell_type in G.nodes(data='cell_type') if
                     cell_type == cell_type_filter[0]]
    # calculate group degree centrality for class 1
    try:
        centrality_measures = nx.group_degree_centrality(G, class_1_nodes)
    except Exception as e:
        centrality_measures = 'NA'
    # ratio is the proportion of class 1 cells in the network
    ratio = n_cells_class_1 / (n_cells_class_2 + n_cells_class_1)
    ccr = cluster_cooccurrence(df, G, critical_distance)
    ccr_between_1_2 = ccr['between_all_ccr']
    ccr_within_1 = ccr['within_1_ccr']
    ccr_within_2 = ccr['within_2_ccr']

    results =  [aac, n_cells_class_1, n_cells_class_2,
           n_islands, centrality_measures, ratio, ccr_between_1_2,
           ccr_within_1, ccr_within_2]
    return results


def network_plot(df, image, tiff_dir, critical_distance, results_dir,
                 prepend='',
                 pair='CD45:PANCK'):
    """
    inner main function for iterating image by image
    """
    cell_type_filter = pair.split(':')
    # class_1 = df.groupby('Class').get_group(cell_type_filter[0])
    # class_2 = df.groupby('Class').get_group(cell_type_filter[1])
    # get image file path
    image_file = tiff_dir + '/' + image
    # filter on cell type included in analysis
    filtered = df[(df['Class'] == cell_type_filter[0]) | (
            df['Class'] == cell_type_filter[1])].reset_index(drop=True)
    # extract the image
    one_pic = filtered.groupby('Image').get_group(image).reset_index(drop=True)
    # make graph
    G, w = make_graph(one_pic, critical_distance)
    # Plotting the graph
    plot_graph(G, pair, results_dir, cell_type_filter, image_file,
               critical_distance, prepend=prepend)
    # create results list
    results = [image, cell_type_filter[0], cell_type_filter[1]]
    # calculate statistics
    results.extend(calculate_statistics(one_pic, G, w,
                                   cell_type_filter, critical_distance))

    return results




if __name__ == "__main__":
    args = parser.parse_args()
    # create the output dictionary
    out = {'image':[],
           'cell_type_1':[],
           'cell_type_2':[],
           'aac':[],
           'n_cells_class_1':[],
           'n_cells_class_2':[],
           'n_islands':[],
           'centrality_measures':[],
           'ratio':[],
           'ccr_between_1_2':[],
           'ccr_within_1':[],
           'ccr_within_2':[]}
    # read the csv file
    df = pd.read_csv(args.file,
                     sep=args.sep,
                     decimal=args.decimal,
                     low_memory=False,
                     skiprows=1,
                     usecols=[0, 1, 2, 3, 4, 5, 6],
                     names=['Image',
                            'Name',
                            'Class',
                            'Parent',
                            'ROI', 'X', 'Y']).dropna(axis=1)
    # list all images
    images = df['Image'].unique()
    # check that all images have both classes
    for image in images:
        if args.pair.split(':')[0] not in df.groupby('Image').get_group(
            image)['Class'].unique():
            print(f"{args.pair.split(':')[0]} not in {image}")
        if args.pair.split(':')[1] not in df.groupby('Image').get_group(
            image)['Class'].unique():
            print(f"{args.pair.split(':')[1]} not in {image}")

    with alive_bar(len(images)*2) as bar:
        # iterate over images
        for image in images:
            try:
            # run the network plot
                results = network_plot(df, image, args.tiff_dir,
                                     args.critical_distance, args.results_dir,
                                        pair=args.pair)
            except Exception as e:
                bar.text(f'Error in {image}: {e}')
            bar()
            # append to output dictionary
            for key, value in zip(out.keys(), results):
                out[key].append(value)
            bar()


    # create dataframe from out dictionary
    out_df = pd.DataFrame(out)
    print(out_df.head())

    # save the output
    nice = args.pair.replace(':', '_')
    out_df.to_csv(
        f'{args.results_dir}/results_{nice}_distance_{args.critical_distance}.csv',
        index=False)
