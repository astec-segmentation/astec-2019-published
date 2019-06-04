import pickle as pkl
import numpy as np
import os


def timeNamed(filename, time):
    time_point = ('00' + str(time))[-3:]  # Format time on 3 digit
    return filename.replace('$TIME', time_point)


def timesNamed(filename, strtime1, time1, strtime2, time2):
    time_point1 = ('00' + str(time1))[-3:]  # Format time on 3 digit
    time_point2 = ('00' + str(time2))[-3:]  # Format time on 3 digit
    return filename.replace(strtime1, time_point1).replace(strtime2, time_point2)


def format_names(names_which_matter):
    """
    Return an ensured formated cell names
    """
    tmp = {}
    for k, v in names_which_matter.iteritems():
        val = v.split('.')[1][:-1]
        tmp[k] = v.split('.')[0] + '.' + ('00' + val)[-3:] + v[-1]
    return tmp


def write_tlp_from_lin_tree(name, lin_tree_information, ListProperties=[]):
    """
    Write a lineage tree into an understable tulip file
    name : path to the tulip file to create
    lin_tree : lineage tree to write
    properties : dictionary of properties { 'Property name': [{c_id: prop_val}, default_val]}
    """

    lin_tree = lin_tree_information['lin_tree']
    if ListProperties is None:
        ListProperties = ['volumes_information', 'h_mins_information', 'sigmas_information']
    properties = []  # Properties Dictionary initialisation
    for l in ListProperties:
        if l is not 'Names':
            properties.append((l, lin_tree_information[l], len(properties) * 100))

    nodes = set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values]))
    do_names = False
    if 'Names' in ListProperties:
        do_names = True
        names = lin_tree_information['Names']
        names_which_matter = {k: v for k, v in names.iteritems() if v != '' and k in nodes}
        names_formated = format_names(names_which_matter)
        order_on_nodes = np.array(names_formated.keys())[np.argsort(names_formated.values())]
        nodes.difference_update(set(order_on_nodes))

    f = open(name, "w")
    inv_lin_tree = {v: k for k, vals in lin_tree.iteritems() for v in vals}
    f.write("(tlp \"2.0\"\n")
    f.write("(nodes ")
    if do_names:
        for n in order_on_nodes:
            f.write(str(n) + " ")
    for n in nodes:
        f.write(str(n) + " ")
    f.write(")\n")

    nodes = set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values]))
    count_edges = 0
    for m, ds in lin_tree.iteritems():
        count_edges += 1
        for d in ds:
            f.write("(edge " + str(count_edges) + " " + str(m) + " " + str(d) + ")\n")
    f.write("(property 0 int \"id\"\n")
    f.write("\t(default \"0\" \"0\")\n")
    for node in nodes:
        f.write("\t(node " + str(node) + str(" \"") + str(node) + "\")\n")
    f.write(")\n")

    if do_names:  # Write Names
        f.write("(property 0 string \"Names\"\n")
        default = ''
        f.write("\t(default \"" + str(default) + "\" \"0\")\n")
        for node in nodes:
            f.write("\t(node " + str(node) + str(" \"") + str(names_formated.get(node, default)) + "\")\n")
        f.write(")\n")

    for property in properties:
        prop_name = property[0]
        vals = property[1]
        default = property[2]
        f.write("(property 0 string \"" + prop_name + "\"\n")
        f.write("\t(default \"" + str(default) + "\" \"0\")\n")
        for node in nodes:
            f.write("\t(node " + str(node) + str(" \"") + str(vals.get(node, default)) + "\")\n")
        f.write(")\n")
    f.write(")")
    f.close()


def read_lineage_tree(filename, t=None):
    """
    Return a lineage tree from a file if exist until time t other return an empty lineage tree
    :param filename:
    :param t:
    :return:
    """
    lin_tree_information = {}
    if os.path.exists(filename):
        f = open(filename)
        lin_tree_information = pkl.load(f)
        f.close()
        if t is not None:
            out = {}
            for k, v in lin_tree_information.iteritems():
                tmp = {}
                if k != 'lin_tree':
                    time = t + 1
                else:
                    time = t
                for key, val in v.iteritems():
                    if key / 10 ** 4 < time:
                        tmp[key] = val
                out[k] = tmp
            lin_tree_information = out
    return lin_tree_information


def write_lineage_tree(filename, lin_tree_information):
    """
    Write in a file the lineage tree structure
    :param filename:
    :param lin_tree_information:
    :return:
    """
    f = open(filename, 'w')
    pkl.dump(lin_tree_information, f)
    f.close()
