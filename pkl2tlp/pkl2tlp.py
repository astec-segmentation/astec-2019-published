#!/usr/bin/env/python

import argparse
import cPickle as pkl
import numpy as np
from numbers import Number
from keys_for_properties import keydictionary

def format_names(names_which_matter):
    """
    Return an ensured formated cell names
    """
    tmp={}
    for k, v in names_which_matter.iteritems():
        try:
            val=v.split('.')[1][:-1]
            tmp[k] = (v.split('.')[0][0] + '%02d'%int(v.split('.')[0][1:]) + '.'
                      + '%04d'%int(v.split('.')[1][:-1]) + v.split('.')[1][-1])
        except Exception as e:
            print v
            exit()
    return tmp


def write_tlp_from_lin_tree(name, lin_tree_information, lin_name):
    """
    Write a lineage tree into an understable tulip file
    name : path to the tulip file to create
    lin_tree : lineage tree to write
    properties : dictionary of properties { 'Property name': [{c_id: prop_val}, default_val]}
    """
    
    lin_tree = lin_tree_information[lin_name]

    nodes = set(lin_tree.keys()).union(set([v for values in lin_tree.values() for v in values]))

    name_key = get_key_from_kdict('name', lin_tree_information)
    inv_lin_tree = {v:k for k, vals in lin_tree.iteritems() for v in vals}
    do_names = name_key is not None
    if do_names:
        names = lin_tree_information[name_key]
        names_which_matter = { k : v for k, v in names.iteritems() if v!='' and v!='NO' and k in nodes}
        names_formated = format_names(names_which_matter)
        order_on_nodes = np.array(names_formated.keys())[np.argsort(names_formated.values())]
        nodes.difference_update(set(order_on_nodes))
        tmp_names = {}
        for k, v in lin_tree_information[name_key].iteritems():
            if len(lin_tree.get(inv_lin_tree.get(k, -1), [])) != 1:
                tmp_names[k] = v
        lin_tree_information[name_key] = tmp_names

    f = open(name, "w")
    f.write("(tlp \"2.0\"\n")
    f.write("(nodes ")
    if do_names:
        for n in order_on_nodes:
           f.write(str(n)+ " ")
    for n in nodes:
        f.write(str(n)+ " ")
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

    for prop_name, property in lin_tree_information.iteritems():
        vals = property
        if (prop_name != lin_name and isinstance(vals, dict)
                and not isinstance(vals.values()[0], dict)):
            if type(vals.values()[0]) == str:
                default=''
                f.write("(property 0 string \""+prop_name+"\"\n")
                f.write("\t(default \""+str(default)+"\" \"0\")\n")
                for node in nodes:
                    f.write("\t(node " + str(node) + str(" \"") + str(vals.get(node, default)) + "\")\n")
                f.write(")\n") 
                if 'Fate' in prop_name or 'fate' in prop_name:
                    fates = set(vals.values())
                    corres = dict(zip(fates, range(len(fates))))
                    corres[''] = len(fates) + 1
                    f.write("(property 0 double \""+prop_name+"_int\"\n")
                    f.write("\t(default \""+str(corres[default])+"\" \"0\")\n")
                    for node in nodes:
                        f.write("\t(node " + str(node) + str(" \"") + str(corres[vals.get(node, default)]) + "\")\n")
                    f.write(")\n")
            elif isinstance(vals.values()[0], Number):
                default = np.median(vals.values())
                f.write("(property 0 double \""+prop_name+"\"\n")
                f.write("\t(default \""+str(default)+"\" \"0\")\n")
                for node in nodes:
                    f.write("\t(node " + str(node) + str(" \"") + str(vals.get(node, default)) + "\")\n")
                f.write(")\n") 
            elif (isinstance(vals.values()[0], np.ndarray) or isinstance(vals.values()[0], list)) and len(vals.values()[0]) == 3:
                f.write("(property 0 layout \"viewLayout\"\n")
                f.write("\t(default \"(0, 0, 0)\" \"()\")\n")
                for n in nodes:
                    f.write("\t(node " + str(n) + str(" \"") + str(tuple(vals.get(n, (0, 0, 0)))) + "\")\n")
                f.write(")\n")
                f.write("(property 0 layout \"Barycenter\"\n")
                f.write("\t(default \"(0, 0, 0)\" \"()\")\n")
                for n in nodes:
                    f.write("\t(node " + str(n) + str(" \"") + str(tuple(vals.get(n, (0, 0, 0)))) + "\")\n")
                f.write(")\n")

    f.write(")")
    f.close()

def get_key_from_kdict(key_name, DATA):
    key = None
    for k in keydictionary[key_name]['input_keys']:
        if k in DATA.keys():
            key = k
    return key


def main():
    parser = argparse.ArgumentParser(description='Convert pkl lineage into tulip lineage.')
    parser.add_argument('-i', '--input', help='input pickle .pkl file', required=True)
    parser.add_argument('-o', '--output', help='output tulip file (has to end with .tlp)', required=True)
    
    args = parser.parse_args()
    with open(args.input) as f:
        DATA = pkl.load(f)

    lin_name = get_key_from_kdict('lineage', DATA)
    if lin_name is None:
        print "There is no known key value for the lineage tree (expected 'lin_tree' or 'Lineage tree').\n\n\tThe script will not output a tulip file."
        exit()
    write_tlp_from_lin_tree(args.output, DATA, lin_name)


if __name__ == '__main__':
    main()
