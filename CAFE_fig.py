#!/usr/bin/python3.4


import ete3
import argparse
import re
import functools




def is_valid_format(line):
    # check if a text line is in the valid format to denote a CAFE family
    values = line.strip().split('\t')
    if len(values) != 4:
        return False
    # the family-wide p-value should be here
    try:
        float(values[2])
    except ValueError:
        return False
    # the trees should start/end with parentheses
    fam_tree, pvalue_tree = values[1], values[3]
    if not fam_tree.startswith('('):
        return False
    if not pvalue_tree.startswith('(') or not pvalue_tree.endswith(')'):
        return False
    return True


class CAFE_fig():
    def __init__(self, report_cafe, families, clades, dump):
        self.graphics_options = {
            '+': '#4dac26',  # expansion
            '=': '#d3d3d3',  # unchanged (remain)
            '-': '#d01c8b',  # contraction (decrease)
            'pixels_per_mya': 0.5,  # pixels per million years (tree width)
            'opacity': 0.7,  # opacity of node circles
            'scale': 1.0,  # size scale factor of node circles
        }
        self.alpha = 0.05  # p-value cutoff
        self.report_path = report_cafe
        self.parse_tree()
        for node in self.tree.traverse():
            print('----------------')
            print('name', node.name)
            print('id', node.id)
            print('lambda', node.lambda_group)
            if not node.is_root():
                print('avg_expansion', node.avg_expansion)
                print('expansion', node.expansion)
                print('remain', node.remain)
                print('decrease', node.decrease)
        if families:
            self.families_of_interest = set(families)
        if clades:
            self.get_clades_of_interest(clades)


    def parse_tree(self):
        '''
        read the first few lines of the CAFE output file to extract the
        phylogeny, the output format/node order and the average expansion
        of all nodes.
        '''
        with open(self.report_path, 'r') as report:
            for line in report:
                if line.startswith('Tree:'):
                    # parse the general phylogeny including branch lengths
                    self.parse_phylo_tree(line)
                if line.startswith('Lambda tree:'):
                    # add the information of the lambda groups to the tree
                    self.parse_lambda_tree(line)
                if line.startswith('# IDs of nodes:'):
                    # add CAFE's numerical ID's to the tree
                    self.parse_node_num_tree(line)
                if line.startswith('# Output format for: '):
                    # find the output format that CAFE used to average expansion etc.
                    self.cafe_node_id_order = [int(i) for i in re.findall(r'\d+', line)]
                if line.startswith('Average Expansion'):
                    self.parse_fam_size_summary_tree(line, 'avg_expansion')
                if line.startswith('Expansion'):
                    self.parse_fam_size_summary_tree(line, 'expansion')
                if line.startswith('Remain'):
                    self.parse_fam_size_summary_tree(line, 'remain')
                if line.startswith('Decrease'):
                    self.parse_fam_size_summary_tree(line, 'decrease')
                if line.startswith('\'ID\''):
                    break  # end of header lines
        return

    def parse_fam_size_summary_tree(self, line, node_attr_name):
        '''
        parses a line that denotes node features (e.g. average expansion)
        and adds the feature as a class attribute to the corresponding node.
        '''
        node_fam_sizes = [float(size) for size in re.findall(r'[\d\.-]+', line)]
        for node_id, node_size in zip(self.cafe_node_id_order, node_fam_sizes):
            nodes = self.tree.search_nodes(id=node_id)
            assert len(nodes) == 1
            node = nodes[0]
            setattr(node, node_attr_name, node_size)
        return

    def parse_phylo_tree(self, line):
        '''
        parse the general phylogeny including branch lengths
        '''
        newick = line[5:].strip() + ';'
        self.tree = ete3.Tree(newick)
        return
        
    def parse_lambda_tree(self, line):
        '''
        find out in which lambda group each node is and add this information
        as a node attribute.
        '''
        lambda_nwk = line[12:].strip() + ';'
        lambda_tree = ete3.Tree(lambda_nwk)
        for node, lambda_node in zip(
            self.tree.traverse(),
            lambda_tree.traverse()
        ):
            if lambda_node.name:
                # ete3 parser calls this info "name" for leaves
                node.lambda_group = int(lambda_node.name)
            else:
                # ete3 parser calls this info "support" for internal nodes
                node.lambda_group = int(lambda_node.support)
            assert node.lambda_group
        return


    def parse_node_num_tree(self, line):
        '''
        parses the numerical ID that CAFE assigned to each node, and adds
        this ID as a node class feature.
        '''
        num_nwk = line[15:].strip() + ';'
        num_tree = ete3.Tree(num_nwk.replace('<', ' ').replace('>', ''))
        for node, num_node in zip(
            self.tree.traverse(),
            num_tree.traverse()
        ):
            if num_node.name:
                node.id = int(num_node.name.split()[-1])
            else:
                node.id = int(num_node.support)
        return


    def summary_tree(self):
        '''
        show a tree that visualizes the total number of expansions and 
        contractions across the whole phylogeny for allgene families.
        '''
        def fam_size_piechart_layout(node):
            '''
            the PieChart layout function, defined in local scope so it
            can access class attributes (graphics options).
            '''
            if not node.is_root():
                n_families = node.expansion + node.remain + node.decrease
                percentages = [
                    (node.expansion / n_families) * 100,
                    (node.remain / n_families) * 100,
                    (node.decrease / n_families) * 100,
                ]
                colors = [self.graphics_options['+'], self.graphics_options['='], self.graphics_options['-']]
                diameter = 13 * (1 + node.avg_expansion) * self.graphics_options['scale']
                piechart = ete3.PieChartFace(percentages, diameter, diameter, colors)
                piechart.opacity = self.graphics_options['opacity']
                ete3.faces.add_face_to_node(piechart, node, 0, position='float')
            nstyle = ete3.NodeStyle()
            nstyle['size'] = 0
            node.set_style(nstyle)
            return

        t = self.tree
        ts = ete3.TreeStyle()
        header = 'frequency of expansions and contractions across the phylogeny'
        ts.title.add_face(ete3.TextFace(header, fsize=8), column=0)
        ts.scale = self.graphics_options['pixels_per_mya']  # pixels per million years
        ts.layout_fn = fam_size_piechart_layout
        t.show(tree_style=ts)


    def get_clades_of_interest(self, clades_of_interest):
        '''
        parses the user-specified "--clades" parameter
        '''
        self.clades_of_interest = {}
        for c in clades_of_interest:
            name, species_str = c.split('=')
            species = self.tree.species_str.split(',')
            if len(species) == 1:
                node = self.tree.search_nodes(name=species[0])
            elif len(species) > 1:
                node = self.tree.get_common_ancestor(species)
            else:
                raise Exception('invalid --clades param')
            self.clades_of_interest.add(node)
        return


    def __iter__(self):
        with open(self.report_path, 'r') as report:
            for line in report:
                if is_valid_format(line):
                    yield Family(line, self)


    def show_fam_size_tree(self, family):
        def fam_size_layout(node):
            try:
                node_color = self.graphics_options[node.event]
            except:
                node_color = self.graphics_options['=']
            relative_size = 8 * (node.fam_size / max(family.fam_sizes))
            cf = ete3.CircleFace(
                radius=relative_size * self.graphics_options['scale'],
                color=node_color,
                style='circle'
            )
            cf.opacity = self.graphics_options['opacity']
            node.add_face(cf, column=10, position='float')
            # remove the silly default blue node dot
            nstyle = ete3.NodeStyle()
            nstyle['size'] = 0
            node.set_style(nstyle)
            return
        
        t = self.tree
        ts = ete3.TreeStyle()
        ts.layout_fn = fam_size_layout
        header = 'Evolution of the gene family "{}" (p={})'.format(
            family.name,
            family.pvalue,
        )
        ts.title.add_face(ete3.TextFace(header, fsize=8), column=0)
        ts.scale = self.graphics_options['pixels_per_mya']  # pixels per million years
        t.show(tree_style=ts)





class Family():
    def __init__(self, txtline, cafe_fig_instance):
        values = txtline.strip().split()
        self.name = values[0]
        self.nwk_famsize_str = values[1].replace(')_', ')') + ';'
        self.pvalue = float(values[2])
        self.branch_pvalue_str = values[3]
        self.c = cafe_fig_instance
        return


    def get_tree_with_famsizes(self):
        self.fam_sizes = []
        size_tree = ete3.Tree(self.nwk_famsize_str)
        phylo_tree = self.c.tree
        # parse family sizes:
        for node, size_tree_node in zip(
            phylo_tree.traverse(),
            size_tree.traverse()
        ):
            if size_tree_node.is_leaf():
                node.fam_size = int(size_tree_node.name.split('_')[1])
            else:
                node.fam_size = int(size_tree_node.support)
            self.fam_sizes.append(node.fam_size)
        # parse family pvalues:
        node_pvalues = re.findall(r'[\d\.]+|-', self.branch_pvalue_str)
        for node_id, node_size in zip(self.c.cafe_node_id_order, node_pvalues):
            nodes = self.c.tree.search_nodes(id=node_id)
            assert len(nodes) == 1
            node = nodes[0]
            if node_size == '-':
                node.pvalue = None
            else:
                node.pvalue = float(node_size)
                if node.pvalue <= self.c.alpha:
                    if node.fam_size > node.up.fam_size:
                        node.event = '+'
                    elif node.fam_size < node.up.fam_size:
                        node.event = '-'
                    else:
                        raise Exception('significant p-value, but no fam size change?')
        return



def main(report_cafe, families, clades, dump):
    # parse initial information (phylogeny and CAFE output formats)
    c = CAFE_fig(report_cafe, families, clades, dump)

    # show a tree that shows how many total expansions/contractions 
    # occured at each node
    c.summary_tree()

    for family in c:
        if family.pvalue > c.alpha:
            continue  # skip family since it's not significant
        if hasattr(c, 'families_of_interest') and family.name not in c.families_of_interest:
            continue  # skip family since the user didn't specifically select it
        fam_size_tree = family.get_tree_with_famsizes()
        if hasattr(c, 'clades_of_interest'):
            clade_pvalues = {n.pvalue for n in c.clades_of_interest}
            if clade_pvalues == {'-'}:  # no p-values were estimated for this family
                continue
            if min(clade_pvalues) > c.alpha:
                continue
        c.show_fam_size_tree(family)
        exit('bye')






if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parses a CAFE output file (.cafe) and plots signi'
                    'ficantly expanded/contracted families')
    parser.add_argument('report_cafe', help='the file report.cafe (or similar'
                        'name)')
    parser.add_argument('-f', '--families', help='only show these families',
                        nargs='+')
    parser.add_argument('-c', '--clades', help='only show families that are '
                        'expanded/contracted at this clade, e.g.: Isoptera=zne'
                        ',mna', nargs='+')
    parser.add_argument('-d', '--dump', action='store_true',
                        help='don\'t open trees in a window, write PDF files instead')
    #parser.add_argument('-a', '--functional_annotation', help='a TSV file that '
                        #'contains a functional annotation of orthology clusters. '
                        #'generated with link_orthoclusters_to_domains.py (optional)')
    #parser.add_argument('-r', '--write_report', help='write a LaTeX report to the '
                        #'specified directory', default=False)
    args = parser.parse_args()
    main(**vars(args))
