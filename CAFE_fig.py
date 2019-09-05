#!/usr/bin/python3


from __future__ import print_function
import ete3
import argparse
import re
import os
import shutil
import copy
import math
from base64 import b16encode


# Python 2 compatability:
try:
    input = raw_input
except NameError:
    pass


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


def get_pvalue_asterisks(node):
    try:
        p = float(node.pvalue)
    except (TypeError, AttributeError):
        # TypeError occurs when CAFE did not compute the pvalue because the whole
        # family did not experience significant size changes (p = None).
        # AttributeError occurs when the node has no pvalue attribute (e.g.
        # because it is the root node of the tree).
        return ''
    if p <= 0.0001:
        return '***'
    if p <= 0.001:
        return '**'
    if p <= 0.01:
        return '*'
    elif p <= 0.05:
        return '.'
    else:
        return ''

def to_rgb(v_abs, min_v, max_v):
    if min_v == max_v:
        hex_color = "#A9A9A9"
    else:
        v = (v_abs - min_v) / (max_v - min_v)  # scaled to a value between 0 and 1
        red = round(255 * v)
        blue = round(255 * (1 - v))
        rgb_triplet = (red, 0, blue)
        hex_color = '#' + b16encode(bytes(rgb_triplet)).decode()
    return hex_color


class CAFE_fig():
    def __init__(self, report_cafe, families, clades, pb, pf,
                 dump, gfx_output_format, count_all_expansions):
        self.graphics_options = {
            '+': '#4dac26',  # expansion
            '=': '#696969',  # unchanged (remain)
            '-': '#d01c8b',  # contraction (decrease)
            'pixels_per_mya': 1.0,  # pixels per million years (tree width)
            'opacity': 1.0,  # opacity of node circles
            'scale': 1.0,  # size scale factor of node circles
        }
        self.count_all_expansions = count_all_expansions
        self.branch_p_cutoff = pb
        self.family_p_cutoff = pf
        self.report_path = report_cafe
        self.prepare_pdf_dump(dump, gfx_format=gfx_output_format)
        self.parse_tree()
        if clades:
            self.get_clades_of_interest(clades)
        if families:
            self.families_of_interest = set(families)
        return

    def prepare_pdf_dump(self, dir_path, gfx_format='pdf'):
        '''
        create a directory to dump the figures (trees) to.
        '''
        self.gfx_format = '.' + gfx_format.lstrip('.')
        if self.gfx_format not in ('.svg', '.pdf', '.png'):
            raise Exception('graphics output format must be one of [svg|pdf|png]')
        if dir_path:
            if os.path.isdir(dir_path):
                answer = ''
                while answer.lower() not in ('y', 'n', 'yes', 'no'):
                    answer = input('The directory "{}" already exists. Overwrite it '
                                   'and delete all its contents? '
                                   '(y/n)? '.format(dir_path))
                if answer.lower() in ('n', 'no'):
                    exit('bye!')
                else:
                    shutil.rmtree(dir_path)
            os.mkdir(dir_path)
            self.out_dir = os.path.abspath(dir_path)
            fam_dir = os.path.join(dir_path, 'families')
            os.mkdir(fam_dir)
            self.out_dir_families = os.path.abspath(fam_dir)
            self.dump = True
        else:
            self.dump = False
        return

    def parse_tree(self):
        '''
        read the first few lines of the CAFE output file to extract the
        phylogeny, the output format/node order and the average expansion
        of all nodes.
        '''
        print('Parsing CAFE report...', end='\r')
        self.multi_lambda = False  # boolean that indicates whether the user ran CAFE
        # with one lambda or multiple lambda values, will be toggled if lambas found
        with open(self.report_path, 'r') as report:
            for line in report:
                if line.startswith('Tree:'):
                    # parse the general phylogeny including branch lengths
                    self.parse_phylo_tree(line)
                if line.startswith('Lambda:'):
                    self.parse_lambdas(line)
                if line.startswith('Lambda tree:'):
                    # add the information of the lambda groups to the tree
                    self.parse_lambda_tree(line)
                    self.multi_lambda = True  # user ran CAFE with more than one lambda
                if line.startswith('# IDs of nodes:'):
                    # add CAFE's numerical ID's to the tree
                    self.parse_node_num_tree(line)
                if line.startswith('# Output format for: '):
                    # find the output format that CAFE used to average expansion etc.
                    self.cafe_node_id_order = [int(i) for i in re.findall(r'\d+', line)]
                if line.startswith('Average Expansion'):
                    self.parse_fam_size_summary_tree(line, 'avg_expansion')
                if line.startswith('Expansion') or line.startswith('nExpansion'):
                    self.parse_fam_size_summary_tree(line, 'expansion')
                if line.startswith('Remain') or line.startswith('nRemain'):
                    self.parse_fam_size_summary_tree(line, 'remain')
                if line.startswith('Decrease') or line.startswith('nDecrease'):
                    self.parse_fam_size_summary_tree(line, 'decrease')
                if line.startswith('\'ID\''):
                    break  # end of header lines

        # count the number of significant expansions
        # and significant contractions per node:
        for family in self:
            if family.pvalue > self.family_p_cutoff:
                continue  # insignificant family
            family.get_tree_with_famsizes()
            for node, fam_tree_node in zip(
                self.tree.traverse(),
                family.tree.traverse()
            ):
                if not fam_tree_node.event:
                    continue
                if fam_tree_node.event == '+':
                    node.sig_expansions += 1
                elif fam_tree_node.event == '-':
                    node.sig_contractions += 1
        print('Parsing CAFE report... done!')
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
        for node in self.tree.traverse():
            node.sig_expansions = 0
            node.sig_contractions = 0
        return

    def parse_lambdas(self, line):
        self.lambdas = {}
        for i, lambda_str in enumerate(line.split()[1:]):
            self.lambdas[i + 1] = float(lambda_str)
        self.lambda_colors = {}
        max_l = math.log(max(self.lambdas.values()))
        min_l = math.log(min(self.lambdas.values()))
        for i, lambda_ in self.lambdas.items():
            self.lambda_colors[i] = to_rgb(math.log(lambda_), min_l, max_l)
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
                if self.count_all_expansions:
                    n_exp = node.expansion
                    n_con = node.decrease
                else:
                    n_exp = node.sig_expansions
                    n_con = node.sig_expansions

                if hasattr(self, 'lambda_colors'):
                    circle_color = self.lambda_colors[node.lambda_group]
                else:
                    circle_color = "blue"

                # add a text that shows expansions & contractions, e.g. +10/-20
                exp_cnt_txt = ete3.TextFace(
                    '+{} -{}\n'.format(int(n_exp), int(n_con)), fsize=6,
                    fgcolor=circle_color
                )
                pos = 'aligned' if node.is_leaf() else 'float'
                # add a circle that shows the average expansion
                ete3.faces.add_face_to_node(exp_cnt_txt, node, 1, position=pos)

                # add average expansion info:
                scale_factor = 1 + node.avg_expansion
                avg_exp = '{:+}'.format(round(node.avg_expansion, 2))
                circle = ete3.CircleFace(radius=9 * scale_factor,
                                         color=circle_color,
                                         label={'text': avg_exp, 'color': 'white',
                                                'fontsize': 3 + (2.25*scale_factor)})
                circle.opacity = self.graphics_options['opacity']
                ete3.faces.add_face_to_node(circle, node, 2, position='float')
            nstyle = ete3.NodeStyle()
            nstyle['size'] = 0
            node.set_style(nstyle)
            return
        t = self.tree
        ts = ete3.TreeStyle()
        header = 'Family expansions and contractions'
        if hasattr(self, "lambdas"):
            header += ('\nmin lambda: {} '
                       '(blue)\nmax lambda: {} (red)').format(
                        min(self.lambdas.values()), max(self.lambdas.values()))
        ts.title.add_face(ete3.TextFace(header, fsize=8), column=0)
        ts.scale = self.graphics_options['pixels_per_mya']  # pixels per million years
        ts.layout_fn = fam_size_piechart_layout
        self.show_or_dump_tree(tree_obj=t, tree_style=ts, fname='summary')
        return

    def get_clades_of_interest(self, clades_of_interest):
        '''
        parses the user-specified "--clades" parameter
        '''
        self.clades_of_interest = set()
        for c in clades_of_interest:
            name, species_str = c.split('=')
            species = species_str.split(',')
            if len(species) == 1:
                search_results = self.tree.search_nodes(name=species[0])
                assert len(search_results) == 1
                node = search_results[0]
            elif len(species) > 1:
                node = self.tree.get_common_ancestor(species)
            else:
                raise Exception('invalid --clades param')
            self.clades_of_interest.add(
                (name, node.id)
            )
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
            # add the family size number and asterisks to the figure
            famsize_str = '{}{}\n'.format(
                get_pvalue_asterisks(node),
                node.fam_size,
            )
            tf = ete3.TextFace(famsize_str, fsize=4, fgcolor=node_color)
            node.add_face(tf, column=1, position='float')
            # remove the silly default blue node dot
            nstyle = ete3.NodeStyle()
            nstyle['size'] = 0
            node.set_style(nstyle)
            return

        t = family.tree
        # write a quick summary of family name, pvalues and what happened
        if hasattr(self, 'clades_of_interest'):
            for clade_name, node_id in self.clades_of_interest:
                search_results = t.search_nodes(id=node_id)
                assert len(search_results) == 1
                clade_node = search_results[0]
                clade_event = getattr(clade_node, 'event', '=')
                clade_pvalue = getattr(clade_node, 'pvalue', '')
                tsv_header = 'family\tfamily_pvalue\tclade\tevent\tclade_pvalue'
                tsv_line = '{}\t{}\t{}\t{}\t{}'.format(
                    family.name, family.pvalue, clade_name, clade_event, clade_pvalue
                )
                if self.dump:
                    outf_name = '{}_summary.tsv'.format(clade_name)
                    outf_path = os.path.join(self.out_dir, outf_name)
                    if not os.path.isfile(outf_path):
                        with open(outf_path, 'a') as outf:
                            outf.write(tsv_header + '\n')
                    with open(outf_path, 'a') as outf:
                        outf.write(tsv_line + '\n')
                print('\n' + tsv_header)
                print(tsv_line)

        ts = ete3.TreeStyle()
        ts.layout_fn = fam_size_layout
        header = 'Evolution of the gene family "{}" (p={})'.format(
            family.name,
            family.pvalue,
        )
        ts.title.add_face(ete3.TextFace(header, fsize=8), column=0)
        ts.scale = self.graphics_options['pixels_per_mya']  # pixels per million years
        self.show_or_dump_tree(tree_obj=t, tree_style=ts,
                               fname=family.name, is_family=True)
        return

    def show_or_dump_tree(self, tree_obj, tree_style, fname, is_family=False):
        '''
        show the tree in a window, or write it to a PDF file if the user used --dump
        '''
        if self.dump:
            if is_family:
                out_dir = self.out_dir_families
            else:
                out_dir = self.out_dir
            out_path = os.path.join(out_dir, fname + self.gfx_format)
            print('\tWriting', os.path.relpath(out_path))
            tree_obj.render(out_path, tree_style=tree_style)
        else:
            tree_obj.show(tree_style=tree_style)
        return


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
        self.tree = copy.deepcopy(self.c.tree)
        # parse family sizes:
        for node, size_tree_node in zip(
            self.tree.traverse(),
            size_tree.traverse()
        ):
            if size_tree_node.is_leaf():
                node.fam_size = int(size_tree_node.name.split('_')[1])
            else:
                node.fam_size = int(size_tree_node.support)
            self.fam_sizes.append(node.fam_size)
            node.event = None
        # parse family pvalues:
        node_pvalues = re.findall(r'[\d\.]+|-', self.branch_pvalue_str)
        for node_id, node_size in zip(self.c.cafe_node_id_order, node_pvalues):
            nodes = self.tree.search_nodes(id=node_id)
            assert len(nodes) == 1
            node = nodes[0]
            if node_size == '-' or self.pvalue > self.c.family_p_cutoff:
                node.pvalue = None
            else:
                node.pvalue = float(node_size)
                if node.pvalue <= self.c.branch_p_cutoff:
                    if node.fam_size > node.up.fam_size:
                        node.event = '+'
                    elif node.fam_size < node.up.fam_size:
                        node.event = '-'
        return


def main(report_cafe, families, clades, pb, pf, dump, gfx_output_format, count_all_expansions):
    # parse initial information (phylogeny and CAFE output formats)
    c = CAFE_fig(report_cafe, families, clades, pb, pf, dump, gfx_output_format, count_all_expansions)

    # show a tree that shows how many total expansions/contractions
    # occured at each node
    c.summary_tree()

    # show a tree for each gene family, unless the user specified a filter
    # rule (specific families, or families that changed in a specific clade)
    for family in c:
        if hasattr(c, 'families_of_interest'):
            if family.name not in c.families_of_interest:
                continue  # skip family since the user didn't specifically select it
            family.get_tree_with_famsizes()  # prepare to plot
        else:
            if family.pvalue > c.family_p_cutoff:
                continue  # skip family since it's not significant
            family.get_tree_with_famsizes()
            if hasattr(c, 'clades_of_interest'):
                for __, node_id in c.clades_of_interest:
                    p_value = family.tree.search_nodes(id=node_id)[0].pvalue
                    if p_value <= c.branch_p_cutoff:
                        break
                else:  # loop wasnt broken = no significant event found
                    continue
        c.show_fam_size_tree(family)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parses a CAFE output file (.cafe) and plots a summary tree '
        'that shows the average expansion/contraction across the phylogeny; a tree '
        'that shows which clades evolved under the same lambda (if available); and '
        'a gene family evolution tree for each user-specified gene family.')
    parser.add_argument('report_cafe', help='the file report.cafe (or similar name)')
    parser.add_argument('-f', '--families', help='only show families with these IDs',
                        nargs='+')
    parser.add_argument('-c', '--clades', help='only show families that are '
                        'expanded/contracted at this clade. Format: [clade]='
                        '[leaf],[leaf] where clade is the name of the last '
                        'common ancestor of the two leaves, e.g.: Isoptera=zne,mna',
                        nargs='+')
    parser.add_argument('-pb', help='branch p-value cutoff (default: 0.05)',
                        default=0.05, type=float)
    parser.add_argument('-pf', help='family p-value cutoff (default: 0.05)',
                        default=0.05, type=float)
    parser.add_argument('-d', '--dump', help='don\'t open trees in a window, write '
                        'them to files in the specified directory instead (default: '
                        'off)', default=None)
    parser.add_argument('-g', '--gfx_output_format', default='.pdf', help='output '
                        'format for the tree figures when using --dump [svg|pdf|png]'
                        ' (default: pdf)')
    parser.add_argument('--count_all_expansions', action='store_true', help='count '
                        'and write down the number of *all* expansions and contrac'
                        'tions (default: only count significant expansions/contrac'
                        'tions)')
    args = parser.parse_args()
    if args.families:
        if args.clades:
            print('\n########\nWarning! "--families" overrides "--clades".\n########\n')
    main(**vars(args))
