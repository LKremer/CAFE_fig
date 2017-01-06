#!/usr/bin/python3.4


import ete3
import re
import os
import argparse
import csv


TEX_HEAD = '''% start of expanded gene family table
\\newcommand{\\domcolwidthA}{3.9cm}
\\newcommand{\\domcolwidthB}{6.7cm}
\\newcommand{\\domcolwidthC}{9cm}
\\newgeometry{left=1.5cm,right=1.5cm,top=0cm,bottom=0cm}
\\begin{landscape}
\\begin{centering}
\\small
\\begin{longtable}[c]{ c p{\\domcolwidthA} c c p{\\domcolwidthB} p{\\domcolwidthC} }
\\caption{Significant gene family contractions and expansions in Isoptera (i.e. the branch leading from the LCA of Bger and the three termites to the LCA of the three termites). "domain\\_content" denotes which domains occur in proteins of the cluster; the two numbers indicate how many of the genes in that cluster have the domain. "min|mean|max" refers to the family sizes of other (non-termite) insects. GO-terms are from Pfam2GO. Note that gene families that do not occur in the LCA of all 19 insect species violate CAFEs model and are thus not included here (e.g. Exo\\_endo phos\\_expansion). }\\\\
% start of header
& \\textbf{contractions:} & & & & \\\\
\\toprule
\\textbf{cluster} & \\textbf{domain content} & \\textbf{size change} & \\textbf{min|mean|max} & \\textbf{domain description} & \\textbf{domain GO-terms (Pfam2GO)} \\\\
& & & \\textbf{in other insects} & & \\\\
\\midrule
% end of header
'''
TEX_TAIL = '''% end of table content
\\bottomrule
\\end{longtable}
\\end{centering}
\\end{landscape}
'''


def get_asterisks(p):
    if p <= 0.001:
        return '***'
    if p <= 0.01:
        return '**'
    elif p <= 0.05:
        return '*'
    else:
        return ''


def get_node_ids(line):
    species2number = {}
    nwk = line.strip().split(':')[1] + ';'
    t = ete3.Tree(nwk, format=1)
    for node in t.traverse():
        if node.is_leaf():
            i = node.name.find('<')
            species = node.name[:i]
            number = int(node.name[i:].lstrip('<').rstrip('>'))
            species2number[species] = number
            node.name = species
        else:
            number = int(node.name.lstrip('<').rstrip('>'))
            node.name = number
    return t, species2number


def store_branch_pvalues(line):
    # parsing the madness that CAFE has brought upon us
    crazy_str = line.split(' (node ID, node ID): ')[1].strip()
    #braceless = [t.strip('()') for t in crazy_str.split()]
    #tuples = [tuple([int(n) for n in t.split(',')]) for t in braceless]
    node_ids = re.split('[ ,()]', crazy_str)
    node_ids_nonempty = [int(node_id) for node_id in node_ids if node_id]  # filtering empty values
    return node_ids_nonempty


def highlight_significantly_evolving_nodes(t, branch_pvalues_raw, node_id_list):
    branch_pvalue_strs = re.split('[ ,()]', branch_pvalues_raw)
    if set(branch_pvalue_strs) == {'', '-'}:  # p-value estimation failed
        branch_pvalues = [1.0 for p in branch_pvalue_strs if p]
    else:
        branch_pvalues = [float(p) for p in branch_pvalue_strs if p]
    for node in t.traverse():
        if node.is_root():
            continue
        node_list_pos = node_id_list.index(node.number)
        branch_pval = branch_pvalues[node_list_pos]
        node.p_value = branch_pval
    return


def get_clade_ids(number_tree, species2number, clades_of_interest):
    print('\nDetermining user-specified clades of interest:')
    clade2id = {}
    for c in clades_of_interest:
        name, species_str = c.split('=')
        species = species_str.split(',')
        if len(species) == 1:
            clade_id = species[0]
        elif len(species) > 1:
            clade_id = number_tree.get_common_ancestor(species).name
        else:
            raise Exception('Malformatted -c param')
        clade2id[name] = clade_id
        print('  {:.<15}node "{}"'.format(name, clade_id))
    return clade2id


def parse_functional_annotation(fpath):
    func_dict = {}
    with open(fpath, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for rowdict in reader:
            cluster = rowdict['cluster']
            if cluster not in func_dict:
                func_dict[cluster] = []
            func_dict[cluster].append(rowdict)
    return func_dict


def get_latex_col_data(cluster, n_before, n_after, leaf_counts, p_value, func_dict):
    domains = []
    domain_names = []
    dom_go_terms = []
    if cluster in func_dict:
        for rowdict in func_dict[cluster]:
            if int(rowdict['n_genes_with_domain']) <= 1:
                continue
            if rowdict['domain_go_terms'] == '':
                dom_go_terms.append('-')
            else:
                dom_go_terms.append(rowdict['domain_go_terms'])
            #dom_str = '{n_genes_with_domain:> 3}/{n_genes_in_cluster:> 3} {domain}'.format(**rowdict)
            dom_str = '\\href{{http://www.ebi.ac.uk/interpro/search?q={}}}{{{}}} ({}/{})'.format(
                rowdict['domain_id'].split('.')[0],
                rowdict['domain'],
                rowdict['n_genes_with_domain'],
                rowdict['n_genes_in_cluster'],
            )
            domains.append(dom_str)
            domain_names.append(rowdict['domain_name'])
        if not domains:
            domains.append('-')
            domain_names.append('-')
            dom_go_terms.append('-')
            
    tex_row = {
        'cluster': cluster.replace('_', '\_'),
        'dom_field': ' \\\\ '.join(domains).replace('_', '\_'),
        'domain_names': ' \\\\ '.join(domain_names).replace('_', '\_'),
        'n_before': n_before,
        'n_after': n_after,
        'min': min(leaf_counts),
        'mean': sum(leaf_counts) / len(leaf_counts),
        'max': max(leaf_counts),
        'p_value': p_value,
        'dom_go_terms': ' \\\\ '.join(dom_go_terms).replace('_', '\_'),
    }
    return tex_row


def main(cafe_out, family_to_plot, clades_of_interest, dont_show_tree, functional_annotation, write_report):
    if functional_annotation:
        cluster_func_dict = parse_functional_annotation(functional_annotation)
        tex_columns = []
    if write_report:
        assert functional_annotation, 'functional annotation is required for a report'
        report_dir = write_report
        os.mkdir(report_dir)
    with open(cafe_out, 'r') as cfout:
        in_data = False
        for line in cfout:
            # store the numerical node IDs:
            if line.startswith('# IDs of nodes:'):
                number_tree, species2number = get_node_ids(line)
                if clades_of_interest:
                    interesting_clades = get_clade_ids(
                        number_tree, species2number, clades_of_interest)
                    fam_change_counter = {c_id: {'name': c_name, 'expn_cnt': 0, 'cont_cnt': 0} for c_name, c_id in interesting_clades.items()}
            # get the order in which branch-specific p-values etc. are stored:
            if '(node ID, node ID):' in line:
                node_id_list = store_branch_pvalues(line)
            if line.strip().startswith('\'ID\'\t\'Newick\''):
                in_data = True
                continue
            if not in_data:
                continue
            # parsing family data:
            family, newick, family_pval, branch_pvalues = line.strip().split()
            if family_to_plot:
                if family != family_to_plot:
                    continue
            if float(family_pval) > 0.05 and not family_to_plot:
                continue
            newick += ';'  # add trailing semicolon (CAFE doesn't use that for some reason)
            newick = newick.replace(')_', ')')

            t = ete3.Tree(newick)
            fam_sizes = []
            for leaf in t:
                species, size_str = leaf.name.split('_')
                leaf.famsize = int(size_str)
                leaf.name = species
                leaf.number = species2number[species]
                fam_sizes.append(leaf.famsize)

            for node in t.traverse():
                if node.is_leaf():
                    continue
                # find out the numerical ID of the current node by getting the 
                # corresponding node in the numerically labelled tree
                children = [leaf.name for leaf in node]
                n_node = number_tree.get_common_ancestor(children)
                node.number = int(n_node.name)
                node.name = str(n_node.name)
                node.famsize = node.support

            highlight_significantly_evolving_nodes(
                t, branch_pvalues, node_id_list
            )

            # if the user specified clades of interest, only show families that changed in these clades!
            if clades_of_interest:
                irrelevant_clade = True
                for clade_name, clade_id in interesting_clades.items():
                    clade = t&clade_id
                    if clade.p_value <= 0.05:
                        irrelevant_clade = False
                if irrelevant_clade:
                    continue  # skip clade
                print('\n{}:'.format(family))

            for node in t.traverse():
                relative_size = 8 * (node.famsize / max(fam_sizes))

                if node.is_root():
                    node.color, asterisks = 'Gray', ''
                else:
                    asterisks = get_asterisks(node.p_value)
                    if node.p_value <= 0.05:
                        if node.famsize > node.up.famsize:
                            node.event = 'expanded'
                            node.color = '#4dac26'
                        elif node.famsize < node.up.famsize:
                            node.event = 'contracted'
                            node.color = '#d01c8b'
                        else:
                            assert False, 'significant p-value, but no fam size change?'

                        if clades_of_interest:
                            # count the number of expansions/contractions in clades of interest:
                            if node.name in fam_change_counter or node.number in fam_change_counter:
                                if node.name in fam_change_counter:
                                    node_identifier = node.name
                                else:
                                    node_identifier = node.number
                                if node.event == 'contracted':
                                    print('There is a {} (node {}) contraction'.format(
                                        fam_change_counter[node_identifier]['name'], node.number))
                                    fam_change_counter[node_identifier]['cont_cnt'] += 1
                                if node.event == 'expanded':
                                    print('There is a {} (node {}) expansion'.format(
                                        fam_change_counter[node_identifier]['name'], node.number))
                                    fam_change_counter[node_identifier]['expn_cnt'] += 1
                                if node.event in ('contracted', 'expanded') and functional_annotation:
                                    if functional_annotation:
                                        for dom_d in cluster_func_dict.get(family, []):
                                            print('  {domain} ({n_genes_with_domain}/{n_genes_in_cluster}): {domain_name}'.format(**dom_d))
                                        texcol = get_latex_col_data(
                                            cluster = family,
                                            n_before = node.up.famsize,
                                            n_after = node.famsize,
                                            leaf_counts = [leaf.famsize for leaf in t if leaf not in node],
                                            p_value = node.p_value,
                                            func_dict = cluster_func_dict,
                                        )
                                        tex_tuple = (node.event, node.p_value, texcol)
                                        tex_columns.append(tex_tuple)
                    else:
                        node.color = 'Gray'
                        node.event = 'unchanged'
                
                if not dont_show_tree or functional_annotation:
                    cf = ete3.CircleFace(
                        radius=relative_size,
                        color=node.color,
                        style='circle'
                    )
                    cf.opacity = 0.7
                    node.add_face(cf, column=10, position='float')
                    nstyle = ete3.NodeStyle()
                    nstyle['size'] = 0  # node.famsize
                    node.set_style(nstyle)
                    famsize_str = '{}{}\n'.format(asterisks, int(node.famsize))
                    tf = ete3.TextFace(famsize_str, fsize=4, fgcolor=node.color)
                    node.add_face(tf, column=1, position='float')
            if not dont_show_tree or functional_annotation:
                ts = ete3.TreeStyle()
                ts.title.add_face(ete3.TextFace(family, fsize=14), column=0)
                ts.scale = 0.5  # pixels per million years
                if not dont_show_tree:
                    t.show(tree_style=ts)
            if functional_annotation and write_report:
                pdf_path = os.path.join(report_dir, '{}.pdf'.format(family))
                t.render(pdf_path, tree_style=ts)

        if clades_of_interest:
            for cnt_dict in fam_change_counter.values():
                print()
                print(cnt_dict['name'] + ':')
                print(cnt_dict['cont_cnt'], 'contractions')
                print(cnt_dict['expn_cnt'], 'expansions')
        if functional_annotation and write_report:
            line = ('\\href{{run:data/{cluster}_IDs_heatmap.pdf}}{{ {cluster} }} & \\parbox[t]{{\\domcolwidthA}}{{ {dom_field} }} & \\href{{run:data/{cluster}.pdf}}{{ {n_before:.0f} $\\rightarrow$ {n_after:.0f} (p={p_value:.4f}) }} & {min:.0f} | {mean:.1f} | {max:.0f} & \\parbox[t]{{\\domcolwidthB}}{{ {domain_names} }} & \\parbox[t]{{\\domcolwidthC}}{{ {dom_go_terms} }} \\\\ \\midrule\n')
            sorted_cols = sorted(tex_columns, key=lambda t: (t[0], t[1]) )
            report_path = os.path.join(report_dir, 'expansions_contractions.tex')
            with open(report_path, 'w') as outfile:
                outfile.write(TEX_HEAD)
                for event, pvalue, rowdict in sorted_cols:
                    outfile.write(line.format(**rowdict))
                outfile.write(TEX_TAIL)
            print('Wrote LaTeX summary table to', report_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Parses a CAFE output file (.cafe) and plots signi'
                    'ficantly expanded/contracted families')
    parser.add_argument('cafe_out', help='the file report.cafe (or similar name)')
    parser.add_argument('-f', '--family_to_plot', help='which family should '
                        'be plotted (optional)')
    parser.add_argument('-c', '--clades_of_interest', nargs='+',
                        help='specifiy clades that should be screened for '
                        'expansions, e.g.: Isoptera=zne,mna Mnat=mna '
                        'Blattodea=bge,cse')
    parser.add_argument('-d', '--dont_show_tree', action='store_true',
                        help='don\'t open windows with trees')
    parser.add_argument('-a', '--functional_annotation', help='a TSV file that '
                        'contains a functional annotation of orthology clusters. '
                        'generated with link_orthoclusters_to_domains.py (optional)')
    parser.add_argument('-r', '--write_report', help='write a LaTeX report to the '
                        'specified directory', default=False)
    args = parser.parse_args()
    main(**vars(args))
