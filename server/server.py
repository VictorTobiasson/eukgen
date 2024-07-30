from flask import Flask, jsonify, request, send_file, make_response
from flask_cors import CORS
from Bio import AlignIO
import pandas as pd
import altair as alt
alt.data_transformers.enable("vegafusion")
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
from helper_functions import column_entropy
from PyQt5.QtWidgets import QApplication
import queue
import threading


app = Flask(__name__)
CORS(app)

# Create a queue for render requests
render_queue = queue.Queue()

MICROCOSM_DIR = 'test_microcosm/'

def create_altair_chart(file, alignment_format="fasta", plot_range=(0, None), seqlimit=None, protein=True, gaptoken='-', label_order=None):

    # read fasta into AlignIO object
    alignment = AlignIO.read(file, alignment_format)
    alignment = alignment[0:seqlimit, plot_range[0]:plot_range[1]]

    # format AlignIO into pd.DataFrame for plotting
    alignmentDF = pd.DataFrame({align.description: list(align.seq) for align in alignment})
    alignmentDF['seqn'] = alignmentDF.index + plot_range[0]

    # calculate columnwise entropies for sequences
    cols = [alignment[:, col] for col in range(alignment.get_alignment_length())]
    if protein:
        alignmentDF['entropy'] = [column_entropy(string, gaptoken='-') for string in cols]
    else:
        alignmentDF['entropy'] = [column_entropy(string, protein=False, gaptoken='-') for string in cols]

    alignmentDF['opacity'] = [np.tanh(x / 2) for x in alignmentDF['entropy']]

    # melt for Altair plots
    alignmentDF_melt = alignmentDF.melt(id_vars=['seqn', 'entropy', 'opacity'],
                                        value_vars=alignmentDF.columns[0:-1],
                                        var_name='sequence', value_name='res')

    # define chart variables
    seq_start = alignmentDF_melt.seqn.min()
    seq_end = alignmentDF_melt.seqn.max()
    seq_len = seq_end - seq_start
    seq_number = alignmentDF_melt.sequence.unique().shape[0]

    # define protein color scheme and style
    if protein:
        resn = ['-', 'R', 'K', 'H', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
        colors = ['white', '#6276ba', '#7297c1', '#7297c1', '#b25652', '#b25652', '#b5b65e', '#94ae57', '#72a551',
                  '#72a551', '#cca389', '#c4ced4', '#95b5c7', '#bfa764', '#b5b65e', '#94ae57', '#72a551', '#cca389',
                  '#d8c7be', '#c4ced4', '#6276ba']
        entropy_domain = [0, 4]
        entropy_axis_values = [1, 2, 3, 4]

    # define nucleic color scheme and style
    else:
        resn = ['-', 'A', 'T', 'U', 'C', 'G', ]
        colors = ['white', '#6276ba', '#7297c1', '#7297c1', '#cca389', '#b25652']
        entropy_domain = [0, 2]
        entropy_axis_values = [0, 0.5, 1, 1.5, 2]

    if label_order != None:
        label_sort = label_order
    else:
        label_sort = None

    # base canvas
    chart_base = alt.Chart(alignmentDF_melt).encode(
        alt.X('seqn:O', axis=alt.Axis(values=list(range(seq_start, seq_end, 5)), grid=False)),
        alt.Y('sequence:O', sort=label_sort, axis=alt.Axis(grid=False, labelLimit=1000, title=None, labelFontSize=8)),
        alt.Opacity('opacity', legend=None)
    ).properties(width=seq_len * 8, height=seq_number * 8)

    # residue labels
    chart_text = chart_base.mark_text(color='black', align='center', fontSize=6.5).encode(
        alt.Text('res')
    )

    # colored boxes
    chart_box = chart_base.mark_rect().encode(
        alt.Color('res', scale=alt.Scale(domain=resn, range=colors)),
        alt.Tooltip(['sequence', 'seqn', 'res', 'entropy'])
    )

    # entropy bars
    bars = alt.Chart(alignmentDF).mark_bar().encode(
        alt.X('seqn:O', axis=alt.Axis(grid=False, labels=False, ticks=False, title=None)),
        alt.Y('entropy', axis=alt.Axis(values=entropy_axis_values, title='bits'),
              scale=alt.Scale(domain=entropy_domain)),
        alt.Color('entropy:Q', legend=None),
        alt.Tooltip('entropy')
    ).properties(width=seq_len * 8, height=40)

    # concat and plot layout
    chart_aln = alt.vconcat(bars, alt.layer(chart_box, chart_text),
                       spacing=8,
                       title=alt.TitleParams(text=file, fontSize=16)).resolve_legend('independent')

    # return json.loads(chart_aln.to_json())
    return chart_aln.to_dict(format="vega")

@app.route('/fetchMSASpec')
def fetchMSASpec():
    system = request.args.get('filename', '')
    if not system:
        return jsonify({"error": "No filename provided"}), 40
    msa = f'{MICROCOSM_DIR}{system}/{system}.merged.fasta.muscle'
    if not os.path.exists(msa):
        return jsonify({"error": "File not found"}), 404
    tree = Tree(f'{MICROCOSM_DIR}{system}/{system}.merged.fasta.muscle.treefile.annot')
    leaf_names = [leaf.name for leaf in tree.get_leaves()]
    try:
        chart_spec = create_altair_chart(msa, seqlimit=100, plot_range=(0,300), label_order=leaf_names)
        chart_spec["background"] = "null"
        return jsonify(chart_spec)
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    
def color_tree(treefile):
    # fixes some rendering crasher as per
    # https://github.com/etetoolkit/ete/issues/296
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    annot_tree = Tree(treefile)

    # set individual node styles
    default_node_style = NodeStyle()
    default_node_style['size'] = 0
    default_node_style['fgcolor'] = 'Black'

    LCA_euk_node_style = NodeStyle()
    LCA_euk_node_style['size'] = 10
    LCA_euk_node_style['fgcolor'] = '#5c7fe1'
    LCA_euk_node_style['bgcolor'] = '#eff3ff'

    LCA_prok_node_style = NodeStyle()
    LCA_prok_node_style['size'] = 8
    LCA_prok_node_style['fgcolor'] = 'Black'
    LCA_prok_node_style['bgcolor'] = 'LightGray'

    for node in annot_tree.traverse():
        node.set_style(default_node_style)

        if node.is_leaf():
            node.add_face(TextFace(node.name, fsize=8), column=1)
            if 'taxa' in node.features:
                node.add_face(TextFace(' ' + str(node.taxa), fsize=8), column=2)
                node.add_face(TextFace(' ' + str(node.counts), fsize=8), column=3)

        if 'LCA' in node.features:
            if node.LCA == 'Eukaryota':
                node.set_style(LCA_euk_node_style)
                #node.add_face(TextFace(node.LCA, fsize=8), column = 1)

                for leaf in node.get_leaves():
                    pass

            else:
                node.set_style(LCA_prok_node_style)
    annot_tree.ladderize()
    return annot_tree


def render_worker():
    while True:
        try:
            # Get render request from queue
            render_job = render_queue.get()
            if render_job is None:
                break  # Exit if None is received

            treefile, tree_img_path, system = render_job
            
            # Perform the actual rendering
            tree_colored = color_tree(treefile)
            ts = TreeStyle()
            ts.mode = 'r'
            ts.show_leaf_name = False
            ts.show_branch_length = False
            ts.show_branch_support = False
            ts.optimal_scale_level = ''
            ts.allow_face_overlap = True
            
            tree_colored.render(tree_img_path, tree_style=ts, units="px")
            
            print(f"Rendered tree for {system}")
        except Exception as e:
            print(f"Error rendering tree: {str(e)}")
        finally:
            render_queue.task_done()
    
@app.route('/fetchTreeImg')
def fetchTreeImg():
    system = request.args.get('filename', '')
    if not system:
        return jsonify({"error": "No filename provided"}), 40
    treefile = f'{MICROCOSM_DIR}{system}/{system}.merged.fasta.muscle.treefile.annot'
    if not os.path.exists(treefile):
        return jsonify({"error": "File not found"}), 404
    tree_img_path = f'{MICROCOSM_DIR}{system}/{system}.tmp_tree.png'
    # Add render job to queue
    render_queue.put((treefile, tree_img_path, system))
    # Wait for the job to complete
    render_queue.join()
    if os.path.exists(tree_img_path):
        return send_file(tree_img_path, mimetype='image/png')
    else:
        return jsonify({"error": "Failed to render tree"}), 500
    # try:

    #     ts = TreeStyle()
    #     ts.mode = 'r'
    #     ts.show_leaf_name = False
    #     ts.show_branch_length = False
    #     ts.show_branch_support = False
    #     ts.optimal_scale_level = ''
    #     ts.allow_face_overlap = True
    #     tree_annot = color_tree(treefile)
    #     tree_img_path = f'{MICROCOSM_DIR}{system}/{system}.tmp_tree.png'
    #     tree_annot.render(tree_img_path, tree_style=ts, units="px")

    #     return send_file(tree_img_path,  mimetype='image/png')
    #     # response = make_response(send_file(tree_img_path, mimetype='image/png'))
    #     # response.headers['Cache-Control'] = 'no-store, no-cache, must-revalidate, max-age=0'
    #     # response.headers['Pragma'] = 'no-cache'    
    #     # response.headers['Expires'] = '0'
    #     # return response
    # except Exception as e:
    #     return jsonify({"error": str(e)}), 500

@app.route('/list_microcosms')
def list_microcosms():
    files = [f for f in os.listdir(MICROCOSM_DIR) if f[0] != "."]
    return jsonify(files)

if __name__ == '__main__':
    # Create QApplication in the main thread
    qapp = QApplication(sys.argv)
    
    # Start the render worker in the main thread
    render_thread = threading.Thread(target=render_worker)
    render_thread.start()
    
    # Run Flask app
    app.run(debug=True, use_reloader=False)
    
    # Clean up
    render_queue.put(None)  # Signal worker to exit
    render_thread.join()
    qapp.quit()