import pandas as pd
import numpy as np
import altair as alt

# configure altair style
from styling_and_visualisation import colorlib, default_theme

# disable altair max rows
alt.data_transformers.disable_max_rows()


def set_alt_theme():
    alt.themes.register("default_theme", default_theme)
    alt.themes.enable("default_theme")


# make a quick plot of a fasta alignment using BIO AlignIO
# slow for large alignments
def plot_alignment(file, alignment_format="fasta", plot_range=(0, None), seqlimit=None, protein=True, gaptoken='-', label_order=None):
    from Bio import AlignIO
    from core_functions.helper_functions import column_entropy

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

    return chart_aln, alignmentDF


# calculate cumulative sum and distribution for pd.Series
# takes pd.Series as input and returns a parsed DF and altair chart object
def plot_cumsum_counts(series, formatted_data=False, plot_type='merged', title='Chart', x_label='value', y_label='count',
                       x_min=1, y_min=1, x_max=None, y_max=None,
                       x_scale_type='log', y_scale_type='log', decimals=2):

    weighted_x_label = 'weighted_' + x_label
    weighted_y_label = 'weighted_' + y_label

    # skip formatting if supplied with preformatted data
    # example if several data output from this function are stacked using pd.concat
    if formatted_data:
        countDF = series

    else:
        # format DF for data handling, filter 0 values for plot
        # round to reduce float data display jaggedness
        series = series[series != 0].round(decimals)

        # format distribution dataframe
        countDF = pd.DataFrame(series.value_counts())
        countDF.columns = ['amount']
        countDF.sort_index(inplace=True)
        countDF.reset_index(inplace=True, names='size')
        countDF['cumsum'] = countDF['amount'].cumsum()
        countDF['frac_cumsum'] = countDF['cumsum'] / countDF['cumsum'].max()

        #calculate weighted amounts for examples the amount of clusters*sequences_in_those_clusters
        #gives ditribution of elements within the sets provided rather than distribution of the set sizes
        countDF['weighted_amount'] = countDF['amount'] * countDF['size']
        countDF['weighted_cumsum'] = countDF['weighted_amount'].cumsum()
        countDF['weighted_frac_cumsum'] = countDF['weighted_cumsum'] / countDF['weighted_cumsum'].max()

        #add title as column for stacking separate runs later on
        countDF['series_label'] = [title for _ in countDF.index]

    # rename columns for plotting
    countDF.columns = [x_label, y_label, 'cumsum', 'frac_cumsum', weighted_y_label, 'weighted_cumsum',
                       'weighted_frac_cumsum', 'series_label']

    # format axis domains
    x_range = [x_min, countDF[x_label].max()]
    y_range = [y_min, countDF[y_label].max()]

    if x_max:
        x_range = [x_min, x_max]
    if y_max:
        y_range = [y_min, y_max]

    # correct axis label placement
    area_axis_labelAlign = 'left'

    if plot_type in ['weighted', 'default']:
        line_axis_labelAlign = 'right'

    else:
        line_axis_labelAlign = 'left'

    # add artifical datapoint to display data in order to have consistent fill color
    # really bad fix but altair has issues with area log plots with multiple series
    countDF_plot = countDF.copy()
    for s in countDF.series_label.unique():
        countDF_plot.loc[countDF_plot.shape[0]+1] = [0,0.0001,0,0,0.00001,0,0,s]
        countDF_plot.loc[countDF_plot.shape[0]+1] = [countDF[x_label].max(),0.0001,0,0,0.00001,0,0,s]

    # plot cumulative distribution
    chart_cumsum = alt.Chart(countDF, title=title).mark_line(color=colorlib['twilight_shifted_r_perm'][2],
                                                             strokeWidth=3).encode(
        x=alt.X(x_label, title=x_label, scale=alt.Scale(domain=x_range, type=x_scale_type)),
        y=alt.Y('frac_cumsum', title=f'Cumulative sum of {y_label}', scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(labelAlign=line_axis_labelAlign)),
        color=alt.Color('series_label'),
        tooltip=alt.Tooltip(['series_label', x_label, y_label, 'frac_cumsum'])
    )

    # plot cumulative distribution
    chart_weighted_cumsum = alt.Chart(countDF, title=title).mark_line(color=colorlib['twilight_shifted_r_perm'][2],
                                                                      strokeWidth=3).encode(
        x=alt.X(x_label, title=x_label, scale=alt.Scale(domain=x_range, type=x_scale_type)),
        y=alt.Y('weighted_frac_cumsum', title=f'Weighted cumulative sum of {y_label}', scale=alt.Scale(domain=[0, 1]),
                axis=alt.Axis(labelAlign=line_axis_labelAlign)),
        color=alt.Color('series_label'),
        tooltip=alt.Tooltip(['series_label',x_label, weighted_y_label, 'weighted_frac_cumsum'])
    )

    # plot value distribution
    chart_bar = alt.Chart(countDF_plot).mark_line(interpolate='step-after',
                                             fillOpacity=0.05, line=True).encode(
        x=alt.X(x_label + ':Q', scale=alt.Scale(domain=x_range, type=x_scale_type),
                axis=alt.Axis(labelAlign=area_axis_labelAlign)),
        y=alt.Y(y_label, scale=alt.Scale(domain=y_range, type=y_scale_type)),
        color = alt.Color('series_label'),
        fill=alt.Fill('series_label'),
        tooltip=alt.Tooltip(['series_label', x_label, y_label, 'frac_cumsum'])
    )

    # plot value distribution
    chart_weighted_bar = alt.Chart(countDF_plot).mark_line(interpolate='step-after',
                                                      fillOpacity=0.05, line=True).encode(
        x=alt.X(x_label + ':Q', scale=alt.Scale(domain=x_range, type=x_scale_type),
                axis=alt.Axis(labelAlign=area_axis_labelAlign)),
        y=alt.Y(weighted_y_label, scale=alt.Scale(domain=y_range, type=y_scale_type)),
        color=alt.Color('series_label'),
        fill=alt.Fill('series_label'),
        tooltip=alt.Tooltip(['series_label', x_label, weighted_y_label, 'weighted_frac_cumsum'])
    )

    # merge and configure
    if plot_type=='weighted':
        chart_merge = alt.layer(chart_weighted_cumsum).resolve_scale(y='independent').interactive()
    if plot_type=='weighted_bar':
        chart_merge = alt.layer(chart_weighted_bar, chart_weighted_cumsum).resolve_scale(y='independent').interactive()
    elif plot_type=='default':
        chart_merge = alt.layer(chart_cumsum).resolve_scale(y='independent').interactive()
    elif plot_type=='default_bar':
        chart_merge = alt.layer(chart_bar, chart_cumsum).resolve_scale(y='independent').interactive()
    elif plot_type=='merged':
        chart_merge = alt.layer(chart_bar, chart_weighted_bar, chart_cumsum, chart_weighted_cumsum).resolve_scale(y='independent').interactive()

    return chart_merge, countDF
