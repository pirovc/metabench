#!/usr/bin/env python3
import argparse
import glob
import sys
import os
import json
import random
import math

import randomname  # pip install randomname https://github.com/beasteers/randomname

import pandas as pd

# Bokeh
from bokeh.io import save
from bokeh.core.enums import MarkerType
from bokeh.plotting import figure, output_file
from bokeh.models import ColumnDataSource, CDSView, Select, CustomJS, CustomJSTransform,CustomJSFilter, FactorRange, MultiSelect, CustomJSHover, Legend, LegendItem, HoverTool, Tabs, Panel, MultiChoice, Button, RadioButtonGroup, CheckboxGroup, CheckboxButtonGroup, Spacer, Range1d, TextInput, Whisker
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.filters import GroupFilter, IndexFilter
from bokeh.palettes import Category10, Category20, Colorblind, linear_palette, Turbo256

from bokeh.transform import linear_cmap, transform, factor_cmap
from bokeh.layouts import column, row


def main():

    parser = argparse.ArgumentParser(
        prog="plot metabench with bokeh", conflict_handler="resolve", add_help=True)
    parser.add_argument("-i", "--input",     metavar="",
                        type=str, nargs="*", required=True, help="json file(s) and/or folder(s) with json reports (recursive search)")
    parser.add_argument("-o", "--output",    metavar="",
                        type=str, help="output html file")
    parser.add_argument("-f", "--font-size",    metavar="",
                        type=int, default=10, help="font size")
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s 1.0.0")

    if len(sys.argv) == 1:
        parser.print_help()
        return False

    args = parser.parse_args()

    default_ranks = ["root",
                     "superkingdom",
                     "phylum",
                     "class",
                     "order",
                     "family",
                     "genus",
                     "species",
                     "assembly"]

    json_files = find_json_files(args.input)
    parsed_json = parse_json(json_files)

    report_types = ["build", "binning", "profiling"]
    report_elements = ["config", "bench", "evals"]
    tables = {}  # tables[type][element]
    for t in report_types:
        tables[t] = {}
        for e in report_elements:
            tables[t][e] = pd.DataFrame()


    for p, pjson in parsed_json.items():
        rep_type = pjson["report"]

        # Parse config element
        tables[rep_type]["config"] = pd.concat([tables[rep_type]["config"], load_config(pjson)], ignore_index=True)

        # Parse bench element
        # concat by rank and get same index
        tables[rep_type]["bench"] = pd.concat([tables[rep_type]["bench"], load_elements(pjson["bench"])], axis=1, ignore_index=True)

        # Load evals if present
        if "evals" in pjson and pjson["evals"]:
            tables[rep_type]["evals"] = pd.concat([tables[rep_type]["evals"], load_elements(pjson["evals"])], axis=1, ignore_index=True)
        else:
            print(p)

    # for t in report_types:
    #     for e in report_elements:
    #         print(t,e)
    #         print(tables[t][e])

    # Make unique random names (same as number of config files)
    rnd_names = set()
    while len(rnd_names) < len(parsed_json):
        rnd_names.add(randomname.get_name())

    main_tabs = []

    tools = "pan,wheel_zoom,box_zoom,box_select,tap,reset,save"

    if not tables["profiling"]["config"].empty:
        main_tabs.append(Panel(child=plot_metabench(
            "profiling", tables, default_ranks, tools, rnd_names, args), title="Profiling"))

    if not tables["binning"]["config"].empty:
        main_tabs.append(Panel(child=plot_metabench(
            "binning", tables, default_ranks, tools, rnd_names,  args), title="Binning"))

    if not tables["build"]["config"].empty:
        main_tabs.append(Panel(child=plot_metabench(
            "build", tables, default_ranks, tools, rnd_names,  args), title="Build"))

    main_layout = Tabs(tabs=main_tabs)

    output_file(args.output, title="Metabench", mode="inline")
    save(main_layout)

    return True


def plot_metabench(report, tables, default_ranks, tools, rnd_names, args):

    #
    # CDS
    #

    df_config = parse_df_config(tables[report]["config"], rnd_names)
    print(df_config)
    cds_config = ColumnDataSource(df_config)

    df_bench = parse_df_data(tables[report]["bench"])
    print(df_bench)
    cds_bench = ColumnDataSource(df_bench)

    df_evals = parse_df_data(tables[report]["evals"])
    print(df_evals)
    cds_evals = ColumnDataSource(df_evals)

    # Check if all jsons are complete (with config, bench and evals)
    if max(df_config.index) != max(df_bench.config) or \
       (not df_evals.empty and max(df_config.index) != max(df_evals.config)):
        raise Exception("Mismatching files")

    #
    # Vars
    #

    metrics = []
    for metric in set(cds_evals.data["metric"]):
        v = metric.split("|")
        metrics.append((metric, v[1] + " (" + v[0] + ")"))
    metrics.sort(key=lambda a: a[1])
    init_metric = next(iter(metrics))[0] if metrics else ""
    init_rank = "species" if "species" in default_ranks else default_ranks[-1]

    metrics_bench = list(sorted(set(cds_bench.data["metric"])))
    init_metric_bench = metrics_bench[0]

    smarkers, scolor = define_markers(df_config.shape[0])

    #
    # Widgets
    #

    # Metric (y)
    select_metric = Select(
        title="Evaluation metric (y):", value=init_metric, options=metrics)

    # Rank
    select_rank = Select(
        title="Rank:", value=init_rank, options=default_ranks)

    # Ranges
    radio_ranges = RadioButtonGroup(labels=["min-max", "0-1", "0-100"], active=0)

    # Group
    multiselect_groups = MultiSelect(title="Group by",
        value=["name"], options=df_config.columns.to_list(), size=9)

    # Marker
    multiselect_markers = MultiSelect(title="Marker",
        value=["name"], options=df_config.columns.to_list(), size=9, width=150)

    # Color
    multiselect_colors = MultiSelect(title="Color",
        value=["name"], options=df_config.columns.to_list(), size=9, width=150)

    # Sort
    select_sort = Select(title="Sort (x):", value="name", options=df_config.columns.to_list())
    
    # Toggle options
    toggle_boxplot = CheckboxGroup(labels=["Show Boxplot"], active=[])
    toggle_legend = CheckboxGroup(labels=["Show legend"], active=[0])
    toggle_label = CheckboxGroup(labels=["Show labels (x)"], active=[0])

    # Bench Metric
    select_metric_bench = Select(
        title="Benchmark metric:", value=init_metric_bench, options=metrics_bench)
    
    # Metric (x)
    select_metric_x_compare = Select(
        title="Evaluation metric (x):", value=init_metric, options=metrics)

    # Rank checkboxes
    checkbox_ranks = CheckboxGroup(
        labels=default_ranks, active=list(range(len(default_ranks))))


    #
    # Filters, View
    #

    # Table
    filter_config_evals = IndexFilter(indices=[])
    filter_config_bench = IndexFilter(indices=[])

    # Ranks
    filter_metric_ranks = GroupFilter(column_name="metric", group=init_metric)
    view_ranks = CDSView(source=cds_evals, filters=[filter_metric_ranks, filter_config_evals])

    # Compare
    filter_rank_compare = GroupFilter(column_name="rank", group=init_rank)
    filter_metric_compare = GroupFilter(column_name="metric", group=init_metric)
    view_compare = CDSView(source=cds_evals, filters=[filter_rank_compare, filter_metric_compare, filter_config_evals])
        
    # Metric
    filter_metric_ranks = GroupFilter(column_name="metric", group=init_metric)
    view_ranks = CDSView(source=cds_evals, filters=[filter_metric_ranks, filter_config_evals])

    # Evals
    filter_rank_evals = GroupFilter(column_name="rank", group=init_rank)
    filter_metric_evals = GroupFilter(column_name="metric", group=init_metric)
    view_evals = CDSView(source=cds_evals, filters=[filter_rank_evals, filter_metric_evals, filter_config_evals])

    # Bench
    filter_metric_bench = GroupFilter(column_name="metric", group=init_metric_bench)
    view_bench = CDSView(source=cds_bench, filters=[filter_config_bench, filter_metric_bench])

    #
    # Plots
    #

    # Table
    table_config, widgets_config = plot_datatable(cds_config, df_config)

    # Ranks
    plot_ranks = figure(title="Ranks", x_range=default_ranks,
                        toolbar_location="above", tools=tools,
                        width=500, height=500)
    plot_ranks.scatter(x="rank", y="value",
                       source=cds_evals,
                       view=view_ranks,
                       size=12,
                       marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                       fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                       line_color="black")
    plot_ranks.add_tools(make_hover(cds_config, exclude_list=["index", "fixed_arguments"]))
    plot_ranks.yaxis.axis_label = init_metric
    plot_ranks.xaxis.major_label_orientation = "vertical"

    # Compare
    plot_compare = figure(title="Compare", toolbar_location="above", tools=tools, width=500, height=500)

    # Select metric for y-axis based on the stacked cds (assumes data is sorted)
    get_metric_x = CustomJSTransform(args=dict(cds_evals=cds_evals,
                                               select_metric=select_metric,
                                               select_metric_x_compare=select_metric_x_compare),
                                     v_func="""
        // find start position of the y-axis metric
        var st_pos = 0;
        for (let i = 0; i < xs.length; i++) {
            if(xs[i]==select_metric_x_compare.value){
                st_pos=i;
                break;
            }
        }

        // Return values starting from position st_pos (assumes data is sorted)
        const y = new Float64Array(xs.length);
        var x = 0;
        for (let i = 0; i < xs.length; i++) {
            if(xs[i]==select_metric.value){
                y[i] = cds_evals.data["value"][st_pos+x];
                x++;
            }else{
                y[i] = -1;
            }
        }
        return y;
    """)

    plot_compare.scatter(x=transform("metric", get_metric_x),
                         y="value",
                         source=cds_evals,
                         view=view_compare,
                         size=12,
                         marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                         fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                         line_color="black")
    plot_compare.add_tools(make_hover(cds_config, exclude_list=["index", "fixed_arguments"]))
    plot_compare.xaxis.axis_label = init_metric
    plot_compare.yaxis.axis_label = init_metric


    # Evals
    plot_evals, legend_plot_evals, cds_boxplot_evals = main_plot(tools, cds_config, cds_evals, view_evals, smarkers, scolor, multiselect_groups, multiselect_markers, multiselect_colors, init_metric, args.font_size)


    # Bench
    plot_bench, legend_plot_bench, cds_boxplot_bench = main_plot(tools, cds_config, cds_bench, view_bench, smarkers, scolor, multiselect_groups, multiselect_markers, multiselect_colors, init_metric, args.font_size)

    #
    # Callbacks
    # 

    # Table
    cb_cds_config = CustomJS(
        args=dict(cds_evals=cds_evals,
                  cds_bench=cds_bench,
                  filter_config_evals=filter_config_evals,
                  filter_config_bench=filter_config_bench),
        code="""
        var indices = [];
        for(let i = 0; i < cds_evals.length; i++){
            if(this.indices.indexOf(cds_evals.data["config"][i]) > -1){
                indices.push(i);
            }
        }
        filter_config_evals.indices = indices;
        cds_evals.change.emit();

        var indices = [];
        for(let i = 0; i < cds_bench.length; i++){
            if(this.indices.indexOf(cds_bench.data["config"][i]) > -1){
                indices.push(i);
            }
        }
        filter_config_bench.indices = indices;
        cds_bench.change.emit();
        """)
    cds_config.selected.js_on_change('indices', cb_cds_config)

    # Ranks
    cb_checkbox_ranks = CustomJS(
        args=dict(plot_ranks=plot_ranks,
                  default_ranks=default_ranks),
        code="""
        plot_ranks.x_range.factors = this.active.map(function(x) { return default_ranks[x]; });
        """)
    checkbox_ranks.js_on_click(cb_checkbox_ranks)

    # Compare

    # Only triggers cds to apply transform on y-axis
    cb_select_metric_x_compare = CustomJS(
        args=dict(cds_evals=cds_evals,
                  xaxis=plot_compare.xaxis[0]),
        code="""
        xaxis.axis_label = this.value;
        cds_evals.change.emit();
        """)
    select_metric_x_compare.js_on_change('value', cb_select_metric_x_compare)

    # Evals
    cb_multiselect_groups_evals, cb_multiselect_markers_color_evals, cb_toggle_boxplot_evals, cb_toggle_label_evals, cb_toggle_legend_evals = main_callbacks(cds_config, cds_evals, plot_evals, view_evals, legend_plot_evals, cds_boxplot_evals, multiselect_groups, multiselect_markers, multiselect_colors, select_sort, toggle_boxplot)

    multiselect_markers.js_on_change('value', cb_multiselect_markers_color_evals, cb_multiselect_groups_evals)
    multiselect_colors.js_on_change('value', cb_multiselect_markers_color_evals, cb_multiselect_groups_evals)
    multiselect_groups.js_on_change('value', cb_multiselect_groups_evals, cb_toggle_boxplot_evals)
    select_sort.js_on_change('value', cb_multiselect_groups_evals)
    cds_config.selected.js_on_change('indices', cb_multiselect_groups_evals, cb_toggle_boxplot_evals)
    toggle_boxplot.js_on_click(cb_toggle_boxplot_evals)
    toggle_label.js_on_click(cb_toggle_label_evals)
    toggle_legend.js_on_click(cb_toggle_legend_evals)


    cb_radio_ranges = CustomJS(
        args=dict(radio_ranges=radio_ranges,
                  cds_evals=cds_evals,
                  y_range_evals=plot_evals.y_range,
                  y_range_ranks=plot_ranks.y_range,
                  y_range_compare=plot_compare.y_range,
                  x_range_compare=plot_compare.x_range),
        code="""
        // workaround found using _initial_start/end: https://discourse.bokeh.org/t/datarange1d-update-problems/9974/4
        
        if (radio_ranges.active==0){
            // auto (reset default)
            y_range_evals._initial_start = null;
            y_range_evals._initial_end = null;

            y_range_ranks._initial_start = null;
            y_range_ranks._initial_end = null;

            y_range_compare._initial_start = null;
            y_range_compare._initial_end = null;
            x_range_compare._initial_start = null;
            x_range_compare._initial_end = null;

            // Call change to update null values 
            cds_evals.change.emit();

        }else if(radio_ranges.active==1){
            // 0-1
            const spacer = 0.05;
            y_range_evals._initial_start = 0-spacer;
            y_range_evals._initial_end = 1+spacer;
            y_range_evals.start = 0-spacer;
            y_range_evals.end = 1+spacer;

            y_range_ranks._initial_start = 0-spacer;
            y_range_ranks._initial_end = 1+spacer;
            y_range_ranks.start = 0;
            y_range_ranks.end = 1+spacer;

            y_range_compare._initial_start = 0-spacer;
            y_range_compare._initial_end = 1+spacer;
            y_range_compare.start = 0-spacer;
            y_range_compare.end = 1+spacer;
            x_range_compare._initial_start = 0-spacer;
            x_range_compare._initial_end = 1+spacer;
            x_range_compare.start = 0-spacer;
            x_range_compare.end = 1+spacer;

        }else if(radio_ranges.active==2){
            // 0-100
            const spacer = 5;
            y_range_evals._initial_start = 0-spacer;
            y_range_evals._initial_end = 100+spacer;
            y_range_evals.start = 0-spacer;
            y_range_evals.end = 100+spacer;

            y_range_ranks._initial_start = 0-spacer;
            y_range_ranks._initial_end = 100+spacer;
            y_range_ranks.start = 0;
            y_range_ranks.end = 100+spacer;

            y_range_compare._initial_start = 0-spacer;
            y_range_compare._initial_end = 100+spacer;
            y_range_compare.start = 0-spacer;
            y_range_compare.end = 100+spacer;
            x_range_compare._initial_start = 0-spacer;
            x_range_compare._initial_end = 100+spacer;
            x_range_compare.start = 0-spacer;
            x_range_compare.end = 100+spacer;
        }
        """)
    radio_ranges.js_on_click(cb_radio_ranges)

    cb_select_metric = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_metric_ranks=filter_metric_ranks,
                  filter_metric_compare=filter_metric_compare,
                  filter_metric_evals=filter_metric_evals,
                  yaxis_ranks=plot_ranks.yaxis[0],
                  yaxis_compare=plot_compare.yaxis[0],
                  yaxis_evals=plot_evals.yaxis[0],
                  select_rank=select_rank),

        code="""
        filter_metric_ranks.group = this.value;
        yaxis_ranks.axis_label = this.value;

        filter_metric_compare.group = this.value;
        yaxis_compare.axis_label = this.value;

        filter_metric_evals.group = this.value;
        yaxis_evals.axis_label = this.value + " (" + select_rank.value + ")";

        cds_evals.change.emit();

        //console.log(Bokeh.documents[0].get_model_by_id('my_select'))
        //radio_ranges.trigger_event(({"event_name": "click"}))
        """)
    select_metric.js_on_change('value', cb_select_metric, cb_toggle_boxplot_evals)

    cb_select_rank = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_rank_compare=filter_rank_compare,
                  filter_rank_evals=filter_rank_evals),
        code="""
        filter_rank_compare.group = this.value;
        filter_rank_evals.group = this.value;
        cds_evals.change.emit();
        """)
    select_rank.js_on_change('value', cb_select_rank, cb_toggle_boxplot_evals)


    # Bench
    cb_multiselect_groups_bench, cb_multiselect_markers_color_bench, cb_toggle_boxplot_bench, cb_toggle_label_bench, cb_toggle_legend_bench = main_callbacks(cds_config, cds_bench, plot_bench, view_bench, legend_plot_bench, cds_boxplot_bench, multiselect_groups, multiselect_markers, multiselect_colors, select_sort, toggle_boxplot)

    multiselect_markers.js_on_change('value', cb_multiselect_markers_color_bench, cb_multiselect_groups_bench)
    multiselect_colors.js_on_change('value', cb_multiselect_markers_color_bench, cb_multiselect_groups_bench)
    multiselect_groups.js_on_change('value', cb_multiselect_groups_bench, cb_toggle_boxplot_bench)
    select_sort.js_on_change('value', cb_multiselect_groups_bench)
    cds_config.selected.js_on_change('indices', cb_multiselect_groups_bench, cb_toggle_boxplot_bench)
    toggle_boxplot.js_on_click(cb_toggle_boxplot_bench)
    toggle_label.js_on_click(cb_toggle_label_bench)
    toggle_legend.js_on_click(cb_toggle_legend_bench)
    
    # Bench
    cb_select_metric_bench = CustomJS(
        args=dict(cds_bench=cds_bench,
                  filter_metric_bench=filter_metric_bench,
                  yaxis_bench=plot_bench.yaxis[0]),

        code="""
        filter_metric_bench.group = this.value;
        yaxis_bench.axis_label = this.value;
        cds_bench.change.emit();
        """)
    select_metric_bench.js_on_change('value', cb_select_metric_bench, cb_toggle_boxplot_bench)

    tabs = []
    if not df_evals.empty:
        evals_layout = column([row([select_metric,select_rank,radio_ranges]),
                               plot_evals,
                               row(column([plot_compare,
                                           select_metric_x_compare]),
                                   plot_ranks,
                                   checkbox_ranks)])
        tabs.append(Panel(child=evals_layout, title="Evaluations"))

    bench_layout = column(select_metric_bench,
                          plot_bench)

    tabs.append(Panel(child=bench_layout, title="Benchmark"))


    layout = column([row(table_config, column(*widgets_config)),
                     row([Tabs(tabs=tabs),
                          column([multiselect_groups,
                                  row([multiselect_markers, multiselect_colors]),
                                  select_sort,
                                  toggle_boxplot,
                                  toggle_legend,
                                  toggle_label]),
                         ])
                    ])


    return layout


def find_json_files(input_list):
    files = []
    for i in input_list:
        if i.endswith(".json"):
            files.append(i)
        elif os.path.exists(i):
            for file in glob.glob(i + '/**/*.json', recursive=True):
                files.append(file)
    return files


def parse_json(json_files):
    parsed = {}
    for file in json_files:
        with open(file, "r") as jfile:
            parsed[file] = json.load(jfile)
    return parsed


def load_config(pjson):
    return pd.DataFrame([pjson["config"].values()], columns=pjson["config"].keys())


def load_elements(pjson_ele):
    cols = []
    vals = []
    for m1, v1 in pjson_ele.items():
        if isinstance(v1, dict):
            for m2, v2 in v1.items():
                if isinstance(v2, dict):
                    for m3, v3 in v2.items():
                        cols.append(("|".join([m1, m2]), m3,))
                        vals.append(v3)
                else:
                    cols.append((m1, m2,))
                    vals.append(v2)
        else:
            cols.append((m1, "",))
            vals.append(v1)

    return pd.DataFrame(vals, index=pd.MultiIndex.from_tuples(cols, names=["metric", "rank"]))


def make_color_palette(n_colors, linear: bool = False, palette: dict = None):
    if isinstance(palette, dict) and n_colors <= max(palette.keys()):
        # Special case for 1 and 2 (not in palettes)
        palette = palette[3 if n_colors < 3 else n_colors]

    if linear or n_colors > 20:
        if not palette:
            palette = Turbo256
        if n_colors <= 256:
            return linear_palette(palette, n_colors)
        else:
            # Repeat colors
            return [palette[int(i * 256.0 / n_colors)] for i in range(n_colors)]
    else:
        # Select color palette based on number of requested colors
        # Return the closest palette with most distinc set of colors
        if not palette:
            if n_colors <= 8:
                palette = Colorblind[8]
            elif n_colors <= 10:
                palette = Category10[10]
            elif n_colors <= 20:
                palette = Category20[20]
            else:
                palette = Turbo256

        return palette[:n_colors]

def parse_df_config(df, rnd_names):
    # Main dataframes and cds
    df_config = pd.DataFrame(df).fillna('')
    # Add toll name to version (in case of same version number among tools)
    df_config["version"] = df_config["version"] + \
        " (" + df_config["tool"] + ")"
    # Create random names for each configuration
    df_config["name"] = [rnd_names.pop()
                         for i in range(df_config.shape[0])]

    return df_config


def parse_df_data(df):
    df_data = pd.DataFrame(columns=["metric", "rank", "config", "value"])
    if not df.empty:
        # Stack metrics and reset index
        df_data = pd.DataFrame(
            df.stack().reset_index().set_index("metric"))
        # Rename column headers
        df_data.rename(columns={'level_2': "config", 0: 'value'}, inplace=True)
        # Make sure is sorted to properly use some js methods
        df_data.sort_values(by=["metric", 'rank', "config"], inplace=True)
    return df_data

def plot_datatable(cds_config, df_config):
    widgets_config = []
    table_config = DataTable(source=cds_config,
                             columns=[TableColumn(field=field, title=field)
                                      for field in df_config.columns],
                             selectable="checkbox",
                             width=1000,
                             height=150)
    

    choices = []
    for h in df_config.columns:
        for u in df_config[h].unique():
            if u:
                choices.append((h+"|"+u, u + " [" + h + "]"))

    multichoice_config = MultiChoice(options=choices, placeholder="Fixed filter")
    textinput_config = TextInput(placeholder="Partial filter e.g. database:abfv|sample:cami")
    radiobutton_config = RadioButtonGroup(labels=["AND", "OR"], active=0)
    button_config = Button(label="Apply", button_type="success")
    widgets_config.extend([multichoice_config, textinput_config, radiobutton_config, button_config])

    cb_button_config = CustomJS(
        args=dict(cds_config=cds_config,
                  multichoice_config=multichoice_config,
                  textinput_config=textinput_config,
                  radiobutton_config=radiobutton_config),
        code="""
        var selected_indices = [];
        
        if (multichoice_config.value.length > 0 || textinput_config.value.length > 0){
            for (var i = 0; i < cds_config.length; i++) {

                var found = 0;
                for (var m=0; m < multichoice_config.value.length; ++m){
                    const val_mc = multichoice_config.value[m].split("|");
                    if(cds_config.data[val_mc[0]][i]==val_mc[1]){
                        found++;
                    }
                }

                var text_filters_len = 0;
                if(textinput_config.value.length > 0){
                    const text_filters = textinput_config.value.split("|");
                    text_filters_len = text_filters.length;
                    for(var t = 0; t < text_filters.length; t++) {
                        const val_text = text_filters[t].split(":");
                        if(val_text.length==2 && val_text[0] in cds_config.data){
                            if(cds_config.data[val_text[0]][i].indexOf(val_text[1]) > -1){
                                found++;
                            }
                        }

                    }
                }

                if (found==0 || (radiobutton_config.active==0 && found<multichoice_config.value.length+text_filters_len)) {
                    continue;
                }

            selected_indices.push(i);
            }
        }
        cds_config.selected.indices = selected_indices;
        """)
    button_config.js_on_click(cb_button_config)

    return table_config, widgets_config

def define_markers(n_elements):

    # Remove markes without fill-color
    #markers = list(MarkerType)
    #for m in ["asterisk", "cross", "dash", "dot", "x", "y"]:
    #    markers.remove(m)

    # Markers in a "meaningful" order, considering their differences
    markers = ["circle", "plus", "square", "triangle", "star", "diamond", "hex", "inverted_triangle", "circle_cross", "square_cross", "triangle_dot", "star_dot", "diamond_cross", "hex_dot", "circle_dot", "diamond_dot", "square_dot", "triangle_pin", "circle_x", "square_pin", "circle_y", "square_x"]
    if n_elements > len(markers):
        markers = markers * math.ceil(n_elements/len(markers))

    smarkers = markers[:n_elements] # random.sample(markers, n_elements)

    # Not enough "different colors", always keep first 20 very different
    if n_elements <= 10:
        scolor = make_color_palette(n_elements)
    else:
        scolor = list(make_color_palette(10))
        scolor.extend(list(random.sample(make_color_palette(n_elements-10, palette=Turbo256), n_elements-10)))

    return smarkers, scolor

def make_hover(cds_config, exclude_list: list = []):
    hover_tool = HoverTool(
        tooltips=[(f, "@config{" + f + "}") for f in cds_config.data.keys() if f not in exclude_list] +
                 [("value", "@value")],
        formatters={"@config": CustomJSHover(args=dict(
            cds_config=cds_config), code="return cds_config.data[format][value]")}
    )
    return hover_tool


def CustomJSTransform_group_x(cds_config, multiselect_groups):
    # Change values on x-axis based on selected configuration groups
    group_x = CustomJSTransform(
        args=dict(cds_config=cds_config,
                  multiselect_groups=multiselect_groups),
        v_func="""
        if(multiselect_groups.value.length>0){
            const y = new Array(xs.length);
            for (let i = 0; i < xs.length; i++) {
                var group = "|";
                for(let v = 0; v < multiselect_groups.value.length; v++){
                    group += cds_config.data[multiselect_groups.value[v]][xs[i]] + "|"
                }
                y[i] = group;
            }
            return y;
        }else{
            return xs;
        }
    """)
    return group_x

def CustomJSTransform_get_markercolor(slist, cds_config, multiselect):
    get_marker = CustomJSTransform(args=dict(slist=slist,
                                             cds_config=cds_config,
                                             multiselect=multiselect),
                                   v_func="""

                                   var groups_dict = {};
                                   var markers_dict = {};
                                   var cnt_groups = 0;
                                   for(let x = 0; x < xs.length; x++){
                                       var group = "|";
                                       for(let v = 0; v < multiselect.value.length; v++){
                                            group += cds_config.data[multiselect.value[v]][x] + "|"
                                       }
                                       if(!(group in groups_dict)){
                                           groups_dict[group] = cnt_groups;
                                           cnt_groups++;
                                       }
                                       markers_dict[x] = slist[groups_dict[group]];
                                   }

                                   return xs.map(function(x) { return markers_dict[x]; });
                                   """)
    return get_marker


def main_plot(tools, cds_config, cds_target, view_target, smarkers, scolor, multiselect_groups, multiselect_markers, multiselect_colors, init_metric, font_size):

    # BOXPLOT
    #        _________
    #       |     |   |
    #  |----|     |   |--|
    #       |_____|___|
    # lower q1   q2   q3 upper
    #
    cds_boxplot = ColumnDataSource(dict(index=[], lower=[], q1=[], q2=[], q3=[], upper=[]))

    plot = figure(title="",
                  x_range=FactorRange(factors=list(cds_config.data["name"])),
                  toolbar_location="above",
                  tools=tools,
                  width=1200, height=600)
                  #width=800, height=500)


    r = plot.scatter(x=transform("config", CustomJSTransform_group_x(cds_config, multiselect_groups)),
                        y="value",
                        source=cds_target,
                        view=view_target,
                        size=font_size,
                        marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                        fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                        line_color="black")
    plot.add_tools(make_hover(cds_config, exclude_list=["index", "fixed_arguments"]))

    # Activate hover only for scatter points (no boxplot)
    plot.hover.renderers = [r]

    plot.xaxis.major_label_orientation = "vertical"
    plot.yaxis.axis_label = init_metric


    # LEGEND
    legend_items = []
    for i in cds_config.data["index"]:
        legend_items.append(LegendItem(label=str(i), renderers=[r], index=i, visible=False))
    legend_plot = Legend(items=legend_items)
    plot.add_layout(legend_plot, 'right')

    # Font sizes
    #plot.yaxis.formatter.use_scientific = False
    plot.xaxis.axis_label_text_font_size = str(font_size) + "pt"
    plot.yaxis.axis_label_text_font_size = str(font_size) + "pt"
    plot.yaxis.major_label_text_font_size = str(font_size) + "pt"
    plot.xaxis.major_label_text_font_size = str(font_size) + "pt"
    plot.legend.label_text_font_size = str(font_size) + "pt"
    plot.legend.glyph_height = font_size*2
    plot.legend.glyph_width = font_size*2

    #plot.legend.background_fill_color = "#bfbfbf"
    #plot.background_fill_color = "#bfbfbf"
    #plot.border_fill_color = "#bfbfbf"

    # Boxplot (vbar + whisker)
    plot.vbar("index", 0.7, "q2", "q3", source=cds_boxplot, line_color="black", line_width=2, fill_color=None, line_alpha=0.5)
    plot.vbar("index", 0.7, "q1", "q2", source=cds_boxplot, line_color="black", line_width=2, fill_color=None, line_alpha=0.5)
    w = Whisker(base="index", upper="upper", lower="lower", source=cds_boxplot, line_color="black", line_width=2, line_alpha=0.5)
    w.upper_head.size = w.lower_head.size = 20
    plot.add_layout(w)

    return plot, legend_plot, cds_boxplot

def main_callbacks(cds_config, cds_target, plot_target, view_target, legend_plot_target, cds_boxplot_target, multiselect_groups, multiselect_markers, multiselect_colors, select_sort, toggle_boxplot):

    cb_multiselect_target = CustomJS(
        args=dict(cds_config=cds_config,
                  plot_target=plot_target,
                  multiselect_groups=multiselect_groups,
                  multiselect_markers=multiselect_markers,
                  multiselect_colors=multiselect_colors,
                  select_sort=select_sort,
                  plot_target_legend=plot_target.legend[0]),
        code="""

        // build index with combination of selected params
        var factors_set = new Set();

        var legend_idx = {};
        for(let i = 0; i < cds_config.selected.indices.length; i++){
            var group = "|";
            for(let v = 0; v < multiselect_groups.value.length; v++){
                group += cds_config.data[multiselect_groups.value[v]][cds_config.selected.indices[i]] + "|";
            }
            factors_set.add(group);


            ///// LEGEND /////
            // make groups out of two multiselect (marker and color)
            var group_legend = "|"
            for(let v = 0; v <= multiselect_markers.size; v++){
                if(multiselect_markers.value.indexOf(multiselect_markers.options[v]) > -1 ||
                   multiselect_colors.value.indexOf(multiselect_colors.options[v]) > -1){
                    group_legend += cds_config.data[multiselect_markers.options[v]][cds_config.selected.indices[i]] + "|";
                }
            }
            legend_idx[plot_target_legend.items[cds_config.selected.indices[i]].index] = group_legend;
            ///// LEGEND /////

        }

        ///// LEGEND /////
        // very strange approach, the index is not followed, so I had to trick the plotting
        for(let i = 0; i < cds_config.length; i++){
            plot_target_legend.items[i].visible=false;
        }
        var i = 0;
        var already_plotted = new Array();
        for(var key in legend_idx) {
            if(already_plotted.indexOf(legend_idx[key]) == -1){
                plot_target_legend.items[i].label = legend_idx[key];
                plot_target_legend.items[i].visible = true;
                already_plotted.push(legend_idx[key]);
            }
            i++;
        }
        ///// LEGEND /////

        var factors = [...factors_set];
        var sorted_factors = new Array();
        const index_sort = multiselect_groups.value.indexOf(select_sort.value);
        // if sort by specific param (should be selected)
        if(index_sort > -1){
            var order = new Array();
            for (var i = 0; i < factors.length; ++i){
                order.push(factors[i].split("|")[index_sort+1]);
            }

            var idx = new Array(order.length);
            for (var i = 0; i < idx.length; ++i){
                idx[i] = i;
            }
            idx.sort((a, b) => order[a].localeCompare(order[b]));

            for (var i = 0; i < idx.length; ++i){
                sorted_factors.push(factors[idx[i]]);
            }

        }else{
            // default sort
            sorted_factors = factors.sort();
        }

        plot_target.x_range.factors = sorted_factors;
        """)

    cb_toggle_boxplot = CustomJS(
        args=dict(multiselect_groups=multiselect_groups,
                  view_target=view_target,
                  cds_target=cds_target,
                  cds_config=cds_config,
                  cds_boxplot_target=cds_boxplot_target,
                  toggle_boxplot=toggle_boxplot,
                  plot_target=plot_target),
        code='''
        cds_boxplot_target.data = {"index": [], "lower": [], "q1": [], "q2": [], "q3": [], "upper": []};

        if(!toggle_boxplot.active.includes(0)){
            // return alpha dots
            plot_target.renderers[0].glyph.fill_alpha = 1;
            plot_target.renderers[0].glyph.line_alpha = 1;
            plot_target.renderers[0].glyph.hatch_alpha = 1;
        }else{
            
            plot_target.renderers[0].glyph.fill_alpha = 0.3;
            plot_target.renderers[0].glyph.line_alpha = 0.3;
            plot_target.renderers[0].glyph.hatch_alpha = 0.3;

            // for each entry on cds_eval/cds_bench, get config + values on multiselect to rebuild
            var groups_values = {};
            for(let i = 0; i < view_target._indices.length; i++){
                const idx_target = view_target._indices[i];
                const val = cds_target.data["value"][idx_target];
                
                var group = "|";
                for(let v = 0; v < multiselect_groups.value.length; v++){
                     group += cds_config.data[multiselect_groups.value[v]][cds_target.data["config"][idx_target]] + "|"
                }
                
                if(!(group in groups_values)){
                    groups_values[group] = new Array();
                }
                groups_values[group].push(val);

            }

            // quantile function
            const quantile = (sorted, q) => {
                // requires sorted array
                const pos = (sorted.length - 1) * q;
                const base = Math.floor(pos);
                const rest = pos - base;
                if (sorted[base + 1] !== undefined) {
                    return sorted[base] + rest * (sorted[base + 1] - sorted[base]);
                } else {
                    return sorted[base];
                }
            };

            for(group in groups_values){
                const sorted_arr = groups_values[group].sort((a, b) => a - b);
                cds_boxplot_target.data["index"].push(group);

                const q1 = quantile(sorted_arr, .25);
                const q3 = quantile(sorted_arr, .75);

                var lower = q1 - (1.5 * (q3-q1)); // considering outliers IQR*1.5
                if(lower < quantile(sorted_arr, 0)){
                    lower = quantile(sorted_arr, 0);
                }

                var upper = q3 + (1.5 * (q3-q1)); // considering outliers IQR*1.5
                if(upper > quantile(sorted_arr, 1)){
                    upper = quantile(sorted_arr, 1)
                }

                cds_boxplot_target.data["lower"].push(lower);
                cds_boxplot_target.data["q1"].push(q1);
                cds_boxplot_target.data["q2"].push(quantile(sorted_arr, .50));
                cds_boxplot_target.data["q3"].push(q3);
                cds_boxplot_target.data["upper"].push(upper);
            }

        }

        cds_boxplot_target.change.emit();

        ''')

    cb_toggle_label = CustomJS(
        args=dict(xaxis=plot_target.xaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }
        ''')
    
    cb_toggle_legend = CustomJS(
        args=dict(legend_plot_target=legend_plot_target),
        code='''
        if(this.active.includes(0)){
            legend_plot_target.visible = true;
        }else{
            legend_plot_target.visible = false;
        }
        ''')

    cb_multiselect_markers_color = CustomJS(
        args=dict(cds_target=cds_target),
        code="""
        cds_target.change.emit();      
        """)

    return cb_multiselect_target, cb_multiselect_markers_color, cb_toggle_boxplot, cb_toggle_label, cb_toggle_legend


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
