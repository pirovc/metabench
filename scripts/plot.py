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
from bokeh.models import ColumnDataSource, CDSView, Select, CustomJS, CustomJSTransform,CustomJSFilter,  FactorRange, MultiSelect, CustomJSHover, HoverTool, Tabs, Panel, MultiChoice, Button, RadioButtonGroup, CheckboxGroup, CheckboxButtonGroup, Spacer, Range1d, TextInput
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
    report_categories = ["benchmark", "evals"]
    report_elements = ["config", "metrics"]
    tables = {}  # tables[type][category][element]
    for t in report_types:
        tables[t] = {}
        for c in report_categories:
            tables[t][c] = {}
            for e in report_elements:
                tables[t][c][e] = pd.DataFrame()


    for pjson in parsed_json.values():
        rep_type = pjson["report"]
        rep_cate = pjson["category"]
        tables[rep_type][rep_cate]["config"] = pd.concat(
            [tables[rep_type][rep_cate]["config"], load_config(pjson)], ignore_index=True)

        # concat by rank and get same index
        tables[rep_type][rep_cate]["metrics"] = pd.concat(
            [tables[rep_type][rep_cate]["metrics"], load_metrics(pjson)], axis=1, ignore_index=True)

    # for t in report_types:
    #     for c in report_categories:
    #         for e in report_elements:
    #             print(t,c,e)
    #             print(tables[t][c][e])


    # Make unique random names (same as number of config files)
    rnd_names = set()
    while len(rnd_names) < len(parsed_json):
        rnd_names.add(randomname.get_name())

    main_tabs = []
    evals_tabs = []
    bench_tabs = []

    tools = "pan,wheel_zoom,box_zoom,box_select,tap,reset,save"

    if not tables["profiling"]["evals"]["config"].empty:
        evals_tabs.append(Panel(child=plot_evals(
            "profiling", tables, default_ranks, tools, rnd_names), title="Profiling"))

    if not tables["binning"]["evals"]["config"].empty:
        evals_tabs.append(Panel(child=plot_evals(
            "binning", tables, default_ranks, tools, rnd_names), title="Binning"))


    if not tables["profiling"]["benchmark"]["config"].empty:
        bench_tabs.append(Panel(child=plot_bench(
            "profiling", tables, default_ranks, tools, rnd_names), title="Profiling"))

    if not tables["binning"]["benchmark"]["config"].empty:
        bench_tabs.append(Panel(child=plot_bench(
            "binning", tables, default_ranks, tools, rnd_names), title="Binning"))
    
    if not tables["build"]["benchmark"]["config"].empty:
        bench_tabs.append(Panel(child=plot_bench(
            "build", tables, default_ranks, tools, rnd_names), title="Build"))


    if evals_tabs:
        main_tabs.append(
            Panel(child=Tabs(tabs=evals_tabs), title="Evaluations"))
    
    if bench_tabs:
        main_tabs.append(Panel(child=Tabs(tabs=bench_tabs), title="Benchamarks"))

    main_layout = Tabs(tabs=main_tabs)

    output_file(args.output, title="Metabench", mode="inline")
    save(main_layout)

    return True


def plot_bench(report, tables, default_ranks, tools, rnd_names):
    df_config = parse_df_config(tables[report]["benchmark"]["config"], rnd_names)
    cds_config = ColumnDataSource(df_config)

    df_bench = parse_df_data(tables[report]["benchmark"]["metrics"])
    cds_bench = ColumnDataSource(df_bench)

    metrics = sorted(set(cds_bench.data["metric"]))
    init_metric = metrics[0]

    smarkers, scolor = define_markers(df_config.shape[0])

    hover_tool = make_hover(cds_config, ["index", "fixed_arguments"])
    #
    # DataTable
    #
    table_config, filter_config, widgets_config = plot_datatable(cds_config, df_config, cds_bench)
    
    # Marker and Color multiselect
    multiselect_groups = MultiSelect(title="Group by",
        value=["name"], options=df_config.columns.to_list(), size=8)
    multiselect_markers = MultiSelect(title="Marker",
        value=["name"], options=df_config.columns.to_list(), size=8)
    multiselect_colors = MultiSelect(title="Color",
        value=["name"], options=df_config.columns.to_list(), size=8)

    #
    # Plot Group
    #
    plot_groups = figure(title="Groups",
                         x_range=FactorRange(
                             factors=list(cds_config.data["name"])),
                         toolbar_location="above",
                         tools=tools)

    select_metric_x_groups = Select(
        title="Metric (x):", value=init_metric, options=metrics)



    filter_metric_x_groups = GroupFilter(column_name="metric", group=init_metric)
    view_groups = CDSView(source=cds_bench, filters=[filter_metric_x_groups, filter_config])

    plot_groups.scatter(x=transform("config", CustomJSTransform_group_x(cds_config, multiselect_groups)),
                        y="value",
                        source=cds_bench,
                        view=view_groups,
                        size=12,
                        marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                        fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                        line_color="black")
    plot_groups.add_tools(hover_tool)

    plot_groups.xaxis.major_label_orientation = "vertical"
    plot_groups.yaxis.axis_label = init_metric

    cb_multiselect_markers_color = CustomJS(
        args=dict(cds_bench=cds_bench),
        code="""
        cds_bench.change.emit();
        """)
    multiselect_markers.js_on_change('value', cb_multiselect_markers_color)
    multiselect_colors.js_on_change('value', cb_multiselect_markers_color)

    cb_select_metric_x_groups = CustomJS(
        args=dict(cds_bench=cds_bench,
                  filter_metric_x_groups=filter_metric_x_groups,
                  yaxis=plot_groups.yaxis[0]),
        code="""
        filter_metric_x_groups.group = this.value;
        yaxis.axis_label = this.value;
        cds_bench.change.emit();
        """)
    select_metric_x_groups.js_on_change('value', cb_select_metric_x_groups)

    cb_multiselect_groups = CustomJS(
        args=dict(cds_config=cds_config,
                  plot_groups=plot_groups,
                  multiselect_groups=multiselect_groups),
        code="""
        var factors = new Set();
        for(let i = 0; i < cds_config.selected.indices.length; i++){
            var group = "|";
            for(let v = 0; v < multiselect_groups.value.length; v++){
                group += cds_config.data[multiselect_groups.value[v]][cds_config.selected.indices[i]] + "|"
            }
            factors.add(group);
        }
        var sorted_factors = [...factors].sort();
        plot_groups.x_range.factors = sorted_factors;
        """)
    multiselect_groups.js_on_change('value', cb_multiselect_groups)
    # trigger changes on checkboxes table
    cds_config.selected.js_on_change('indices', cb_multiselect_groups)

    layout = column([row(table_config, column(*widgets_config)),
                     row(plot_groups, column([select_metric_x_groups,
                                              multiselect_groups,
                                              multiselect_markers,
                                              multiselect_colors]))
                     ])
    return layout

def plot_evals(report, tables, default_ranks, tools, rnd_names):

    df_config = parse_df_config(tables[report]["evals"]["config"], rnd_names)
    cds_config = ColumnDataSource(df_config)

    df_evals = parse_df_data(tables[report]["evals"]["metrics"])
    cds_evals = ColumnDataSource(df_evals)

    metrics = []
    for metric in set(cds_evals.data["metric"]):
        v = metric.split("|")
        metrics.append((metric, v[1] + " (" + v[0] + ")"))
    metrics.sort(key=lambda a: a[1])
    init_metric = next(iter(metrics))[0]
    init_rank = "species" if "species" in default_ranks else default_ranks[-1]

    smarkers, scolor = define_markers(df_config.shape[0])

    # Hover tool for tooltips (all plots)
    hover_tool = make_hover(cds_config, ["index"])

    # General widgets

    select_metric = Select(
        title="Metric:", value=init_metric, options=metrics)

    select_rank = Select(
        title="Rank:", value=init_rank, options=default_ranks)

    radio_ranges = RadioButtonGroup(labels=["min-max", "0-1", "0-max"], active=0)

    #
    # DataTable
    #
    table_config, filter_config, widgets_config = plot_datatable(cds_config, df_config, cds_evals)

    # Marker and Color multiselect
    multiselect_groups = MultiSelect(title="Group by",
        value=["name"], options=df_config.columns.to_list(), size=8)
    toggle_label_groups = CheckboxGroup(labels=["Hide samples labels"], active=[])
    sort_groups = Select(title="Sort y:", value="name", options=df_config.columns.to_list())
    multiselect_markers = MultiSelect(title="Marker",
        value=["name"], options=df_config.columns.to_list(), size=8)
    multiselect_colors = MultiSelect(title="Color",
        value=["name"], options=df_config.columns.to_list(), size=8)

    #
    # Plot Ranks
    #
    plot_ranks = figure(title="Ranks", x_range=default_ranks,
                        toolbar_location="above", tools=tools,
                        width=400, height=400)

    checkbox_ranks = CheckboxGroup(
        labels=default_ranks, active=list(range(len(default_ranks))))

    filter_metric_ranks = GroupFilter(column_name="metric", group=init_metric)
    view_ranks = CDSView(source=cds_evals, filters=[
        filter_metric_ranks, filter_config])

    plot_ranks.scatter(x="rank", y="value",
                       source=cds_evals,
                       view=view_ranks,
                       size=12,
                       marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                       fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                       line_color="black")
    plot_ranks.add_tools(hover_tool)
    plot_ranks.yaxis.axis_label = init_metric
    plot_ranks.xaxis.major_label_orientation = "vertical"

    cb_checkbox_ranks = CustomJS(
        args=dict(plot_ranks=plot_ranks,
                  default_ranks=default_ranks),
        code="""
        plot_ranks.x_range.factors = this.active.map(function(x) { return default_ranks[x]; });
        """)
    checkbox_ranks.js_on_click(cb_checkbox_ranks)

    #
    # Plot Compare
    #
    plot_compare = figure(
        title="Compare", toolbar_location="above", tools=tools,
        width=400, height=400)

    select_metric_x_compare = Select(
        title="Metric (x):", value=init_metric, options=metrics)
    
    filter_rank_compare = GroupFilter(column_name="rank", group=init_rank)
    filter_metric_compare = GroupFilter(
        column_name="metric", group=init_metric)
    view_compare = CDSView(source=cds_evals, filters=[
        filter_rank_compare, filter_metric_compare, filter_config])

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

    plot_compare.add_tools(hover_tool)
    plot_compare.xaxis.axis_label = init_metric
    plot_compare.yaxis.axis_label = init_metric

    # Only triggers cds to apply transform on y-axis
    cb_select_metric_x_compare = CustomJS(
        args=dict(cds_evals=cds_evals,
                  xaxis=plot_compare.xaxis[0]),
        code="""
        xaxis.axis_label = this.value;
        cds_evals.change.emit();
        """)
    select_metric_x_compare.js_on_change('value', cb_select_metric_x_compare)

    #
    # Plot Group
    #
    plot_groups = figure(title="Groups",
                         x_range=FactorRange(
                             factors=list(cds_config.data["name"])),
                         #y_range=Range1d(start=0),
                         toolbar_location="above",
                         tools=tools,
                         width=800, height=800)



    filter_rank_groups = GroupFilter(column_name="rank", group=init_rank)
    filter_metric_groups = GroupFilter(column_name="metric", group=init_metric)
    view_groups = CDSView(source=cds_evals, filters=[
        filter_rank_groups, filter_metric_groups, filter_config])

    plot_groups.scatter(x=transform("config", CustomJSTransform_group_x(cds_config, multiselect_groups)),
                        y="value",
                        source=cds_evals,
                        view=view_groups,
                        size=12,
                        marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                        fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                        line_color="black")
    plot_groups.add_tools(hover_tool)

    plot_groups.xaxis.major_label_orientation = "vertical"
    plot_groups.yaxis.axis_label = init_metric

    cb_multiselect_markers_color = CustomJS(
        args=dict(cds_evals=cds_evals),
        code="""
        cds_evals.change.emit();
        """)
    multiselect_markers.js_on_change('value', cb_multiselect_markers_color)
    multiselect_colors.js_on_change('value', cb_multiselect_markers_color)

    cb_multiselect_groups = CustomJS(
        args=dict(cds_config=cds_config,
                  plot_groups=plot_groups,
                  multiselect_groups=multiselect_groups,
                  sort_groups=sort_groups),
        code="""

        // build index with combination of selected params
        var factors_set = new Set();
        for(let i = 0; i < cds_config.selected.indices.length; i++){
            var group = "|";
            for(let v = 0; v < multiselect_groups.value.length; v++){
                group += cds_config.data[multiselect_groups.value[v]][cds_config.selected.indices[i]] + "|"
            }
            factors_set.add(group);
        }
        var factors = [...factors_set];

        
        var sorted_factors = new Array();
        const index_sort = multiselect_groups.value.indexOf(sort_groups.value);
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

        plot_groups.x_range.factors = sorted_factors;
        """)
    multiselect_groups.js_on_change('value', cb_multiselect_groups)
    sort_groups.js_on_change('value', cb_multiselect_groups)

    cb_toggle_label_groups = CustomJS(
        args=dict(xaxis=plot_groups.xaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }else{
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
        }
        ''')
    toggle_label_groups.js_on_click(cb_toggle_label_groups)

    # trigger changes on checkboxes table
    cds_config.selected.js_on_change('indices', cb_multiselect_groups)

    cb_radio_ranges = CustomJS(
        args=dict(radio_ranges=radio_ranges,
                  cds_evals=cds_evals,
                  y_range_groups=plot_groups.y_range,
                  y_range_ranks=plot_ranks.y_range,
                  y_range_compare=plot_compare.y_range,
                  x_range_compare=plot_compare.x_range),
        code="""
        if (radio_ranges.active==0){
            // auto
            //y_range_groups.start = NaN;
            //y_range_groups.end = NaN;
            cds_evals.change.emit();
        }else if(radio_ranges.active==1){
            // 0-1
            y_range_groups.start = 0;
            y_range_groups.end = 1;

            y_range_ranks.start = 0;
            y_range_ranks.end = 1;

            y_range_compare.start = 0;
            y_range_compare.end = 1;
            x_range_compare.start = 0;
            x_range_compare.end = 1;
        }else if(radio_ranges.active==2){
            // 0-max
            y_range_groups.start = 0;

            y_range_ranks.start = 0;

            y_range_compare.start = 0;
            x_range_compare.start = 0;
        }
        """)
    radio_ranges.js_on_click(cb_radio_ranges)

    cb_select_metric = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_metric_ranks=filter_metric_ranks,
                  filter_metric_compare=filter_metric_compare,
                  filter_metric_groups=filter_metric_groups,
                  yaxis_ranks=plot_ranks.yaxis[0],
                  yaxis_compare=plot_compare.yaxis[0],
                  yaxis_groups=plot_groups.yaxis[0]),

        code="""
        filter_metric_ranks.group = this.value;
        yaxis_ranks.axis_label = this.value;

        filter_metric_compare.group = this.value;
        yaxis_compare.axis_label = this.value;

        filter_metric_groups.group = this.value;
        yaxis_groups.axis_label = this.value;

        cds_evals.change.emit();

        //console.log(Bokeh.documents[0].get_model_by_id('my_select'))
        //radio_ranges.trigger_event(({"event_name": "click"}))
        """)
    select_metric.js_on_change('value', cb_select_metric)

    cb_select_rank = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_rank_compare=filter_rank_compare,
                  filter_rank_groups=filter_rank_groups),
        code="""
        filter_rank_compare.group = this.value;
        filter_rank_groups.group = this.value;
        cds_evals.change.emit();
        """)
    select_rank.js_on_change('value', cb_select_rank)


    layout = column([row(table_config, column(*widgets_config)),
                     row([plot_groups, 
                          column([select_metric,
                                  select_rank,
                                  radio_ranges,
                                  multiselect_groups,
                                  sort_groups,
                                  multiselect_markers,
                                  multiselect_colors,
                                  toggle_label_groups
                                 ]),
                          column(row([plot_compare,select_metric_x_compare]),
                                 row([plot_ranks, checkbox_ranks])),
                    ])])

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


def load_metrics(pjson):
    cols = []
    vals = []
    for m1, v1 in pjson["metrics"].items():
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
    # Stack metrics and reset index
    df_data = pd.DataFrame(
        df.stack().reset_index().set_index("metric"))
    # Rename column headers
    df_data.rename(columns={'level_2': "config", 0: 'value'}, inplace=True)
    # Make sure is sorted to properly use some js methods
    df_data.sort_values(by=["metric", 'rank', "config"], inplace=True)
    return df_data

def plot_datatable(cds_config, df_config, cds_data):
    widgets_config = []
    table_config = DataTable(source=cds_config,
                             columns=[TableColumn(field=field, title=field)
                                      for field in df_config.columns],
                             selectable="checkbox",
                             width=1200,
                             height=150)
    

    choices = []
    for h in df_config.columns:
        for u in df_config[h].unique():
            if u:
                print(h,u)
                choices.append((h+"|"+u, u + " [" + h + "]"))

    multichoice_config = MultiChoice(options=choices, placeholder="Fixed filter")
    textinput_config = TextInput(placeholder="Partial filter e.g. database:abfv|sample:cami")
    radiobutton_config = RadioButtonGroup(labels=["AND", "OR"], active=0)
    button_config = Button(label="Apply", button_type="success")
    widgets_config.extend([multichoice_config, textinput_config, radiobutton_config, button_config])

    filter_config = IndexFilter(indices=[])

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

    cb_cds_config = CustomJS(
        args=dict(cds_data=cds_data, filter_config=filter_config),
        code="""
        const indices = [];
        for(let i = 0; i < cds_data.length; i++){
            if(this.indices.indexOf(cds_data.data["config"][i]) > -1){
                indices.push(i);
            }
        }
        filter_config.indices = indices;
        cds_data.change.emit();
        """)
    cds_config.selected.js_on_change('indices', cb_cds_config)

    return table_config, filter_config, widgets_config

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

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
