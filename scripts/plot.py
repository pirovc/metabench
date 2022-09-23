#!/usr/bin/env python3
import argparse
import glob
import sys
import os
import json
import random

import randomname  # pip install randomname https://github.com/beasteers/randomname

import pandas as pd

# Bokeh
from bokeh.io import save
from bokeh.core.enums import MarkerType
from bokeh.plotting import figure, output_file
from bokeh.models import ColumnDataSource, CDSView, Select, CustomJS, CustomJSTransform, FactorRange, MultiSelect, CustomJSHover, HoverTool, Tabs, Panel, MultiChoice, Button, RadioButtonGroup, CheckboxGroup, Spacer
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.filters import GroupFilter, IndexFilter
from bokeh.palettes import Category10, Category20, Colorblind, linear_palette, Turbo256

from bokeh.transform import linear_cmap, transform, factor_mark
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
    report_categories = ["benchmark", "size", "stats", "evals"]
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

        # Need to be to concat by rank
        tables[rep_type][rep_cate]["metrics"] = pd.concat(
            [tables[rep_type][rep_cate]["metrics"], load_metrics(pjson)], axis=1, ignore_index=True)

    # for t in report_types:
    #     for c in report_categories:
    #         for e in report_elements:
    #             print(t,c,e)
    #             print(tables[t][c][e])

    main_tabs = []
    building_tabs = []
    profiling_tabs = []
    binning_tabs = []
    if not tables["profiling"]["evals"]["config"].empty:
        profiling_tabs.append(Panel(child=plot_evals(
            "profiling", tables, default_ranks), title="Evaluations"))

    if not tables["binning"]["evals"]["config"].empty:
        binning_tabs.append(Panel(child=plot_evals(
            "binning", tables, default_ranks), title="Evaluations"))

    if not tables["profiling"]["stats"]["config"].empty:
            profiling_tabs.append(Panel(child=plot_stats(
                "profiling", tables, default_ranks), title="Stats"))

    if profiling_tabs:
        main_tabs.append(
            Panel(child=Tabs(tabs=profiling_tabs), title="Profiling"))
    if binning_tabs:
        main_tabs.append(Panel(child=Tabs(tabs=binning_tabs), title="Binning"))

    main_layout = Tabs(tabs=main_tabs)

    output_file(args.output, title="Metabench", mode="inline")
    save(main_layout)

    return True

def plot_stats(report, tables, default_ranks):

    df_config = parse_df_config(tables[report]["stats"]["config"])
    cds_config = ColumnDataSource(df_config)

    df_stats = parse_df_data(tables[report]["stats"]["metrics"])
    cds_stats = ColumnDataSource(df_stats)

    print(df_config)
    print(df_stats)

    #
    # DataTable
    #
    table_config, filter_config, widgets_config = plot_datatable(cds_config, df_config, cds_stats)


    # #
    # # Stats
    # #
    # plot_stats = figure(title="Stats",
    #                     x_range=FactorRange(
    #                          factors=cds_config.data["name"]),
    #                     y_range=Range1d(start=0, end=100)

    # # Change values on x-axis based on selected configuration groups
    # name_x = CustomJSTransform(
    #     args=dict(cds_config=cds_config),
    #     v_func="""
    #     const y = new Array(xs.length);
    #     for (let i = 0; i < xs.length; i++) {
    #         y[i] = cds_config.data["name"][xs[i]];
    #     }
    #     return y;
    # """)

    # plot_stats.scatter(default_ranks, x=transform("config", name_x),
    #                       source=cds_stats,
    #                       width=1,
    #                       #line_color=None,  # to avoid printing small border for zeros
    #                       #color=make_color_palette(top_obs_bars, linear=True) + ("#868b8e", "#eeede7"))
    #                       color=make_color_palette(len(default_ranks)))


    # plot_groups.scatter(x=transform("config", name_x), y="value",
    #                     source=cds_evals,
    #                     view=view_groups,
    #                     size=12,
    #                     marker=transform("config", CustomJSTransform(
    #                         args=dict(e=smarkers), v_func="return xs.map(function(x) { return e[x]; });")),
    #                     fill_color=transform("config", CustomJSTransform(args=dict(
    #                         e=sfillcolor), v_func="return xs.map(function(x) { return e[x]; });")),
    #                     line_color=transform("config", CustomJSTransform(args=dict(e=slinecolor), v_func="return xs.map(function(x) { return e[x]; });")))

    layout = column([row(table_config, column(*widgets_config)),
                     row(Spacer())
                     ])


    return layout

def plot_evals(report, tables, default_ranks):

    df_config = parse_df_config(tables[report]["evals"]["config"])
    cds_config = ColumnDataSource(df_config)

    df_evals = parse_df_data(tables[report]["evals"]["metrics"])
    cds_evals = ColumnDataSource(df_evals)

    #print(df_config)
    #print(df_evals)

    # General variables
    tools = "pan,wheel_zoom,box_zoom,box_select,tap,reset,save"

    metrics = []
    for metric in set(cds_evals.data["metric"]):
        v = metric.split("|")
        metrics.append((metric, v[1] + " (" + v[0] + ")"))
    metrics.sort(key=lambda a: a[1])
    init_metric = next(iter(metrics))[0]
    init_rank = "species" if "species" in default_ranks else default_ranks[-1]

    # Remove markes without fill-color
    markers = list(MarkerType)
    for m in ["asterisk", "cross", "dash", "dot", "x", "y"]:
        markers.remove(m)
    n_config = df_config.shape[0]
    if n_config > len(markers):
        smarkers = random.choices(markers, k=n_config)
    else:
        smarkers = random.sample(markers, n_config)
    sfillcolor = random.sample(make_color_palette(n_config), n_config)
    slinecolor = random.sample(make_color_palette(n_config), n_config)

    # Hover tool for tooltips (all plots)
    hover_tool = HoverTool(
        tooltips=[("Tool", "@config{tool}"),
                  ("Version", "@config{version}"),
                  ("Sample", "@config{sample}"),
                  ("DB", "@config{database}"),
                  ("DB Args.", "@config{database_arguments}"),
                  ("Args. ", "@config{arguments}"),
                  ("Name", "@config{name}")],
        formatters={"@config": CustomJSHover(args=dict(
            cds_config=cds_config), code="return cds_config.data[format][value]")}
    )

    #
    # DataTable
    #
    table_config, filter_config, widgets_config = plot_datatable(cds_config, df_config, cds_evals)

    #
    # Plot Ranks
    #
    plot_ranks = figure(title="Ranks", x_range=default_ranks,
                        toolbar_location="above", tools=tools)

    select_metric_ranks = Select(
        title="Metric:", value=init_metric, options=metrics)
    checkbox_ranks = CheckboxGroup(
        labels=default_ranks, active=list(range(len(default_ranks))))

    filter_metric_ranks = GroupFilter(column_name="metric", group=init_metric)
    view_ranks = CDSView(source=cds_evals, filters=[
        filter_metric_ranks, filter_config])

    plot_ranks.scatter(x="rank", y="value",
                       source=cds_evals,
                       view=view_ranks,
                       size=12,
                       marker=transform("config", CustomJSTransform(
                           args=dict(e=smarkers), v_func="return xs.map(function(x) { return e[x]; });")),
                       fill_color=transform("config", CustomJSTransform(args=dict(
                           e=sfillcolor), v_func="return xs.map(function(x) { return e[x]; });")),
                       line_color=transform("config", CustomJSTransform(args=dict(e=slinecolor), v_func="return xs.map(function(x) { return e[x]; });")))
    plot_ranks.add_tools(hover_tool)
    plot_ranks.yaxis.axis_label = init_metric

    cb_select_metric_ranks = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_metric_ranks=filter_metric_ranks,
                  yaxis=plot_ranks.yaxis[0]),
        code="""
        filter_metric_ranks.group = this.value;
        yaxis.axis_label = this.value;
        cds_evals.change.emit();
        """)
    select_metric_ranks.js_on_change('value', cb_select_metric_ranks)

    cb_checkbox_ranks = CustomJS(
        args=dict(plot_ranks=plot_ranks,
                  default_ranks=default_ranks),
        code="""
        plot_ranks.x_range.factors = this.active.map(function(x) { return default_ranks[x]; });;
        """)
    checkbox_ranks.js_on_click(cb_checkbox_ranks)

    #
    # Plot Compare
    #
    plot_compare = figure(
        title="Compare", toolbar_location="above", tools=tools)

    select_metric_x_compare = Select(
        title="Metric (x):", value=init_metric, options=metrics)
    select_metric_y_compare = Select(
        title="Metric (y):", value=init_metric, options=metrics)
    select_rank_compare = Select(
        title="Rank:", value=init_rank, options=default_ranks)

    filter_rank_compare = GroupFilter(column_name="rank", group=init_rank)
    filter_metric_compare = GroupFilter(
        column_name="metric", group=init_metric)
    view_compare = CDSView(source=cds_evals, filters=[
        filter_rank_compare, filter_metric_compare, filter_config])

    # Select metric for y-axis based on the stacked cds (assumes data is sorted)
    get_metric_y = CustomJSTransform(args=dict(cds_evals=cds_evals,
                                               select_rank_compare=select_rank_compare,
                                               select_metric_x_compare=select_metric_x_compare,
                                               select_metric_y_compare=select_metric_y_compare),
                                     v_func="""
        // find start position of the y-axis metric
        var st_pos = 0;
        for (let i = 0; i < xs.length; i++) {
            if(xs[i]==select_metric_y_compare.value){
                st_pos=i;
                break;
            }
        }

        // Return values starting from position st_pos (assumes data is sorted)
        const y = new Float64Array(xs.length);
        var x = 0;
        for (let i = 0; i < xs.length; i++) {
            if(xs[i]==select_metric_x_compare.value){
                y[i] = cds_evals.data["value"][st_pos+x];
                x++;
            }else{
                y[i] = -1;
            }
        }
        return y;
    """)

    plot_compare.scatter(x="value",
                         y=transform("metric", get_metric_y),
                         source=cds_evals,
                         view=view_compare,
                         size=12,
                         marker=transform("config", CustomJSTransform(
                             args=dict(e=smarkers), v_func="return xs.map(function(x) { return e[x]; });")),
                         fill_color=transform("config", CustomJSTransform(args=dict(
                             e=sfillcolor), v_func="return xs.map(function(x) { return e[x]; });")),
                         line_color=transform("config", CustomJSTransform(args=dict(e=slinecolor), v_func="return xs.map(function(x) { return e[x]; });")))
    plot_compare.add_tools(hover_tool)
    plot_compare.xaxis.axis_label = init_metric
    plot_compare.yaxis.axis_label = init_metric

    cb_select_rank_compare = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_rank_compare=filter_rank_compare),
        code="""
        filter_rank_compare.group = this.value;
        cds_evals.change.emit();
        """)
    select_rank_compare.js_on_change('value', cb_select_rank_compare)

    cb_select_metric_x_compare = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_metric_compare=filter_metric_compare,
                  yaxis=plot_compare.yaxis[0]),
        code="""
        filter_metric_compare.group = this.value;
        yaxis.axis_label = this.value;
        cds_evals.change.emit();
        """)
    select_metric_x_compare.js_on_change("value", cb_select_metric_x_compare)

    # Only triggers cds to apply transform on y-axis
    cb_select_metric_y_compare = CustomJS(
        args=dict(cds_evals=cds_evals,
                  xaxis=plot_compare.xaxis[0]),
        code="""
        xaxis.axis_label = this.value;
        cds_evals.change.emit();
        """)
    select_metric_y_compare.js_on_change('value', cb_select_metric_y_compare)

    #
    # Plot Group
    #
    plot_groups = figure(title="Groups",
                         x_range=FactorRange(
                             factors=list(cds_config.data["name"])),
                         toolbar_location="above",
                         tools=tools)

    select_metric_groups = Select(
        title="Metric:", value=init_metric, options=metrics)
    select_rank_groups = Select(
        title="Rank:", value=init_rank, options=default_ranks)
    multiselect_groups = MultiSelect(
        value=["name"], options=df_config.columns.to_list())

    filter_rank_groups = GroupFilter(column_name="rank", group=init_rank)
    filter_metric_groups = GroupFilter(column_name="metric", group=init_metric)
    view_groups = CDSView(source=cds_evals, filters=[
        filter_rank_groups, filter_metric_groups, filter_config])

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

    plot_groups.scatter(x=transform("config", group_x), y="value",
                        source=cds_evals,
                        view=view_groups,
                        size=12,
                        marker=transform("config", CustomJSTransform(
                            args=dict(e=smarkers), v_func="return xs.map(function(x) { return e[x]; });")),
                        fill_color=transform("config", CustomJSTransform(args=dict(
                            e=sfillcolor), v_func="return xs.map(function(x) { return e[x]; });")),
                        line_color=transform("config", CustomJSTransform(args=dict(e=slinecolor), v_func="return xs.map(function(x) { return e[x]; });")))
    plot_groups.add_tools(hover_tool)
    plot_groups.xaxis.major_label_orientation = "vertical"
    plot_groups.yaxis.axis_label = init_metric

    cb_select_rank_groups = CustomJS(
        args=dict(cds_evals=cds_evals, filter_rank_groups=filter_rank_groups),
        code="""
        filter_rank_groups.group = this.value;
        cds_evals.change.emit();
        """)
    select_rank_groups.js_on_change('value', cb_select_rank_groups)

    cb_select_metric_groups = CustomJS(
        args=dict(cds_evals=cds_evals,
                  filter_metric_groups=filter_metric_groups,
                  yaxis=plot_groups.yaxis[0]),
        code="""
        filter_metric_groups.group = this.value;
        yaxis.axis_label = this.value;
        cds_evals.change.emit();
        """)
    select_metric_groups.js_on_change('value', cb_select_metric_groups)

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
                     row([
                         column([plot_ranks, select_metric_ranks, checkbox_ranks]),
                         column(
                             [plot_compare, select_rank_compare, select_metric_x_compare, select_metric_y_compare]),
                         column([plot_groups, select_rank_groups, select_metric_groups, multiselect_groups])]
                         )
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

def parse_df_config(df):
    # Main dataframes and cds
    df_config = pd.DataFrame(df)
    # Add toll name to version (in case of same version number among tools)
    df_config["version"] = df_config["version"] + \
        " (" + df_config["tool"] + ")"
    # Create random names for each configuration
    df_config["name"] = [randomname.get_name()
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
                             height=300)
    

    choices = []
    for h in df_config.columns:
        for u in df_config[h].unique():
            if u:
                choices.append((h+"|"+u, u + " [" + h + "]"))

    multichoice_config = MultiChoice(title="Select", options=choices)
    radiobutton_config = RadioButtonGroup(labels=["AND", "OR"], active=0)
    button_config = Button(label="Apply", button_type="default")
    widgets_config.extend([multichoice_config, radiobutton_config, button_config])

    filter_config = IndexFilter(indices=[])

    cb_button_config = CustomJS(
        args=dict(cds_config=cds_config,
                  multichoice_config=multichoice_config,
                  radiobutton_config=radiobutton_config),
        code="""
        var selected_indices = [];
        for (var i = 0; i < cds_config.length; i++) {

            if (multichoice_config.value.length > 0 ){
                var found = 0;
                for (var m=0; m < multichoice_config.value.length; ++m){
                    const val = multichoice_config.value[m].split("|");
                    if(cds_config.data[val[0]][i]==val[1]){
                        found++;
                    }
                }
                if (found==0 || (radiobutton_config.active==0 && found<multichoice_config.value.length)) {
                    continue;
                }
            }
            selected_indices.push(i);
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

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
