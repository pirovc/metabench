#!/usr/bin/env python3
import argparse
import glob
import sys
import os
import json

import pandas as pd

# Bokeh
from bokeh.io import save
from bokeh.plotting import figure, output_file
from bokeh.models import ColumnDataSource, CDSView, Select, CustomJS, CustomJSTransform, FactorRange, MultiSelect

from bokeh.models.widgets import DataTable, TableColumn
from bokeh.models.filters import GroupFilter, IndexFilter
from bokeh.palettes import Category10, Category20, Colorblind, linear_palette, Turbo256
    
from bokeh.transform import linear_cmap, transform   
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

    default_ranks = ["root",
                     "superkingdom",
                     "phylum",
                     "class",
                     "order",
                     "family",
                     "genus",
                     "species"]

    if len(sys.argv) == 1:
        parser.print_help()
        return False

    args = parser.parse_args()

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
        tables[rep_type][rep_cate]["config"] =  pd.concat(
                    [tables[rep_type][rep_cate]["config"], load_config(pjson)], ignore_index=True)

        # Need to be to concat by rank
        tables[rep_type][rep_cate]["metrics"] =  pd.concat(
                    [tables[rep_type][rep_cate]["metrics"], load_metrics(pjson)], axis=1, ignore_index=True)
    

    # for t in report_types:
    #     for c in report_categories:
    #         for e in report_elements:
    #             print(t,c,e)
    #             print(tables[t][c][e])


    table = tables["profiling"]["evals"]["config"]
    table["version"] = table["version"] + " (" + table["tool"] + ")"

    cols = [TableColumn(field=field, title=field) for field in table.columns]
    cds_table = ColumnDataSource(table)
    data_table = DataTable(source=cds_table, columns=cols, selectable="checkbox")

    s = pd.DataFrame(tables["profiling"]["evals"]["metrics"].stack().reset_index().set_index("metric"))
    s.rename(columns = {'level_2':'config', 0:'value'}, inplace=True)
    s.sort_values(by=['metric', 'rank', 'config'], inplace=True)

    print(s)
    cds = ColumnDataSource(s)



    metrics = []
    for metric in set(cds.data["metric"]):
        v = metric.split("|")
        metrics.append((metric, v[1] + " (" + v[0] + ")"))
    metrics.sort(key=lambda a: a[1])

    ########################################### p1
    p = figure(title="ranks", x_range=default_ranks)
    config_filter = IndexFilter(indices=[])
    metric_filter = GroupFilter(column_name='metric', group="rank-based|f1_score")
    cds_view = CDSView(source=cds, filters=[metric_filter, config_filter])
    p.scatter(x="rank", y="value", 
        legend_field="config",
        source=cds,
        view=cds_view,
        color=linear_cmap(field_name="config", palette=make_color_palette(len(set(cds.data["config"]))), low=0, high=len(set(cds.data["config"]))))

    metric_select = Select(title="Metric:", value="", options=metrics)

    callback_metrics = CustomJS(
        args=dict(cds=cds, metric_filter=metric_filter),
        code='''
        metric_filter.group = this.value;
        cds.change.emit();
        ''')
    metric_select.js_on_change('value', callback_metrics)
    
    callback_metrics = CustomJS(
        args=dict(cds=cds,config_filter=config_filter),
        code='''
        const indices = [];
        for(let i = 0; i < cds.length; i++){
            if(this.indices.indexOf(cds.data["config"][i]) > -1){
                indices.push(i);
            }
        }
        config_filter.indices = indices;
        cds.change.emit();
        ''') 
    cds_table.selected.js_on_change('indices', callback_metrics)



    ########################################### p2
    rank_filter = GroupFilter(column_name='rank', group="species")
    metric1_select = Select(title="Metric:", value="", options=metrics)
    metric2_select = Select(title="Metric:", value="", options=metrics)
    rank_select = Select(title="Rank:", value="species", options=default_ranks)
    callback_rank = CustomJS(
        args=dict(cds=cds, rank_filter=rank_filter),
        code='''
        rank_filter.group = this.value;
        cds.change.emit();
        ''')
    rank_select.js_on_change('value', callback_rank)


    v_func = """
        var shift = 0;
        for (let i = 0; i < xs.length; i++) {
            if(xs[i]==metric2_select.value){
                shift=i;
                break;
            }
        }

        const y = new Float64Array(xs.length);
        var x = 0;
        for (let i = 0; i < xs.length; i++) {
            if(xs[i]==metric1_select.value){
                y[i] = cds.data["value"][shift+x];
                x++;
            }else{
                y[i] = -1;
            }
        }

        return y;
    """
    normalize = CustomJSTransform(args=dict(cds=cds, rank_select=rank_select, metric1_select=metric1_select, metric2_select=metric2_select), v_func=v_func)

    metric1_filter = GroupFilter(column_name='metric', group="rank-based|f1_score")
    cds_view2 = CDSView(source=cds, filters=[rank_filter, metric1_filter, config_filter])
    p2 = figure(title="compare")
    p2.scatter(x="value", y=transform('metric', normalize), 
        legend_field="config",
        source=cds,
        view=cds_view2,
        color=linear_cmap(field_name="config", palette=make_color_palette(len(set(cds.data["config"]))), low=0, high=len(set(cds.data["config"]))))

    callback_metrics1 = CustomJS(
        args=dict(cds=cds, metric1_filter=metric1_filter),
        code='''
        metric1_filter.group = this.value;
        cds.change.emit();
        ''')
    metric1_select.js_on_change('value', callback_metrics1)

    callback_metrics2 = CustomJS(
        args=dict(cds=cds),
        code='''
        cds.change.emit();
        ''')
    metric2_select.js_on_change('value', callback_metrics2)

 





    ########################################### p3
    p3 = figure(title="group", x_range=FactorRange(factors=list(map(str,cds_table.data["index"]))))

    rank_filter_group = GroupFilter(column_name='rank', group="species")
    metric_filter_group = GroupFilter(column_name='metric', group="rank-based|f1_score")
    cds_view_group = CDSView(source=cds, filters=[rank_filter_group, metric_filter_group, config_filter])

    metric_select_group = Select(title="Metric:", value="rank-based|f1_score", options=metrics)
    rank_select_group = Select(title="Rank:", value="species", options=default_ranks)
    group_multi_select = MultiSelect(value=[], options=table.columns.to_list())

    group_x = CustomJSTransform(
        args=dict(cds=cds, cds_table=cds_table, group_multi_select=group_multi_select), 
        v_func="""
        if(group_multi_select.value.length>0){
            const y = new Array(xs.length);
            for (let i = 0; i < xs.length; i++) {
                var group = "|";
                for(let v = 0; v < group_multi_select.value.length; v++){
                    group += cds_table.data[group_multi_select.value[v]][xs[i]] + "|"
                }
                y[i] = group;
            }
            return y;
        }else{
            return xs;
        }
    """)


    p3.scatter(x=transform('config', group_x), y="value", 
        legend_field="config",
        source=cds,
        view=cds_view_group,
        color=linear_cmap(field_name="config", palette=make_color_palette(len(set(cds.data["config"]))), low=0, high=len(set(cds.data["config"]))))


    callback_rank_group = CustomJS(
        args=dict(cds=cds, rank_filter_group=rank_filter_group),
        code='''
        rank_filter_group.group = this.value;
        cds.change.emit();
        ''')
    rank_select_group.js_on_change('value', callback_rank_group)

    callback_metric_group = CustomJS(
        args=dict(cds=cds, metric_filter_group=metric_filter_group),
        code='''
        metric_filter_group.group = this.value;
        cds.change.emit();
        ''')
    metric_select_group.js_on_change('value', callback_metric_group)


    callback_group = CustomJS(
        args=dict(cds=cds, cds_table=cds_table, p3=p3),
        code='''
        var factors = new Set();

        console.log(cds_table.data["index"]);
        console.log(cds_table.selected.indices);
        for(let i = 0; i < cds_table.selected.indices.length; i++){
            var group = "|";
            for(let v = 0; v < this.value.length; v++){
                group += cds_table.data[this.value[v]][cds_table.selected.indices[i]] + "|"
            }
            factors.add(group);
        }
        var sorted_factors = [...factors].sort();
        p3.x_range.factors = sorted_factors;
        ''')
    group_multi_select.js_on_change('value', callback_group)

    layout = column([data_table, 
                    row([
                        column([p, metric_select]),
                        column([p2, rank_select, metric1_select, metric2_select]),
                        column([p3, rank_select_group, metric_select_group, group_multi_select])]
                        )
                    ])
    output_file(args.output, title="Metabench", mode="inline")
    save(layout)

    return True


def plot_config_datatable(table):
    cols = [TableColumn(field=field, title=field) for field in table.columns]
    cds = ColumnDataSource(table)
    cds.selected.indices = list(range(len(cds.data["index"])))
    return DataTable(source=cds, columns=cols, selectable="checkbox")

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
            cols.append((m1,"",))
            vals.append(v1)

    return pd.DataFrame(vals, index=pd.MultiIndex.from_tuples(cols, names=["metric", "rank"]))



def make_color_palette(n_colors, linear: bool=False, palette: dict=None):
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

if __name__ == "__main__":
    sys.exit(0 if main() else 1)
