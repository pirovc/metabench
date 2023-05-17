def plot_bench(report, tables, default_ranks, tools, rnd_names):
    
    df_config = parse_df_config(tables[report]["benchmark"]["config"], rnd_names)
    cds_config = ColumnDataSource(df_config)

    df_bench = parse_df_data(tables[report]["benchmark"]["metrics"])
    cds_bench = ColumnDataSource(df_bench)

    # BOXPLOT
    #        _________
    #       |     |   |
    #  |----|     |   |--|
    #       |_____|___|
    # lower q1   q2   q3 upper
    #
    cds_boxplot = ColumnDataSource(dict(index=[], lower=[], q1=[], q2=[], q3=[], upper=[]))

    metrics = sorted(set(cds_bench.data["metric"]))
    init_metric = metrics[0]

    smarkers, scolor = define_markers(df_config.shape[0])

    hover_tool = make_hover(cds_config, exclude_list=["index", "fixed_arguments"])

    #
    # Widgets
    #
    # Metric (y)
    select_metric = Select(
        title="Metric (y):", value=init_metric, options=metrics)

    # Group
    multiselect_groups = MultiSelect(title="Group by",
        value=["name"], options=df_config.columns.to_list(), size=8)

    # Marker
    multiselect_markers = MultiSelect(title="Marker",
        value=["name"], options=df_config.columns.to_list(), size=8)

    # Color
    multiselect_colors = MultiSelect(title="Color",
        value=["name"], options=df_config.columns.to_list(), size=8)

    # Sort
    sort_groups = Select(title="Sort (x):", value="name", options=df_config.columns.to_list())
    
    # Toggle options
    toggle_boxplot = CheckboxGroup(labels=["Show Boxplot"], active=[])
    toggle_legend = CheckboxGroup(labels=["Show legend"], active=[0])
    toggle_label = CheckboxGroup(labels=["Show labels (x)"], active=[0])


    #
    # DataTable
    #
    table_config, filter_config, widgets_config = plot_datatable(cds_config, df_config, cds_bench)
    

    #
    # Plot Group
    #
    plot_groups = figure(title="Groups",
                         x_range=FactorRange(factors=list(cds_config.data["name"])),
                         toolbar_location="above",
                         tools=tools,
                         width=1000, height=800)

    filter_metric_groups = GroupFilter(column_name="metric", group=init_metric)
    view_groups = CDSView(source=cds_bench, filters=[filter_metric_groups, filter_config])

    r = plot_groups.scatter(x=transform("config", CustomJSTransform_group_x(cds_config, multiselect_groups)),
                        y="value",
                        source=cds_bench,
                        view=view_groups,
                        size=12,
                        marker=transform("config", CustomJSTransform_get_markercolor(smarkers, cds_config, multiselect_markers)),
                        fill_color=transform("config", CustomJSTransform_get_markercolor(scolor, cds_config, multiselect_colors)),
                        line_color="black")
    plot_groups.add_tools(hover_tool)

    # Activate hover only for scatter points (no boxplot)
    plot_groups.hover.renderers = [r]

    plot_groups.xaxis.major_label_orientation = "vertical"
    plot_groups.yaxis.axis_label = init_metric

    # LEGEND
    legend_items = []
    for i in cds_config.data["index"]:
        legend_items.append(LegendItem(label=str(i), renderers=[r], index=i, visible=False))
    legend_plot_groups = Legend(items=legend_items)
    plot_groups.add_layout(legend_plot_groups, 'right')

    # Boxplot (vbar + whisker)
    plot_groups.vbar("index", 0.7, "q2", "q3", source=cds_boxplot, line_color="black", line_width=2, fill_color=None, line_alpha=0.5)
    plot_groups.vbar("index", 0.7, "q1", "q2", source=cds_boxplot, line_color="black", line_width=2, fill_color=None, line_alpha=0.5)
    w = Whisker(base="index", upper="upper", lower="lower", source=cds_boxplot, line_color="black", line_width=2, line_alpha=0.5)
    w.upper_head.size = w.lower_head.size = 20
    plot_groups.add_layout(w)
    
    cb_toggle_boxplot = CustomJS_toggle_boxplot(multiselect_groups, view_groups, cds_bench, cds_config, cds_boxplot, toggle_boxplot)
    toggle_boxplot.js_on_click(cb_toggle_boxplot)

    cb_multiselect_groups = CustomJS_multiselect(cds_config, plot_groups, multiselect_groups, multiselect_markers, multiselect_colors, sort_groups)
    multiselect_groups.js_on_change('value', cb_multiselect_groups, cb_toggle_boxplot)
    
    cb_multiselect_markers_color = CustomJS(
        args=dict(cds_bench=cds_bench),
        code="""
        cds_bench.change.emit();
        """)
    multiselect_markers.js_on_change('value', cb_multiselect_markers_color, cb_multiselect_groups)
    multiselect_colors.js_on_change('value', cb_multiselect_markers_color, cb_multiselect_groups)

    sort_groups.js_on_change('value', cb_multiselect_groups)

    cb_select_metric = CustomJS(
        args=dict(cds_bench=cds_bench,
                  filter_metric_groups=filter_metric_groups,
                  yaxis=plot_groups.yaxis[0]),
        code="""
        filter_metric_groups.group = this.value;
        yaxis.axis_label = this.value;
        cds_bench.change.emit();
        """)
    select_metric.js_on_change('value', cb_select_metric)

    # trigger changes on checkboxes table
    cds_config.selected.js_on_change('indices', cb_multiselect_groups, cb_toggle_boxplot)

    cb_toggle_label = CustomJS(
        args=dict(xaxis=plot_groups.xaxis[0]),
        code='''
        if(this.active.includes(0)){
            xaxis.major_label_text_font_size = "10px";
            xaxis.major_tick_line_color="black";
        }else{
            xaxis.major_label_text_font_size = "0px";
            xaxis.major_tick_line_color=null;
        }
        ''')
    toggle_label.js_on_click(cb_toggle_label)

    cb_toggle_legend = CustomJS(
        args=dict(legend_plot_groups=legend_plot_groups),
        code='''
        if(this.active.includes(0)){
            legend_plot_groups.visible = true;
        }else{
            legend_plot_groups.visible = false;
        }
        ''')
    toggle_legend.js_on_click(cb_toggle_legend)

    layout = column([row(table_config, column(*widgets_config)),
                     row(plot_groups, column([select_metric,
                                              multiselect_groups,
                                              multiselect_markers,
                                              multiselect_colors,
                                              sort_groups,
                                              toggle_boxplot,
                                              toggle_legend,
                                              toggle_label]))
                     ])
    return layout
