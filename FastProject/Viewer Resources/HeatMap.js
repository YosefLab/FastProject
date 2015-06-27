/**
 * Created by David DeTomaso on 6/24/2015.
 */

function HeatMap(parent)
{
    self = this;
    this.h = 1;  //height of row
    this.w = 4;  //width of column

    this.width = 600;
    this.height = 450;

    this.svg = d3.select(parent).append("svg")
        .attr("width", self.width)
        .attr("height", self.height);

    this.grid = this.svg.append("g");
    this.cluster_bar = this.svg.append("g");
    this.cluster_bar.on("mouseleave", function(){
        self.setHovered(-1);
    });

    this.cluster_bar_props = {
        height: 10,
        margin: 10,
        colors: d3.scale.category10().domain(d3.range(10))
    };

    this.grid
        .attr("transform", "translate(0," +
        (self.cluster_bar_props.margin + self.cluster_bar_props.height)+")");

    //define a color scale using the min and max expression values
    this.colorScale = d3.scale.linear()
        .domain([-2, 0, 2])
        .range(["steelblue", "lightgreen", "lightcoral"]);

    this.data = [];
    this.selected = -1;
    this.selected_links = [];

    this.hover_cols = -1;
    this.hovered_links = [];

    this.last_event = 0;

    this.col_clusters = null; //Cluster assignment for each column in data matrix
    this.col_order = null;  //Position for each column in data matrix



}

HeatMap.prototype.cluster_columns = function(assignments)
{
    self = this;
    this.col_clusters = assignments;

    //Argsort the col_clusters
    rr = d3.range(0,this.col_clusters.length);
    rr.sort(function(a,b){return self.col_clusters[a] - self.col_clusters[b];});
    //Argsort again to get rank
    rr2 = d3.range(0,rr.length);
    rr2.sort(function(a,b){return rr[a] - rr[b];});
    
    this.col_order = rr2;
    this.redraw()();
};

HeatMap.prototype.setData = function(data, render)
{
    this.data = data;
    N_ROWS = data.length;
    N_COLS = data[0].length;

    this.w = Math.floor(this.width/N_COLS);
    if(this.w == 0) {this.w = 1;}
    this.h = Math.floor(this.height/N_ROWS);
    if(this.h == 0) {this.h = 1;}

    this.col_clusters = Array.apply(null, Array(N_COLS)).map(Number.prototype.valueOf,0);
    this.col_order = d3.range(0, N_COLS);  //Initial sorting
    
    if(render){
        this.redraw(true)();
    }
};

HeatMap.prototype.setSelected = function(selected_index, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event != event_id) {
        this.last_event = event_id;
        this.selected = selected_index;
        this.redraw()();
        this.selected_links.forEach(function (e, i) {
            e.setSelected(selected_index, event_id);
        });
    }
};

HeatMap.prototype.setHovered = function(hovered_indices, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //test for single index, and wrap in list
    if(typeof(hovered_indices) == "number"){hovered_indices = [hovered_indices];}

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event != event_id) {
        this.last_event = event_id;
        this.hover_cols = hovered_indices;
        this.grid.selectAll("g").selectAll("rect")
            .classed("heatmap-hover", function (d, i) {
                return hovered_indices.indexOf(i) > -1;
            });
        this.hovered_links.forEach(function (e, i) {
            e.setHovered(hovered_indices, event_id);
        });
    }
};

HeatMap.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
        //self.svg.select(".x.axis").call(self.xAxis);
        //self.svg.select(".y.axis").call(self.yAxis);

        //Draw cluster-bar
        var clusterRects = self.cluster_bar.selectAll("rect")
            .data(self.col_clusters);

        clusterRects.enter()
            .append("rect")
            .attr('width', self.w)
            .attr('height', self.cluster_bar_props.height)
            .attr('y', 0);

        clusterRects.style('fill',function(d) {
                return self.cluster_bar_props.colors(d);})
            .attr('x', function(d,i){
                return (self.col_order[i] * self.w); })
            .on("mouseover", function(d) {
            ii = d3.range(self.col_clusters.length);
            selected_i = ii.filter(function(e,j){
                return self.col_clusters[j] == d;});
            self.setHovered(selected_i);
        });

        clusterRects.exit().remove();


        //generate heatmap rows
        var heatmapRow = self.grid.selectAll("g")
            .data(self.data);

        heatmapRow.enter()
            .append("g");

        heatmapRow.attr("transform", function(d, j){
                return "translate(0," + (j*self.h)+ ")"});
                
        heatmapRow.exit().remove();

        //generate heatmap columns
        var heatmapRects = heatmapRow
            .selectAll("rect")
            .data(function(d) {
                return d;
            });

        heatmapRects.enter().append("rect")
            .attr('y',0)
            .on("mouseover", function(d,i){self.setHovered(i);});

        self.svg
            .on("mouseleave", function(d,i){self.setHovered(-1);});

        heatmapRects.style('fill',function(d) {
                return self.colorScale(d);})
            .attr('width',self.w)
            .attr('height',self.h)
            .attr('x', function(d,i) {
                return (self.col_order[i] * self.w);});

        heatmapRects.exit().remove();

    }
};
