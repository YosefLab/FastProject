/**
 * Created by David DeTomaso on 6/24/2015.
 */

function HeatMap(parent)
{
    self = this;
    this.h = 1;  //height of row
    this.w = 4;  //width of column

    this.svg = d3.select(parent).append("svg")
        .attr("width", 600)
        .attr("height", 450);

    //define a color scale using the min and max expression values
    this.colorScale = d3.scale.linear()
        .domain([-2, 0, 2])
        .range(["steelblue", "lightgreen", "lightcoral"]);

    this.data = [];
    this.selected = -1;
    this.selected_links = [];

    this.hover_col = -1;
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

HeatMap.prototype.setHovered = function(hovered_index, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event != event_id) {
        this.last_event = event_id;
        this.hover_col = hovered_index;
        this.svg.selectAll("g").selectAll("rect")
            .classed("heatmap-hover", function (d, i) {
                return i == hovered_index
            });
        this.hovered_links.forEach(function (e, i) {
            e.setHovered(hovered_index, event_id);
        });
    }
};

HeatMap.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
        //self.svg.select(".x.axis").call(self.xAxis);
        //self.svg.select(".y.axis").call(self.yAxis);

        //generate heatmap rows
        var heatmapRow = self.svg.selectAll("g")
            .data(self.data);

        heatmapRow.enter()
            .append("g")
            .attr("transform", function(d, j){
                return "translate(0," + (j*self.h) + ")"});
                
        heatmapRow.exit().remove();

        //generate heatmap columns
        var heatmapRects = heatmapRow
            .selectAll("rect")
            .data(function(d) {
                return d;
            });

        heatmapRects.enter().append("rect")
            .attr('width',self.w)
            .attr('height',self.h)
            .attr('y',0)
            .on("mouseover", function(d,i){self.setHovered(i);});

        self.svg
            .on("mouseleave", function(d,i){self.setHovered(-1);});

        heatmapRects.style('fill',function(d) {
            return self.colorScale(d);})
        .attr('x', function(d,i) {
                return (self.col_order[i] * self.w);})

        heatmapRects.exit().remove();

    }
};
