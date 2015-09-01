/**
 * Created by David DeTomaso on 6/24/2015.
 */

function HeatMap(parent)
{
    var self = this;
    this.h = 1;  //height of row

    this.width = 500;
    this.height = 450;

    this.svg = d3.select(parent).append("svg")
        .attr("width", self.width)
        .attr("height", self.height);

    var offset = 0;
        
    this.labels = this.svg.append("g");
    
    this.labels.append("text")
        .classed("col_label", true)
        .attr("x", 0)
        .attr("y", offset + 20)
        .attr("font-size", "20px");
        
    this.labels.append("text")
        .classed("row_label", true)
        .attr("x", 0)
        .attr("y", offset + 40)
        .attr("font-size", "20px");

    this.labels.append("rect")
        .attr("x", this.width-40)
        .attr("y", offset)
        .attr("width", 40)
        .attr("height", 40)
        .style("fill", "white");

    this.labels.append("text")
        .classed("rect_label", true)
        .attr("text-anchor", "middle")
        .attr("x", this.width-20)
        .attr("y", offset+23)
        .attr("font-size", "10px");
        
    offset += 40;
    offset += 10; //Some margin
        
    this.grid = this.svg.append("g");
    this.grid
        .attr("transform", "translate(0," +
        (offset)+")");
    
    offset += this.height;
    this.svg.attr("height", offset);

    //define a color scale using the min and max expression values
    this.colorScale = d3.scale.linear()
        .domain([-.6, 0, .6])
        .range(["steelblue", "white", "lightcoral"]);

    this.data = [];
    this.selected = -1;
    this.selected_links = [];

    this.hover_cols = -1;
    this.hovered_links = [];
    
    this.hover_rows = -1;

    this.last_event = 0;

    this.row_labels = [];
    this.col_labels = [];

    this.cluster_assignments = []; //Used in hoverCol.  Denotes which cluster each sample is assigned to

}


HeatMap.prototype.setData = function(data, cluster_assignments)
{
    //Data is an array of rows, each containing an array of values for each col
    //cluster_assignments is an array of numbers indicating assignment

    this.cluster_assignments = cluster_assignments;
    var dataT = d3.transpose(data);
    var TOTAL_SAMPLES = cluster_assignments.length;
    var clusters = {};
    for(var j = 0; j < cluster_assignments.length; j++)
    {
        if(clusters[cluster_assignments[j]] === undefined){
            clusters[cluster_assignments[j]] = {'weight':0, 'data':[]};
        }

        var clust = clusters[cluster_assignments[j]];
        clust.weight = clust.weight + 1;
        clust.data.push(dataT[j]);
    }

    //Now clusters is a dict of "cluster index" -> ('weight' -> # of cols in cluster, 'data' -> N_cols x N_row for cluster)
    for (var clust_i in clusters)
    {
        var clust = clusters[clust_i];
        var cdT = d3.transpose(clust.data);
        clust.data = cdT.map(function(e){return d3.mean(e);});
    }
    //Now, each clusters data is just a length N_genes array of values for that cluster.

    //Convert to a list for easier d3 binding
    cluster_list = [];
    for (var clust_i in clusters)
    {
        var clust = clusters[clust_i];
        clust['index'] = clust_i;
        cluster_list.push(clust);
    }

    //Add in information on the width of each cluster
    var x_offset = 0;
    for (var j = 0; j < cluster_list.length; j++)
    {
        var clust = cluster_list[j];
        var width = clust['weight'] / TOTAL_SAMPLES * this.width;
        
        clust.data = clust.data.map(function(e,i){
            return {"value":e, "x":x_offset, "width":width, "index": clust['index']};
        });

        x_offset = x_offset + width;
    }

    this.data = cluster_list;
    N_ROWS = this.data[0].data.length;

    this.h = Math.floor(this.height/N_ROWS);
    if(this.h == 0) {this.h = 1;}

    this.redraw()();
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

};

HeatMap.prototype.setHoveredCol = function(hovered_col_index)
{
    hovered_indices = [];
    if(hovered_col_index !== undefined) {
        for (var i = 0; i < this.cluster_assignments.length; i++) {
            if (this.cluster_assignments[i] == hovered_col_index) {
                hovered_indices.push(i);
            }
        }
    }

    this.hovered_links.forEach(function (e, i) {
        e.setHovered(hovered_indices);
    });
};

HeatMap.prototype.setHoveredRow = function(hovered_row_indices)
{
    if(typeof(hovered_row_indices) == "number"){
        hovered_row_indices = [hovered_row_indices];
    };
    
    this.hover_rows = hovered_row_indices;
    
    if(this.hover_rows.length == 1)
    {
        var ii = this.hover_rows[0];
        this.labels.select(".row_label").text(this.row_labels[ii]);
    }
    
}

//Sets the text and color for the square upper-right indicator
HeatMap.prototype.setHoveredIndicator = function(data_val)
{
    if(data_val !== undefined)
    {
        this.labels.select("rect")
            .style("fill", this.colorScale(data_val));
        this.labels.select(".rect_label")
            .text(data_val.toFixed(3));
    }
    else
    {
        this.labels.select("rect")
            .style("fill", "white");
        this.labels.select(".rect_label")
            .text("");
    }

}

HeatMap.prototype.redraw = function() {
    var self = this;
    return function(){

        //generate heatmap columns
        var heatmapCols = self.grid.selectAll("g")
            .data(self.data);

        heatmapCols.enter()
            .append("g");

        heatmapCols.exit().remove();

        //generate heatmap rows
        var heatmapRects = heatmapCols
            .selectAll("rect")
            .data(function(d) {
                return d.data;
            });

        heatmapRects.enter().append("rect")
            .on("mouseover", function(d,i){ self.setHoveredRow(i); self.setHoveredCol(d.index); self.setHoveredIndicator(d.value);});
        //    .on("mouseover", function(d){self.setHovered(d.col); self.setHoveredRow(d.row); self.setHoveredIndicator(d.value);});

        self.svg
            .on("mouseleave", function(d){ self.setHoveredRow(-1); self.setHoveredCol(); self.setHoveredIndicator();});
        //    .on("mouseleave", function(d){self.setHovered(-1); self.setHoveredRow(-1); self.setHoveredIndicator();});

        heatmapRects.style('fill',function(d) {
                return self.colorScale(d.value);})
            .attr('width',function(d){return d.width;})
            .attr('height',self.h)
            .attr('y',function(d,i){ return i*self.h;})
            .attr('x', function(d) {
                return d.x;});

        heatmapRects.exit().remove();
    }
};
