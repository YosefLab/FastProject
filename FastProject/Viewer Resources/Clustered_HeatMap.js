/**
 * Created by David DeTomaso on 6/24/2015.
 */

function HeatMap(parent)
{
    var self = this;
    this.h = 1;  //height of row

    this.width = $(parent).width();
    this.height = $(parent).height();

    var otherHeight = 0;
    //Subtract height of anything else
    $(parent).children().each(function(i,e){ otherHeight += $(e).outerHeight(true);});
    this.height -= otherHeight;

    this.svg = d3.select(parent).append("svg")
        .attr("width", self.width)
        .attr("height", self.height);

    var offset = 0;
        
    this.labels = this.svg.append("g");
        
    this.labels.append("text")
        .classed("row_label", true)
        .attr("x", 0)
        .attr("y", offset + 28)
        .attr("font-size", "20px");

    this.labels.append("rect")
        .attr("x", this.width-40)
        .attr("y", offset)
        .attr("width", 40)
        .attr("height", 30)
        .style("fill", "white");

    this.labels.append("text")
        .classed("rect_label", true)
        .attr("text-anchor", "middle")
        .attr("x", this.width-20)
        .attr("y", offset+17)
        .attr("font-size", "10px");
        
    offset += 30;
    offset += 10; //Some margin
        
    this.heat_height = this.height - offset;
    this.grid_gap = 10; //# of pixels between the plus and minus heat maps
    this.grid_start = offset;
    this.grid_xoffset = 50; //# of pixels to the left of the grids for the + and - labels

    this.grid_plus = this.svg.append("g");
    this.grid_minus = this.svg.append("g");

    this.grid_plus_label = this.svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", this.grid_xoffset/2)
        .attr("y", 0)
        .attr("font-size", "25px")
        .text("");

    this.grid_minus_label = this.svg.append("text")
        .attr("text-anchor", "middle")
        .attr("x", this.grid_xoffset/2)
        .attr("y", 0)
        .attr("font-size", "25px")
        .text("");
    

    //define a color scale using the min and max expression values
    this.colorScale = d3.scale.linear()
        .domain([-0.6, 0, 0.6])
        .range(["steelblue", "white", "lightcoral"]);

    this.data_plus = [];
    this.data_minus = [];
    this.selected = -1;
    this.selected_links = [];

    this.hover_cols = -1;
    this.hovered_links = [];
    
    this.hovered_row = "";

    this.last_event = 0;

    this.row_labels = [];
    this.col_labels = [];

    this.cluster_assignments = []; //Used in hoverCol.  Denotes which cluster each sample is assigned to

}


HeatMap.prototype.setData = function(data, cluster_assignments, gene_labels, gene_signs, sample_labels)
{
    //Data is an array of rows, each containing an array of values for each col
    //cluster_assignments is an array of numbers indicating assignment
    
    this.row_labels = gene_labels;
    this.col_labels = sample_labels;

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
    var cluster_list = [];
    for (var clust_i in clusters)
    {
        var clust = clusters[clust_i];
        clust.index = parseInt(clust_i);
        cluster_list.push(clust);
    }

    //Add in information on the width of each cluster
    var x_offset = this.grid_xoffset;
    for (var j = 0; j < cluster_list.length; j++)
    {
        var clust = cluster_list[j];
        var width = clust.weight / TOTAL_SAMPLES * (this.width - this.grid_xoffset);
        
        clust.data = clust.data.map(function(e,i){
            return {"value":e, "x":x_offset, "width":width, "index": clust.index, "gene": gene_labels[i]};
        });

        x_offset = x_offset + width;
    }

    var N_ROWS = cluster_list[0].data.length;

    this.h = Math.floor((this.heat_height - this.grid_gap)/N_ROWS);
    if(this.h === 0) {this.h = 1;}

    //Split data into data_plus and data_minus
    //Unsigned sigs (sign = 0) go in data_plus

    var cluster_list_plus = [];
    var cluster_list_minus = [];
    for(var i = 0; i < cluster_list.length; i++)
    {
        var clust = cluster_list[i];

        //get positive rows
        var clust_plus = {};
        clust_plus.data = clust.data.filter(function(e,i){
            return gene_signs[i] === 0 || gene_signs[i] === 1;
        });

        //get negative rows
        var clust_minus = {};
        clust_minus.data = clust.data.filter(function(e,i){
            return gene_signs[i] === -1;
        });

        //copy over all other properties
        for(var key in clust){
            if(key !== "data"){
                clust_plus[key] = clust[key];
                clust_minus[key] = clust[key];
            }
        }

        cluster_list_plus.push(clust_plus);
        cluster_list_minus.push(clust_minus);

    }

    this.data_plus = cluster_list_plus;
    this.data_minus = cluster_list_minus;

    this.redraw()();
};

HeatMap.prototype.setSelected = function(selected_index, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event !== event_id) {
        this.last_event = event_id;
        this.selected = selected_index;
        this.redraw()();
        this.selected_links.forEach(function (e) {
            e.setSelected(selected_index, event_id);
        });
    }
};

HeatMap.prototype.setHovered = function(hovered_indices, event_id)
{

};

HeatMap.prototype.setHoveredCol = function(hovered_col_index)
{
    var hovered_indices = [];
    if(hovered_col_index !== undefined) {
        for (var i = 0; i < this.cluster_assignments.length; i++) {
            if (this.cluster_assignments[i] === hovered_col_index) {
                hovered_indices.push(i);
            }
        }
    }

    this.hovered_links.forEach(function (e) {
        e.setHovered(hovered_indices);
    });
};

HeatMap.prototype.setHoveredRowLabel = function(row_label)
{

    this.hovered_row = row_label;
    this.labels.select(".row_label").text(this.hovered_row);
    
};

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

};

HeatMap.prototype.redraw = function() {
    var self = this;
    return function(){

        //generate heatmap columns for plus grid
        
        var pos_grid_start = self.grid_start;
        self.grid_plus.attr("transform",
                "translate(0," + (pos_grid_start)+")"
                );

        var heatmapCols_plus = self.grid_plus.selectAll("g")
            .data(self.data_plus);

        heatmapCols_plus.enter()
            .append("g");

        heatmapCols_plus.exit().remove();

        //generate heatmap rows
        var heatmaRects_plus = heatmapCols_plus
            .selectAll("rect")
            .data(function(d) {
                return d.data;
            });

        heatmaRects_plus.enter().append("rect")
            .on("mouseover", function(d,i){ self.setHoveredRowLabel(d.gene); self.setHoveredCol(d.index); self.setHoveredIndicator(d.value);});

        heatmaRects_plus.style('fill',function(d) {
                return self.colorScale(d.value);})
            .attr('width',function(d){return d.width;})
            .attr('height',self.h)
            .attr('y',function(d,i){ return i*self.h;})
            .attr('x', function(d) {
                return d.x;});

        heatmaRects_plus.exit().remove();

        var num_pos_rects;
        if(self.data_plus.length === 0)
        {
            num_pos_rects = 0;
        } else {
            num_pos_rects = self.data_plus[0].data.length;
        }

        var neg_grid_start = 0;
        neg_grid_start += pos_grid_start; // Start location for all grids
        neg_grid_start += self.h*num_pos_rects; // Offset for space taken up by positive grid
        if(num_pos_rects > 0)
        {
            neg_grid_start += self.grid_gap;  // If there's a positive grid, introduce a gap between them
        }

        self.grid_minus.attr("transform",
                "translate(0," + (neg_grid_start) + ")"
                );

        //generate heatmap columns for minus grid
        var heatmapCols_minus = self.grid_minus.selectAll("g")
            .data(self.data_minus);

        heatmapCols_minus.enter()
            .append("g");

        heatmapCols_minus.exit().remove();

        //generate heatmap rows
        var heatmapRects_minus = heatmapCols_minus
            .selectAll("rect")
            .data(function(d) {
                return d.data;
            });

        heatmapRects_minus.enter().append("rect")
            .on("mouseover", function(d,i){ self.setHoveredRowLabel(d.gene); self.setHoveredCol(d.index); self.setHoveredIndicator(d.value);});

        self.svg
            .on("mouseleave", function(d){ self.setHoveredRowLabel(""); self.setHoveredCol(); self.setHoveredIndicator();});

        heatmapRects_minus.style('fill',function(d) {
                return self.colorScale(d.value);})
            .attr('width',function(d){return d.width;})
            .attr('height',self.h)
            .attr('y',function(d,i){ return i*self.h;})
            .attr('x', function(d) {
                return d.x;});

        heatmapRects_minus.exit().remove();

        var num_neg_rects;
        if(self.data_minus.length === 0)
        {
            num_neg_rects = 0;
        } else {
            num_neg_rects = self.data_minus[0].data.length;
        }


        //Plus and minus labels for the grids
        if(num_pos_rects > 0){
            var pos_grid_center = Math.floor(pos_grid_start + (self.h*num_pos_rects) / 2);
            self.grid_plus_label
                .attr("y", pos_grid_center+9) // Y lines up with baseline, need to offset to vertically center
                .text("+");                   // Some browsers have an attribute to do this, doesn't work on FireFox
        }
        else
        {
            self.grid_plus_label.text("");
        }

        if(num_neg_rects > 0){
            var neg_grid_center = Math.floor(neg_grid_start + (self.h*num_neg_rects) / 2);
            self.grid_minus_label
                .attr("y", neg_grid_center+9) // See note for plus label above
                .text("-");
        }
        else
        {
            self.grid_minus_label.text("");
        }
    };
};
