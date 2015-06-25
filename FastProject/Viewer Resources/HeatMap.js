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
        .range(["blue", "yellow", "red"]);

    this.data = [];

}

HeatMap.prototype.setData = function(data)
{
    this.data = data;
    this.redraw(true)();
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
            .attr('x', function(d,i) {
                return (i * self.w);
            });

        heatmapRects.style('fill',function(d) {
            return self.colorScale(d);
        });

        heatmapRects.exit().remove();

    }
};
