/*
   Initializes a zoomable scatter plot in the element "parent"
   parent = , for example, "#chart_div"
*/
function ColorScatter(parent)
{
    self = this;
    var xdomain = [-2, 2];
    var ydomain = [-2, 2];

    this.margin = {top: 20, right: 20, bottom: 30, left: 40};
    this.width = 500 - this.margin.left - this.margin.right;
    this.height = 375 - this.margin.top - this.margin.bottom;

    this.x = d3.scale.linear()
        .domain(xdomain)
        .range([0, self.width]);

    this.y = d3.scale.linear()
        .domain(ydomain)
        .range([self.height, 0]);

    this.xAxis = d3.svg.axis()
        .scale(self.x)
        .orient("bottom")
        .tickSize(-self.height);

    this.yAxis = d3.svg.axis()
        .scale(self.y)
        .orient("left")
        .ticks(5)
        .tickSize(-self.width);

    this.zoom = d3.behavior.zoom()
        .x(self.x)
        .y(self.y)
        .scaleExtent([.2, 32])
        .on("zoom", self.redraw());

    //This gets overwritten when you set data
    this.colorScale = d3.scale.linear()
        .domain([0,.5,1])
        .range(["blue", "green", "red"]);

    this.svg = d3.select(parent).append("svg")
        .attr("width", self.width + self.margin.left + self.margin.right)
        .attr("height", self.height + self.margin.top + self.margin.bottom)
        .append("g")
        .attr("transform", "translate(" + self.margin.left + "," + self.margin.top + ")")
        .call(self.zoom);

    this.svg.append("rect")
        .attr("width", self.width)
        .attr("height", self.height);

    this.svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + self.height + ")")
        .call(self.xAxis);

    this.svg.append("g")
        .attr("class", "y axis")
        .call(self.yAxis);

    this.points = [];

    this.selected = -1;
    this.selected_links = [];

    this.hovered = -1;
    this.hovered_links = [];

    this.last_event = -1;
}

ColorScatter.prototype.setData = function(points)
{
    this.points = points;

    cvals = points.map(function(e){return e[2];}); //extract 3rd column
    low = d3.min(cvals);
    high = d3.max(cvals);
    mid = (low+high)/2;

    this.colorScale = d3.scale.linear()
        .domain([low, mid, high])
        .range(["blue", "green", "red"]);

    this.redraw(true)();
};

ColorScatter.prototype.setSelected = function(selected_index, event_id)
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

ColorScatter.prototype.setHovered = function(hovered_index, event_id)
{
    if(event_id === undefined){
        event_id = Math.random();
    }

    //Needed to prevent infinite loops with linked hover and select events
    if(this.last_event != event_id) {
        this.last_event = event_id;
        this.hover_col = hovered_index;
        this.svg.selectAll("circle")
            .classed("point-hover", function (d, i) {
                return i == hovered_index;
            });

        this.hovered_links.forEach(function (e, i) {
            e.setHovered(hovered_index, event_id);
        });
    }
};

ColorScatter.prototype.redraw = function(performTransition) {
    var self = this;
    return function(){
    self.svg.select(".x.axis").call(self.xAxis);
    self.svg.select(".y.axis").call(self.yAxis);
    var circles = self.svg.selectAll("circle")
        .data(self.points);

    circles.enter().append("circle").attr("r",4.5);
    circles.style("fill", function(d){return self.colorScale(d[2]);})
        .on("click", function(d,i){self.setSelected(i)})
        .on("mouseover", function(d,i){self.setHovered(i);})
        .on("mouseout", function(d,i){self.setHovered(-1);})
        .classed("point-selected", function(d, i){return i == self.selected});

    if(performTransition !== undefined && performTransition == true)
    {
        circles
            .transition()
            .duration(1000)
            .attr("cx", function(d){return self.x(d[0]);})
            .attr("cy", function(d){return self.y(d[1]);});
    }
    else
    {
        circles
            .attr("cx", function(d){return self.x(d[0]);})
            .attr("cy", function(d){return self.y(d[1]);})
    }

    circles.exit().remove();
    }
};