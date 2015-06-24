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
	this.width = 800 - this.margin.left - this.margin.right;
	this.height = 600 - this.margin.top - this.margin.bottom;

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
	
	this.colorScale = d3.scale.linear()
	.domain([0, 60, 130])
	.range(["red", "green", "blue"]);
	
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
}

ColorScatter.prototype.setData = function(points)
{
	this.points = points;
	this.redraw(true)();
}

ColorScatter.prototype.redraw = function(performTransition) {
	  var self = this;
	  return function(){
		self.svg.select(".x.axis").call(self.xAxis);
		self.svg.select(".y.axis").call(self.yAxis);
		var circles = self.svg.selectAll("circle")
			.data(self.points);

		circles.enter().append("circle").attr("r",4.5);
		if(performTransition !== undefined && performTransition == true)
		{
			circles
				.style("fill", function(d){return self.colorScale(d[2]);})
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
				.style("fill", function(d){return self.colorScale(d[2]);});
		}
		
		circles.exit().remove();
	  }
	}