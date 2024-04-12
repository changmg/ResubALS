module deviation(a, b, f);
parameter width = 25;
input [width - 1: 0] a;
input [width - 1: 0] b;
output [width * 2 - 1: 0] f;
wire [width - 1: 0] diff;
assign diff = (a > b)? (a - b): (b - a);
assign f = diff * diff;
endmodule
