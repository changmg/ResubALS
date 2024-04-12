module deviation(a, b, f);
parameter width = 9;
input [width - 1: 0] a;
input [width - 1: 0] b;
output [3: 0] f;
wire [width - 1: 0] diff;
assign diff = a ^ b;
assign f = 1'b0 + diff[0] + diff[1] + diff[2] + diff[3] + diff[4] + diff[5] + diff[6] + diff[7] + diff[8];
endmodule
