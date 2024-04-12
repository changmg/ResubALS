module deviation(a, b, f);
parameter width = 32;
input [width - 1: 0] a;
input [width - 1: 0] b;
output [5: 0] f;
wire [width - 1: 0] diff;
assign diff = a ^ b;
assign f = 1'b0 + diff[0] + diff[1] + diff[2] + diff[3] + diff[4] + diff[5] + diff[6] + diff[7] + diff[8] + diff[9] + diff[10] + diff[11] + diff[12] + diff[13] + diff[14] + diff[15] + diff[16] + diff[17] + diff[18] + diff[19] + diff[20] + diff[21] + diff[22] + diff[23] + diff[24] + diff[25] + diff[26] + diff[27] + diff[28] + diff[29] + diff[30] + diff[31];
endmodule
