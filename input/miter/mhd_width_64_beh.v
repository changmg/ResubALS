module deviation(a, b, f);
parameter width = 64;
input [width - 1: 0] a;
input [width - 1: 0] b;
output [6: 0] f;
wire [width - 1: 0] diff;
assign diff = a ^ b;
assign f = 1'b0 + diff[0] + diff[1] + diff[2] + diff[3] + diff[4] + diff[5] + diff[6] + diff[7] + diff[8] + diff[9] + diff[10] + diff[11] + diff[12] + diff[13] + diff[14] + diff[15] + diff[16] + diff[17] + diff[18] + diff[19] + diff[20] + diff[21] + diff[22] + diff[23] + diff[24] + diff[25] + diff[26] + diff[27] + diff[28] + diff[29] + diff[30] + diff[31] + diff[32] + diff[33] + diff[34] + diff[35] + diff[36] + diff[37] + diff[38] + diff[39] + diff[40] + diff[41] + diff[42] + diff[43] + diff[44] + diff[45] + diff[46] + diff[47] + diff[48] + diff[49] + diff[50] + diff[51] + diff[52] + diff[53] + diff[54] + diff[55] + diff[56] + diff[57] + diff[58] + diff[59] + diff[60] + diff[61] + diff[62] + diff[63];
endmodule
