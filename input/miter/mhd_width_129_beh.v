module deviation(a, b, f);
parameter width = 129;
input [width - 1: 0] a;
input [width - 1: 0] b;
output [7: 0] f;
wire [width - 1: 0] diff;
assign diff = a ^ b;
assign f = 1'b0 + diff[0] + diff[1] + diff[2] + diff[3] + diff[4] + diff[5] + diff[6] + diff[7] + diff[8] + diff[9] + diff[10] + diff[11] + diff[12] + diff[13] + diff[14] + diff[15] + diff[16] + diff[17] + diff[18] + diff[19] + diff[20] + diff[21] + diff[22] + diff[23] + diff[24] + diff[25] + diff[26] + diff[27] + diff[28] + diff[29] + diff[30] + diff[31] + diff[32] + diff[33] + diff[34] + diff[35] + diff[36] + diff[37] + diff[38] + diff[39] + diff[40] + diff[41] + diff[42] + diff[43] + diff[44] + diff[45] + diff[46] + diff[47] + diff[48] + diff[49] + diff[50] + diff[51] + diff[52] + diff[53] + diff[54] + diff[55] + diff[56] + diff[57] + diff[58] + diff[59] + diff[60] + diff[61] + diff[62] + diff[63] + diff[64] + diff[65] + diff[66] + diff[67] + diff[68] + diff[69] + diff[70] + diff[71] + diff[72] + diff[73] + diff[74] + diff[75] + diff[76] + diff[77] + diff[78] + diff[79] + diff[80] + diff[81] + diff[82] + diff[83] + diff[84] + diff[85] + diff[86] + diff[87] + diff[88] + diff[89] + diff[90] + diff[91] + diff[92] + diff[93] + diff[94] + diff[95] + diff[96] + diff[97] + diff[98] + diff[99] + diff[100] + diff[101] + diff[102] + diff[103] + diff[104] + diff[105] + diff[106] + diff[107] + diff[108] + diff[109] + diff[110] + diff[111] + diff[112] + diff[113] + diff[114] + diff[115] + diff[116] + diff[117] + diff[118] + diff[119] + diff[120] + diff[121] + diff[122] + diff[123] + diff[124] + diff[125] + diff[126] + diff[127] + diff[128];
endmodule
