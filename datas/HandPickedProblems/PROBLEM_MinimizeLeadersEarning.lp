\ENCODING=ISO-8859-1
\Problem name: IloCplex

Minimize
 obj1: x_0_9 - x_1_4 - x_2_5 + x_0_6 - x_0_5 - x_0_4 + x_1_6 - x_1_5 + x_1_9
       - x_2_4 + x_2_6 + x_2_9
Subject To
 c1:   u_00 <= 2
 c2:   u_01 <= 0
 c3:   u_02 <= 3
 c4:   u_03 <= 1
 c5:   u_04 <= 2
 c6:   u_05 <= 0
 c7:   u_06 <= 2
 c8:   u_07 <= 3
 c9:   u_08 <= 0
 c10:  u_09 <= 2
 c11:  u_010 <= 3
 c12:  u_10 <= 0
 c13:  u_11 <= 2
 c14:  u_12 <= 0
 c15:  u_13 <= 3
 c16:  u_14 <= 3
 c17:  u_15 <= 0
 c18:  u_16 <= 3
 c19:  u_17 <= 2
 c20:  u_18 <= 1
 c21:  u_19 <= 2
 c22:  u_110 <= 1
 c23:  u_20 <= 0
 c24:  u_21 <= 0
 c25:  u_22 <= 2
 c26:  u_23 <= 3
 c27:  u_24 <= 2
 c28:  u_25 <= 0
 c29:  u_26 <= 3
 c30:  u_27 <= 1
 c31:  u_28 <= 1
 c32:  u_29 <= 1
 c33:  u_210 <= 3
 c34:  u_02 + u_04 <= 6.27251
 c35:  u_05 + u_07 <= 4.07877
 c36:  u_05 + u_010 <= 4.19895
 c37:  u_04 + u_05 <= 5.18422
 c38:  u_06 + u_08 <= 5.54063
 c39:  u_03 + u_07 <= 3.84333
 c40:  u_01 + u_04 <= 5.26472
 c41:  u_02 + u_05 + u_06 + u_08 <= 4.38171
 c42:  u_09 + u_010 <= 6.23465
 c43:  u_04 + u_010 <= 4.60373
 c44:  u_00 + u_05 + u_06 <= 4.84786
 c45:  u_00 + u_02 <= 3.28148
 c46:  u_01 + u_04 + u_07 <= 5.44064
 c47:  u_02 + u_06 + u_07 + u_09 <= 5.47817
 c48:  u_03 + u_07 + u_08 + u_010 <= 2.27754
 c49:  u_00 + u_05 <= 3.58434
 c50:  u_03 + u_05 <= 3.76305
 c51:  u_05 + u_08 + u_010 <= 6.16789
 c52:  u_06 + u_07 <= 5.28225
 c53:  u_01 + u_08 <= 5.98274
 c54:  u_00 + u_01 <= 5.7473
 c55:  u_02 + u_07 <= 5.30214
 c56:  u_08 + u_010 <= 5.89941
 c57:  u_02 + u_03 + u_04 + u_09 <= 5.2585
 c58:  u_04 + u_06 <= 5.3526
 c59:  u_06 + u_08 <= 6.01278
 c60:  u_00 + u_02 + u_07 + u_08 + u_010 <= 5.54811
 c61:  u_04 + u_09 <= 6.23738
 c62:  u_00 + u_01 + u_06 <= 4.40089
 c63:  u_00 + u_04 <= 4.14016
 c64:  x_0_9 >= 1
 c65:  - 1000 hap_1_0 + ap_0 <= 0
 c66:  x_0_9 + 2 hap_1_0 <= 3
 c67:  - x_0_3 - x_0_2 >= -1
 c68:  - 1000 han_2_0 + an_0 <= 0
 c69:  - x_0_3 - x_0_2 + 2 han_2_0 <= 1
 c70:  x_1_4 + x_1_3 + x_1_0 >= 1
 c71:  - 1000 hap_1_1 + ap_1 <= 0
 c72:  x_1_4 + x_1_3 + x_1_0 + 2 hap_1_1 <= 3
 c73:  - x_1_10 >= -1
 c74:  - 1000 han_2_1 + an_1 <= 0
 c75:  - x_1_10 + 2 han_2_1 <= 1
 c76:  x_2_7 + x_2_5 + x_2_1 >= 1
 c77:  - 1000 hap_1_2 + ap_2 <= 0
 c78:  x_2_7 + x_2_5 + x_2_1 + 2 hap_1_2 <= 3
 c79:  - x_2_1 - x_2_0 >= -1
 c80:  - 1000 han_2_2 + an_2 <= 0
 c81:  - x_2_1 - x_2_0 + 2 han_2_2 <= 1
 c82:  x_0_6 - x_0_1 - x_0_0 >= 0
 c83:  - 1000 hb_00 + b_0_0 <= 0
 c84:  x_0_6 - x_0_1 - x_0_0 + 2 hb_00 <= 2
 c85:  x_0_2 - x_0_5 - x_0_4 >= 0
 c86:  - 1000 hb_02 + b_0_2 <= 0
 c87:  x_0_2 - x_0_5 - x_0_4 + 2 hb_02 <= 2
 c88:  x_0_3 - x_0_6 + x_0_0 + x_0_4 - x_0_7 >= 0
 c89:  - 1000 hb_03 + b_0_3 <= 0
 c90:  x_0_3 - x_0_6 + x_0_0 + x_0_4 - x_0_7 + 2 hb_03 <= 2
 c91:  - x_0_9 + x_0_10 - x_0_8 >= 0
 c92:  - 1000 hb_04 + b_0_4 <= 0
 c93:  - x_0_9 + x_0_10 - x_0_8 + 2 hb_04 <= 2
 c94:  x_0_1 + x_0_5 + x_0_7 - x_0_10 >= 0
 c95:  - 1000 hb_06 + b_0_6 <= 0
 c96:  x_0_1 + x_0_5 + x_0_7 - x_0_10 + 2 hb_06 <= 2
 c97:  - x_1_0 + x_1_6 - x_1_1 >= 0
 c98:  - 1000 hb_10 + b_1_0 <= 0
 c99:  - x_1_0 + x_1_6 - x_1_1 + 2 hb_10 <= 2
 c100: - x_1_3 + x_1_8 - x_1_2 >= 0
 c101: - 1000 hb_11 + b_1_1 <= 0
 c102: - x_1_3 + x_1_8 - x_1_2 + 2 hb_11 <= 2
 c103: - x_1_4 + x_1_2 - x_1_5 >= 0
 c104: - 1000 hb_12 + b_1_2 <= 0
 c105: - x_1_4 + x_1_2 - x_1_5 + 2 hb_12 <= 2
 c106: x_1_10 - x_1_8 - x_1_9 >= 0
 c107: - 1000 hb_14 + b_1_4 <= 0
 c108: x_1_10 - x_1_8 - x_1_9 + 2 hb_14 <= 2
 c109: x_1_9 >= 0
 c110: - 1000 hb_15 + b_1_5 <= 0
 c111: x_1_9 + 2 hb_15 <= 2
 c112: x_2_8 - x_2_3 - x_2_2 >= 0
 c113: - 1000 hb_21 + b_2_1 <= 0
 c114: x_2_8 - x_2_3 - x_2_2 + 2 hb_21 <= 2
 c115: - x_2_5 + x_2_2 - x_2_4 >= 0
 c116: - 1000 hb_22 + b_2_2 <= 0
 c117: - x_2_5 + x_2_2 - x_2_4 + 2 hb_22 <= 2
 c118: - x_2_7 + x_2_0 + x_2_3 + x_2_4 - x_2_6 >= 0
 c119: - 1000 hb_23 + b_2_3 <= 0
 c120: - x_2_7 + x_2_0 + x_2_3 + x_2_4 - x_2_6 + 2 hb_23 <= 2
 c121: - x_2_8 + x_2_10 - x_2_9 >= 0
 c122: - 1000 hb_24 + b_2_4 <= 0
 c123: - x_2_8 + x_2_10 - x_2_9 + 2 hb_24 <= 2
 c124: x_2_9 >= 0
 c125: - 1000 hb_25 + b_2_5 <= 0
 c126: x_2_9 + 2 hb_25 <= 2
 c127: - u_00 - b_0_0 + b_0_3 <= 0
 c128: x_0_0 - 2 hdu_0_0 <= 0
 c129: u_00 + b_0_0 - b_0_3 + 1000 hdu_0_0 <= 1000
 c130: - u_01 - b_0_0 + b_0_6 <= 1
 c131: x_0_1 - 2 hdu_0_1 <= 0
 c132: u_01 + b_0_0 - b_0_6 + 1000 hdu_0_1 <= 999
 c133: - u_02 - an_0 + b_0_2 <= 0
 c134: x_0_2 - 2 hdu_0_2 <= 0
 c135: u_02 + an_0 - b_0_2 + 1000 hdu_0_2 <= 1000
 c136: - u_03 - an_0 + b_0_3 <= 0
 c137: x_0_3 - 2 hdu_0_3 <= 0
 c138: u_03 + an_0 - b_0_3 + 1000 hdu_0_3 <= 1000
 c139: - u_04 - b_0_2 + b_0_3 <= 0
 c140: x_0_4 - 2 hdu_0_4 <= 0
 c141: u_04 + b_0_2 - b_0_3 + 1000 hdu_0_4 <= 1000
 c142: - u_05 - b_0_2 + b_0_6 <= 0
 c143: x_0_5 - 2 hdu_0_5 <= 0
 c144: u_05 + b_0_2 - b_0_6 + 1000 hdu_0_5 <= 1000
 c145: - u_06 + b_0_0 - b_0_3 <= 1
 c146: x_0_6 - 2 hdu_0_6 <= 0
 c147: u_06 - b_0_0 + b_0_3 + 1000 hdu_0_6 <= 999
 c148: - u_07 - b_0_3 + b_0_6 <= 0
 c149: x_0_7 - 2 hdu_0_7 <= 0
 c150: u_07 + b_0_3 - b_0_6 + 1000 hdu_0_7 <= 1000
 c151: - u_08 - b_0_4 <= 0
 c152: x_0_8 - 2 hdu_0_8 <= 0
 c153: u_08 + b_0_4 + 1000 hdu_0_8 <= 1000
 c154: - u_09 + ap_0 - b_0_4 <= 2
 c155: x_0_9 - 2 hdu_0_9 <= 0
 c156: u_09 - ap_0 + b_0_4 + 1000 hdu_0_9 <= 998
 c157: - u_010 + b_0_4 - b_0_6 <= 1
 c158: x_0_10 - 2 hdu_0_10 <= 0
 c159: u_010 - b_0_4 + b_0_6 + 1000 hdu_0_10 <= 999
 c160: - u_10 + ap_1 - b_1_0 <= 0
 c161: x_1_0 - 2 hdu_1_0 <= 0
 c162: u_10 - ap_1 + b_1_0 + 1000 hdu_1_0 <= 1000
 c163: - u_11 - b_1_0 <= 1
 c164: x_1_1 - 2 hdu_1_1 <= 0
 c165: u_11 + b_1_0 + 1000 hdu_1_1 <= 999
 c166: - u_12 - b_1_1 + b_1_2 <= 0
 c167: x_1_2 - 2 hdu_1_2 <= 0
 c168: u_12 + b_1_1 - b_1_2 + 1000 hdu_1_2 <= 1000
 c169: - u_13 + ap_1 - b_1_1 <= 0
 c170: x_1_3 - 2 hdu_1_3 <= 0
 c171: u_13 - ap_1 + b_1_1 + 1000 hdu_1_3 <= 1000
 c172: - u_14 + ap_1 - b_1_2 <= 0
 c173: x_1_4 - 2 hdu_1_4 <= 0
 c174: u_14 - ap_1 + b_1_2 + 1000 hdu_1_4 <= 1000
 c175: - u_15 - b_1_2 <= 0
 c176: x_1_5 - 2 hdu_1_5 <= 0
 c177: u_15 + b_1_2 + 1000 hdu_1_5 <= 1000
 c178: - u_16 + b_1_0 <= 1
 c179: x_1_6 - 2 hdu_1_6 <= 0
 c180: u_16 - b_1_0 + 1000 hdu_1_6 <= 999
 c181: - u_17 <= 0
 c182: - 2 hdu_1_7 + x_1_7 <= 0
 c183: u_17 + 1000 hdu_1_7 <= 1000
 c184: - u_18 + b_1_1 - b_1_4 <= 0
 c185: x_1_8 - 2 hdu_1_8 <= 0
 c186: u_18 - b_1_1 + b_1_4 + 1000 hdu_1_8 <= 1000
 c187: - u_19 - b_1_4 + b_1_5 <= 2
 c188: x_1_9 - 2 hdu_1_9 <= 0
 c189: u_19 + b_1_4 - b_1_5 + 1000 hdu_1_9 <= 998
 c190: - u_110 - an_1 + b_1_4 <= 1
 c191: x_1_10 - 2 hdu_1_10 <= 0
 c192: u_110 + an_1 - b_1_4 + 1000 hdu_1_10 <= 999
 c193: - u_20 - an_2 + b_2_3 <= 0
 c194: x_2_0 - 2 hdu_2_0 <= 0
 c195: u_20 + an_2 - b_2_3 + 1000 hdu_2_0 <= 1000
 c196: - u_21 + ap_2 - an_2 <= 1
 c197: x_2_1 - 2 hdu_2_1 <= 0
 c198: u_21 - ap_2 + an_2 + 1000 hdu_2_1 <= 999
 c199: - u_22 - b_2_1 + b_2_2 <= 0
 c200: x_2_2 - 2 hdu_2_2 <= 0
 c201: u_22 + b_2_1 - b_2_2 + 1000 hdu_2_2 <= 1000
 c202: - u_23 - b_2_1 + b_2_3 <= 0
 c203: x_2_3 - 2 hdu_2_3 <= 0
 c204: u_23 + b_2_1 - b_2_3 + 1000 hdu_2_3 <= 1000
 c205: - u_24 - b_2_2 + b_2_3 <= 0
 c206: x_2_4 - 2 hdu_2_4 <= 0
 c207: u_24 + b_2_2 - b_2_3 + 1000 hdu_2_4 <= 1000
 c208: - u_25 + ap_2 - b_2_2 <= 0
 c209: x_2_5 - 2 hdu_2_5 <= 0
 c210: u_25 - ap_2 + b_2_2 + 1000 hdu_2_5 <= 1000
 c211: - u_26 - b_2_3 <= 1
 c212: x_2_6 - 2 hdu_2_6 <= 0
 c213: u_26 + b_2_3 + 1000 hdu_2_6 <= 999
 c214: - u_27 + ap_2 - b_2_3 <= 0
 c215: x_2_7 - 2 hdu_2_7 <= 0
 c216: u_27 - ap_2 + b_2_3 + 1000 hdu_2_7 <= 1000
 c217: - u_28 + b_2_1 - b_2_4 <= 0
 c218: x_2_8 - 2 hdu_2_8 <= 0
 c219: u_28 - b_2_1 + b_2_4 + 1000 hdu_2_8 <= 1000
 c220: - u_29 - b_2_4 + b_2_5 <= 2
 c221: x_2_9 - 2 hdu_2_9 <= 0
 c222: u_29 + b_2_4 - b_2_5 + 1000 hdu_2_9 <= 998
 c223: - u_210 + b_2_4 <= 1
 c224: x_2_10 - 2 hdu_2_10 <= 0
 c225: u_210 - b_2_4 + 1000 hdu_2_10 <= 999
Bounds
 0 <= x_0_9 <= 1
 0 <= hap_1_0 <= 1
 0 <= x_0_3 <= 1
 0 <= x_0_2 <= 1
 0 <= han_2_0 <= 1
 0 <= x_1_4 <= 1
 0 <= x_1_3 <= 1
 0 <= x_1_0 <= 1
 0 <= hap_1_1 <= 1
 0 <= x_1_10 <= 1
 0 <= han_2_1 <= 1
 0 <= x_2_7 <= 1
 0 <= x_2_5 <= 1
 0 <= x_2_1 <= 1
 0 <= hap_1_2 <= 1
 0 <= x_2_0 <= 1
 0 <= han_2_2 <= 1
 0 <= x_0_6 <= 1
 0 <= x_0_1 <= 1
 0 <= x_0_0 <= 1
 0 <= hb_00 <= 1
 0 <= x_0_5 <= 1
 0 <= x_0_4 <= 1
 0 <= hb_02 <= 1
 0 <= x_0_7 <= 1
 0 <= hb_03 <= 1
 0 <= x_0_10 <= 1
 0 <= x_0_8 <= 1
 0 <= hb_04 <= 1
 0 <= hb_06 <= 1
 0 <= x_1_6 <= 1
 0 <= x_1_1 <= 1
 0 <= hb_10 <= 1
 0 <= x_1_8 <= 1
 0 <= x_1_2 <= 1
 0 <= hb_11 <= 1
 0 <= x_1_5 <= 1
 0 <= hb_12 <= 1
 0 <= x_1_9 <= 1
 0 <= hb_14 <= 1
 0 <= hb_15 <= 1
 0 <= x_2_8 <= 1
 0 <= x_2_3 <= 1
 0 <= x_2_2 <= 1
 0 <= hb_21 <= 1
 0 <= x_2_4 <= 1
 0 <= hb_22 <= 1
 0 <= x_2_6 <= 1
 0 <= hb_23 <= 1
 0 <= x_2_10 <= 1
 0 <= x_2_9 <= 1
 0 <= hb_24 <= 1
 0 <= hb_25 <= 1
 0 <= hdu_0_0 <= 1
 0 <= hdu_0_1 <= 1
 0 <= hdu_0_2 <= 1
 0 <= hdu_0_3 <= 1
 0 <= hdu_0_4 <= 1
 0 <= hdu_0_5 <= 1
 0 <= hdu_0_6 <= 1
 0 <= hdu_0_7 <= 1
 0 <= hdu_0_8 <= 1
 0 <= hdu_0_9 <= 1
 0 <= hdu_0_10 <= 1
 0 <= hdu_1_0 <= 1
 0 <= hdu_1_1 <= 1
 0 <= hdu_1_2 <= 1
 0 <= hdu_1_3 <= 1
 0 <= hdu_1_4 <= 1
 0 <= hdu_1_5 <= 1
 0 <= hdu_1_6 <= 1
 0 <= hdu_1_7 <= 1
 0 <= x_1_7 <= 1
 0 <= hdu_1_8 <= 1
 0 <= hdu_1_9 <= 1
 0 <= hdu_1_10 <= 1
 0 <= hdu_2_0 <= 1
 0 <= hdu_2_1 <= 1
 0 <= hdu_2_2 <= 1
 0 <= hdu_2_3 <= 1
 0 <= hdu_2_4 <= 1
 0 <= hdu_2_5 <= 1
 0 <= hdu_2_6 <= 1
 0 <= hdu_2_7 <= 1
 0 <= hdu_2_8 <= 1
 0 <= hdu_2_9 <= 1
 0 <= hdu_2_10 <= 1
Binaries
 hap_1_0  han_2_0  hap_1_1  han_2_1  hap_1_2  han_2_2  hb_00  hb_02  hb_03 
 hb_04  hb_06  hb_10  hb_11  hb_12  hb_14  hb_15  hb_21  hb_22  hb_23  hb_24 
 hb_25  hdu_0_0  hdu_0_1  hdu_0_2  hdu_0_3  hdu_0_4  hdu_0_5  hdu_0_6  hdu_0_7 
 hdu_0_8  hdu_0_9  hdu_0_10  hdu_1_0  hdu_1_1  hdu_1_2  hdu_1_3  hdu_1_4 
 hdu_1_5  hdu_1_6  hdu_1_7  hdu_1_8  hdu_1_9  hdu_1_10  hdu_2_0  hdu_2_1 
 hdu_2_2  hdu_2_3  hdu_2_4  hdu_2_5  hdu_2_6  hdu_2_7  hdu_2_8  hdu_2_9 
 hdu_2_10 
End
