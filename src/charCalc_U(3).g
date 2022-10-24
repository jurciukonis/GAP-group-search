# In this program we look for characters of the groups listed below, we find the highest powers (group centers), we identify classes, 
# and we find faithful irreducible representations (irrep) of the groups.
# The group series counted here are classified in ref. [2].
# For more information about centers and computation methods see ref. [1]. 
# For more information on each discrete subgroups of U(3) and SU(3) see ref. [2].

# Author: Darius Jurciukonis
# Date: 2022-07

# If you use this code, please cite these references:
# [1] D. Jurciukonis and L. Lavoura, The centers of discrete groups as stabilizers of Dark Matter, arXiv:2210.12133.
# [2] D. Jurciukonis and L. Lavoura, GAP listing of the finite subgroups of U(3) of order smaller than 2000,  PTEP 2017 (2017) 5, 053A03, arXiv:1702.00005.

#------------------------------------------------------------------
# this program can be run in GAP with the command
# Read("~/GAP/CharCalc/charCalc_U(3).g");
#------------------------------------------------------------------

# List that contains series of groups that have faithful irreducible representations of dimension 3 with matrices with determinant different from 1.
nbGlistU3:=[
# 1) T(m)
[[63, 1], [189, 1], [567, 1], [1701, 68], [117, 1], [351, 1], [1053, 16], [171, 1], [513, 1], [1539, 16], [279, 1], [837, 1], [333, 1], [999, 1], [387, 1], [1161, 6], [441, 1], [1323, 1], [549, 1], [1647, 6], [603, 1], [1809, 6], [657, 1], [1971, 6], [711, 1], [819, 4], [819, 3], [873, 1], [927, 1], [981, 1], [1143, 1], [1197, 3], [1197, 4], [1251, 1], [1359, 1], [1413, 1], [1467, 1], [1521, 1], [1629, 1], [1737, 1], [1791, 1], [1899, 1], [1953, 3], [1953, 4]],
# 2) Delta(3n,m)
[[36, 3], [108, 3], [324, 3], [972, 3], [144, 3], [432, 3], [1296, 3], [225, 3], [675, 5], [441, 7], [1323, 14], [576, 3], [1728, 3], [900, 66], [1089, 3], [1521, 7], [1764, 91]],
# 3) S4(j)
[[48, 30], [96, 65], [192,186], [384,581], [768, 1085351], [1536, 408544687]],
# 4) Delta(6n,j)
[[108, 11], [216, 17], [432, 33], [864, 69], [1728, 185], [192, 182], [384, 571], [768, 1085333], [1536, 408544678], [300, 13], [600, 45], [1200, 183], [432, 260], [864, 703], [1728, 2855], [588, 16], [1176, 57], [768, 1085335], [1536, 408544641], [972, 64], [1944, 70], [1200, 682], [1452, 11], [1728, 2847]],
# 5) DeltaPrime(6n,m,j)
[[162, 44], [324, 102], [648, 244], [1296, 647], [486, 164], [972, 348], [1944, 746], [1458, 1354], [648, 563], [1296, 2113], [1944,   2415], [1458, 1371]],
# 6) L(n,m)
[[252, 11], [468, 14], [684, 11], [756, 11], [1008, 57], [1116, 11], [1332, 14], [1404, 14], [1548, 11], [1575, 7], [1764,11], [1872, 60]],
# 7) P(m)
[[189, 7], [351, 7], [513, 8], [837, 7], [999, 8], [1161, 12], [1323, 7], [1647, 12], [1809, 12], [1971, 12], [567, 7], [1053, 27], [1539, 27], [1701, 128]],
# 8) Q(m)
[[189, 4], [351, 4], [513, 5], [567, 4], [837, 5], [999, 6], [1053, 26], [1161, 10], [1323, 4], [1539, 26], [1647, 10], [1701, 127], [1809, 10], [1971, 11]],
# 9) Qprime(m)
[[189, 5], [351, 5], [513, 6], [567, 5], [837, 4], [999, 5], [1053, 25], [1161, 11], [1323, 5], [1539, 25], [1647, 11], [1701, 126], [1809, 11], [1971, 10]],
# 10) X(n)
[[27, 4], [108, 21], [243, 27], [432, 102], [675, 11], [972, 123], [1323, 42], [1728, 1290]],
# 11) S(m)
[[567, 36], [1053, 47], [1539, 47], [1701, 240]],
# 12) Sprime(m)
[[567, 12], [1053, 32], [1539, 32], [1701, 115]],
# 13) Y(m)
[[567, 23], [1053, 29], [1539, 29], [1701, 261]],
# 14) V(m)
[[567, 14], [1053, 37], [1539, 37], [1701, 138]],
# 15) M
[[756, 113], [1404, 137]],
# 16) Mprime
[[756, 114], [1404, 138]],
# 17) J
[[756, 116], [1404, 140]],
# 18) W(n,m)
[[27, 4], [81, 6], [243, 24], [729, 94], [108, 19], [324, 43], [972, 117], [432, 100], [1296, 220], [675, 9], [1323, 40], [1728, 1286]],
# 19) Z(n,m)
[[81, 14], [324, 128], [729, 397], [1296, 1499], [243, 50], [972, 520], [729, 393]],
# 20) Zprime(n,m)
[[81, 8], [324, 49], [729, 397], [1296, 227], [243, 20], [972, 152], [729, 64]],
# 21) ZprimePrime(n,m)
[[81, 9], [324, 50], [243, 26], [1296, 228], [243, 19], [972, 153], [729, 63]],
# 22) Z(n,m,j)
[[162, 12], [324, 15], [648, 21], [1296, 37], [486, 28], [972, 31], [1944, 37], [1458, 618], [648, 260], [1296, 689], [1944, 833]],
# 23) Zprime(n,m,j)
[[162, 14], [324, 17], [648, 23], [1296, 39], [486, 26], [972, 29], [1944, 35], [1458, 615], [648, 259], [1296, 688], [1944, 832]],
# 24) H(n,m,j)
[[486, 125], [972, 309], [1944, 707], [1458, 1095], [1944, 2363]],
# 25) G(m,j)
[[324, 13], [972, 309], [648, 19], [1944, 707], [1296, 35], [1296, 699]],
# 26) Y(m,j)
[[324, 45], [972, 147], [1296, 222]],
# 27) Y(j)
[[81, 7], [324, 60], [1296, 237]],
# 28) Ytilde(j)
[[162, 10], [648, 266]],
# 29) U(n,m,j)
[[243, 55], [729, 86], [729, 284], [972, 550]],
# 30) L(m)
[[243, 16], [729, 62], [1701, 102]],
# 31) V(j)
[[81, 10], [324, 51], [1296, 226]],
# 32) D(j)
[[243, 25], [972, 121]],
# 33) J(m)
[[243, 27], [729, 80]],
# 34) 729
[[729, 96], [729, 97], [729, 98]],
# 35) 972
[[972, 170]],
# 36) 1458
[[1458, 663], [1458, 666]],
# 37) 1701
[[1701, 112], [1701, 130], [1701, 131]],
# 38) Xi(m,j)
[[108, 15], [216, 25], [432, 57], [864, 194], [1728, 953], [324, 111], [648, 352], [1296, 1239], [972, 411], [1944, 1123]],
# 39) Xihat(m,j)
[[432, 273], [864, 737], [1728, 2929], [1296, 2203]],
# 40) Pi(m,j)
[[432, 239], [864, 675], [1728, 2785], [1296, 1995]],
# 41) Theta(m)
[[216, 88], [648, 551], [1944, 2333]],
# 42) Upsilon(m)
[[648, 531], [1944, 2293]],
# 43) UpsilonPrime(m)
[[648, 532], [1944, 2294]],
# 44) Omega(m)
[[648, 533], [1944, 3448]]
];

# List that contains series of subgroups of SU(3).
nbGlistSU3:=[
# 1) Delta(3n)
[[ 12, 3 ],[ 27, 3 ],[ 48, 3 ],[ 75, 2 ],[ 108, 22 ],[ 147, 5 ],[ 192, 3 ],[ 243, 26 ],[ 300, 43 ],[ 363, 2 ],[ 432, 103 ],[ 507, 5 ],[ 588, 60 ],[ 675, 12 ],[ 768, 1083477 ],[ 867, 2 ],[ 972, 122 ],[ 1083, 5 ],[ 1200, 384 ],[ 1323, 43 ],[ 1452, 34 ],[ 1587, 2 ],[ 1728, 1291 ],[ 1875, 16 ]],
# 2) Delta(6n)
[[ 24, 12 ],[ 54, 8 ],[ 96, 64 ],[ 150, 5 ],[ 216, 95 ],[ 294, 7 ],[ 384, 568 ],[ 486, 61 ],[ 600, 179 ],[ 726, 5 ],[ 864, 701 ],[ 1014, 7 ],[ 1176, 243 ],[ 1350, 46 ],[ 1536, 408544632 ],[ 1734, 5 ],[ 1944, 849 ]],
# 3) C(r,k,l)
[[21, 1], [39, 1], [57, 1], [81, 9], [84, 11], [93, 1], [111, 1], [129, 1], [147, 1], [156, 14], [183, 1], [189, 8], [201, 1], [219, 1], [228, 11], [237, 1], [273, 3], [273, 4], [291, 1], [309, 1], [324, 50], [327, 1], [336, 57], [351, 8], [372, 11], [381, 1], [399, 3], [399, 4], [417, 1], [444, 14], [453, 1], [471, 1], [489, 1], [507, 1], [513, 9], [516, 11], [525, 5], [543, 1], [567, 13], [579, 1], [588, 11], [597, 1], [624, 60], [633, 1], [651, 3], [651, 4], [669, 1], [687, 1], [723, 1], [729, 95], [732, 14], [741, 3], [741, 4], [756, 117], [777, 3], [777, 4], [804, 11], [813, 1], [831, 1], [837, 8], [849, 1], [876, 14], [903, 5], [903, 6], [912, 57], [921, 1], [939, 1], [948, 11], [975, 5], [993, 1], [999, 9], [1011, 1], [1029, 6], [1029, 9], [1047, 1], [1053, 35], [1083, 1], [1092, 68], [1092, 69], [1101, 1], [1119, 1], [1137, 1], [1161, 9], [1164, 14], [1191, 1], [1209, 3], [1209, 4], [1227, 1], [1236, 11], [1263, 1], [1281, 3], [1281, 4], [1296, 228], [1299, 1], [1308, 14], [1317, 1], [1323, 8], [1344, 393], [1371, 1], [1389, 1], [1404, 141], [1407, 3], [1407, 4], [1425, 5], [1443, 3], [1443, 4], [1461, 1], [1488, 57], [1497, 1], [1524, 11], [1533, 3], [1533, 4], [1539, 35], [1569, 1], [1596, 55], [1596, 56], [1623, 1], [1641, 1], [1647, 9], [1659, 3], [1659, 4], [1668, 11], [1677, 3], [1677, 4], [1701, 135], [1713, 1], [1731, 1], [1767, 3], [1767, 4], [1776, 60], [1803, 1], [1809, 9], [1812, 11], [1821, 1], [1839, 1], [1857, 1], [1884, 14], [1893, 1], [1911, 3], [1911, 4], [1911, 14], [1929, 1], [1956, 11], [1971, 9], [1983, 1]],
# 4) D(1,3l,l)
[[162, 14], [648, 259], [1458, 659]],
# 5) Sigma(m)
[[60, 5], [108, 15], [168, 42], [216, 88], [648, 532], [1080, 260]]
];

groupNamesU3:=["T(m)", "Delta(3n,m)", "S4(j)", "Delta(6n,j)", "DeltaPrime(6n,m,j)", "L(n,m)", "P(m)", "Q(m)", "Qprime(m)", "X(n)", "S(m)", "Sprime(m)", "Y(m)", "V(m)", "M", "Mprime", "J", "W(n,m)", "Z(n,m)", "Zprime(n,m)", "ZprimePrime(n,m)", "Z(n,m,j)", "Zprime(n,m,j)", "H(n,m,j)", "G(m,j)", "Y(m,j)", "Y(j)", "Ytilde(j)", "U(n,m,j)", "L(m)", "V(j)", "D(j)", "J(m)", "729", "972", "1458", "1701", "Xi(m,j)", "Xihat(m,j)", "Pi(m,j)", "Theta(m)", "Upsilon(m)", "UpsilonPrime(m)", "Omega(m)"];

groupNamesSU3:=["Delta(3n)", "Delta(6n)", "C(r,k,l)", "D(1,3l,l)", "Sigma(m)"];


SizeScreen([1000,20]);

# nbGlist:=[[48, 30], [96, 65], [192,186], [384,581], [768, 1085351], [1536, 408544687]];
# G:=SmallGroup(48,30);

# output folder
folder:="~/GAP/Data/charCalcU3/";

filename1:=JoinStringsWithSeparator([folder,"Grp_Test.dat"], "");
filename2:=JoinStringsWithSeparator([folder,"OutputLog.dat"], "");
filename3:=JoinStringsWithSeparator([folder,"Grp_Series_Summ.dat"], "");


# we choose which group list we will examine
nbGlistFull:=nbGlistSU3;
groupNames:=groupNamesSU3;
# nbGlistFull:=nbGlistU3;
# groupNames:=groupNamesU3;

ii2:=1;
# loop through group sets:
for ii2 in [1 .. 5] do 
filename1:=JoinStringsWithSeparator([folder,"Group_", String(groupNames[ii2]), ".dat"], "");
AppendTo(filename1, JoinStringsWithSeparator(["Analysis of groups ", String(groupNames[ii2])], ""), "\n");
AppendTo(filename1, "------------------------", "\n\n");

Print("\n", "Nb.:  ", ii2, ", group: ", groupNames[ii2], "\n"); 

nbGlist:=nbGlistFull[ii2];
#-----------------------------
# nbGlist:=[[60, 5]];

time0:=Runtime();
infoList:=[]; dimCalcListFull:=[]; maxPowListFull:=[];

# loop through groups
for ii1 in [1 .. Length(nbGlist)] do 

G:=SmallGroup(nbGlist[ii1]);
ordG:=Order(G);

#-----------------------------
# OutputLogTo(filename1);   
# Display(CharacterTable(G));
# PrintCharTable( charTab );
# InfoCharacterTable(charTab);
#----------------------------

charTab:=CharacterTable(G);
repG:=IrreducibleRepresentations(G);
irr:=Irr(charTab);    # irrep with characters

# let's make matrix of characters:
charMatr:=[];
for i1 in [1 .. Size(repG)] do 
    chi:= ValuesOfClassFunction( irr[i1] );    # we take only the character list in each irrep
    charMatr[i1]:=chi;

    # AppendTo(filename2, chi, "\n" );         # we write character matrix to file, if necessary
od;

# we transpose the character matrix so that the rows correspond to the classes
charMatrTr:=TransposedMat(charMatr);

# let's note the first line, with identical characters that indicate the irrep dimension
irrepDimList := charMatrTr[1]; 

# we count the dimensions of irreps
dimCalcListFull[ii1]:=Collected(irrepDimList);

# we select rows without zeros
charMatrTr1 := Filtered(charMatrTr, x -> ForAll(x, i -> i<>0));  

# we divide all rows from the irrep dimension:
charMatrTr2:=NullMat(DimensionsMat(charMatrTr1)[1], DimensionsMat(charMatrTr1)[2]);
for i1 in [1 .. Size(charMatrTr1)] do 
    for i2 in [1 .. Length(irrepDimList)] do 
          charMatrTr2[i1,i2]:=charMatrTr1[i1,i2]/irrepDimList[i2];
    od;
od;

# we discard the first line (trivial case)
# charMatrTr3 := Filtered(charMatrTr2, x -> ForAny(x, i -> i<>1));
# actually we don't need to throw out the first row, because if the rest of the classes have zeros, 
# then charMatrTr3 is empty (this happens in the case of SU(3))
charMatrTr3 := charMatrTr2;

# we select rows that have all elements with abs() = 1
charMatrTr4 := Filtered(charMatrTr3, x -> ForAll(x, i -> i*ComplexConjugate(i)=1));  

# we restore the original characters in the charMatrTr4 matrix
charMatrTr4org:=NullMat(DimensionsMat(charMatrTr4)[1], DimensionsMat(charMatrTr4)[2]);
for i1 in [1 .. Size(charMatrTr4)] do 
    for i2 in [1 .. Length(irrepDimList)] do 
          charMatrTr4org[i1,i2]:=charMatrTr4[i1,i2]*irrepDimList[i2];
    od;
od;

# we determine the positions of the classes in the full character matrix charMatrTr
posClassList:=[];
for i1 in [1 .. Size(charMatrTr4org)] do 
      posClassList[i1]:=Positions(charMatrTr,charMatrTr4org[i1]);
od;
posClassList:=Flat(posClassList);

# we make a power matrix
charMatrTr5:=NullMat(DimensionsMat(charMatrTr4)[1], DimensionsMat(charMatrTr4)[2]);
for i1 in [1 .. Size(charMatrTr4)] do 
    for i2 in [1 .. Length(charMatrTr4[i1])] do 
          pow:=0; sqVal:=0;
          while sqVal <> 1 do
              pow:=pow+1;
              sqVal:=charMatrTr4[i1,i2]^pow;
          od;
    charMatrTr5[i1,i2]:=pow;  
    od;
od;

# we select max powers
maxPowList:=[];
for i1 in [1 .. Size(charMatrTr5)] do 
      maxPowList[i1]:=Maximum(charMatrTr5[i1]);
od;
maxPow:=Maximum(maxPowList);
maxPowListFull[ii1]:=maxPow;

#----------- test 1 ---------------
# groups are divided by max powers
maxPowOrdGList:=[];
for i1 in [1 .. Length(maxPowList)] do 
      maxPowOrdGList[i1]:=ordG/maxPowList[i1];
od;

# we find positions of max powers
maxPowPosList:=[];
for i1 in [1 .. Size(charMatrTr5)] do 
      maxPowPosList[i1]:=Positions(charMatrTr5[i1],maxPowList[i1]);
od;

# we find characters with max powers
maxPowCharList:=[];
for i1 in [1 .. Size(charMatrTr4)] do 
      maxPowCharList[i1]:=DuplicateFreeList(charMatrTr4[i1]{maxPowPosList[i1]});
od;

# we find positions for characters with max powers
maxPowCharPosList:=[];
for i1 in [1 .. Size(charMatrTr4)] do 
tmpList:=[];
    for i2 in [1 .. Length(maxPowCharList[i1])] do 
      tmpList[i2]:=Positions(charMatrTr4[i1],maxPowCharList[i1,i2]);
    od;
    maxPowCharPosList[i1]:=tmpList;
od;

# sum the squares of irrep dimensions
irrepDimSqSumList:=[];
for i1 in [1 .. Size(charMatrTr4)] do 
tmpList:=[];
    for i2 in [1 .. Length(maxPowCharPosList[i1])] do 
      tmpList[i2]:=Sum(irrepDimList{maxPowCharPosList[i1,i2]}, x -> x^2);    
    od;
    irrepDimSqSumList[i1]:=DuplicateFreeList(tmpList);
od;
irrepDimSqSumList:=Flat(irrepDimSqSumList);

# we check whether the squares of irrep dimensions coincide with order of G divided by max power
irrepDimSqTest1:=irrepDimSqSumList = maxPowOrdGList;

# some checks:
# Print("Test:  ", irrepDimSqTest1, "\n");
# Print("order test:  ", maxPowList, "\n");
# Print("order test:  ", maxPowOrdGList, "\n");
# Print("order test:  ", irrepDimSqSumList, "\n");


#-----------------------------
# We are looking for irreps with the same characters across all selected groups
charMatr4:=TransposedMat(charMatrTr4);
equalCharList:=[];
for i1 in [1 .. Size(charMatr4)] do 
     equalCharList[i1]:=Length(DuplicateFreeList(charMatr4[i1]));
od;

# we find irrep positions with identical characters
posEqualCharList:=Positions(equalCharList,1);

# we select irreps with identical characters through all selected groups
irrepEqualCharList:=irrepDimList{posEqualCharList};

# we select the same characters through all the selected groups
charEqualCharList:=charMatrTr4[1]{posEqualCharList};

# we combine irrep dimensions with identical characters
irrepCharEqualCharList:=[];
for i1 in [1 .. Length(irrepEqualCharList)] do 
     irrepCharEqualCharList[i1]:=[irrepEqualCharList[i1],charEqualCharList[i1]];
od;

#----------- test 2 ---------------
# we select distinct characters
uniqueCharList:=DuplicateFreeList(charEqualCharList);

# we collect irreps dimensions with identical characters and calculate their sum of squares
tmpList:=[];
for i1 in [1 .. Size(uniqueCharList)] do 
    tmp:=Filtered(irrepCharEqualCharList, x -> x[2]=uniqueCharList[i1]);
    tmpList[i1]:=Sum(tmp, x -> x[1]^2);
od;
irrepDimSqSumList2:=DuplicateFreeList(tmpList);

# we check whether the squares of dimensions of irrep Equal coincide with the orders of G divided by max powers
irrepDimSqTest2:=Length(irrepDimSqSumList2)=1 and irrepDimSqSumList2[1]=ordG/maxPow;

# we make a list for summary output
infoList[ii1]:=[IdGroup(G), irrepDimSqTest1 and irrepDimSqTest2, maxPow];

# output to the terminal
Print( IdGroup(G), "    ", irrepDimSqTest1 and irrepDimSqTest2, "      ",  maxPow, "\n");

# output to the files
AppendTo(filename1, "Group:", "\n");
AppendTo(filename1, IdGroup(G), "\n");
AppendTo(filename1, "[nb. classes, nb. irreps]:", "\n");
AppendTo(filename1, DimensionsMat(charMatrTr), "\n");
AppendTo(filename1, "number of inequivalent irreps [dim(irrep), quantity]:", "\n");
AppendTo(filename1, dimCalcListFull[ii1], "\n");
AppendTo(filename1, "{max(p_1),..,max(p_n)} in each class:", "\n");
AppendTo(filename1, maxPowList, "\n");
AppendTo(filename1, "max(p):", "\n");
AppendTo(filename1, maxPow, "\n");
AppendTo(filename1, "[nb. classes with max(p_n), nb. classes with max(p)]:", "\n");
AppendTo(filename1, [Length(maxPowList), Length(Filtered(maxPowList, x -> x = maxPow))], "\n");
AppendTo(filename1, "positions of classes in characters table for max(p_n):", "\n");
AppendTo(filename1, posClassList, "\n");
AppendTo(filename1, "positions of irreps in characters table for max(p):", "\n");
AppendTo(filename1, maxPowPosList[Positions(maxPowList, maxPow)[1]], "\n");
AppendTo(filename1, "dimensions of irreps of max(p):", "\n");
AppendTo(filename1, irrepDimList{maxPowPosList[Positions(maxPowList, maxPow)[1]]}, "\n");
AppendTo(filename1, "inequivalent irreps which have the same character for all the classes [dim(irrep), sigma_{ij}]:", "\n");
AppendTo(filename1, irrepCharEqualCharList, "\n");
AppendTo(filename1, "test1 [sum(dim(irrep)^2) = order(G)/max(p)]:", "\n");
AppendTo(filename1, irrepDimSqTest1, "\n");
AppendTo(filename1, "test2 [sum(dim(irrep ineq.)^2) = order(G)/max(p)]:", "\n");
AppendTo(filename1, irrepDimSqTest2, "\n");

AppendTo(filename1, "---------------------------------------------------------------------", "\n\n");

od; # cycle through groups
#------------------------------------------------

AppendTo(filename1, JoinStringsWithSeparator(["Summary of groups ", String(groupNames[ii2]), ":"], ""), "\n");
AppendTo(filename1, "[Group,   test,   max(p)]", "\n");
AppendTo(filename1, "-------------------------------", "\n");
for i1 in [1 .. Length(nbGlist)] do 
     AppendTo(filename1, infoList[i1], "\n");
od;
#-----------------------------

AppendTo(filename3, "Nb. ",ii2, JoinStringsWithSeparator([",  groups: ", String(groupNames[ii2]), ""], ""), "\n");
AppendTo(filename3, "orders of groups:", "\n");
AppendTo(filename3, List(nbGlistFull[ii2], x -> x[1]), "\n");
AppendTo(filename3, "quantity of 1 dimensional irreps:", "\n");
AppendTo(filename3, List(dimCalcListFull, x -> x[1,2]), "\n");
AppendTo(filename3, "max(p) list:", "\n");
AppendTo(filename3, maxPowListFull, "\n");
AppendTo(filename3, "--------------------------", "\n\n");
#-----------------------------

# let's count how long it takes:
time3:=Runtime();
timeF:=time3-time0;
timeH:=Int(timeF/1000/60/60);
timeMin:=Int((timeF/1000-timeH*60*60)/60);
timeSec:=Int((timeF/1000-timeH*60*60-timeMin*60));
timeMsec:=Int((timeF/1000-timeH*60*60-timeMin*60-timeSec)*1000);

AppendTo(filename1, "\n");
AppendTo(filename1, "runTime: ", timeH, " h ", timeMin, " min ", timeSec, " sec ", timeMsec, " msec", "\n");
Print("runTime: ", timeH, " h ", timeMin, " min ", timeSec, " sec ", timeMsec, " msec", "\n");

od; # main loop through group sets

#------------------------------------------------------------------------------




