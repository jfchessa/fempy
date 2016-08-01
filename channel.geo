// Gmsh project created on Mon Jul 18 12:12:00 2016
lc = DefineNumber[ 0.25, Name "Parameters/lc" ];
ri = 5;
tw = 0.25;
ch = 1.75;
ow = 0.5;
phi = 15*Pi/180;
cratio = 0.65;

Point(1) = {0, 0, 0, lc};
Point(2) = {ri, 0, 0, lc};
Point(3) = {ri+tw, 0, 0, lc};
Point(4) = {ri+tw+ch, 0, 0, lc};
Point(5) = {ri+tw+ch+ow, 0, 0, lc};

Point(6) = { ri*Cos(phi), ri*Sin(phi), 0, lc };
Point(7) = { (ri+tw+ch+ow)*Cos(phi), (ri+tw+ch+ow)*Sin(phi), 0, lc };

Point(12) = { (ri+tw)*Cos(phi), (ri+tw)*Sin(phi), 0, lc };
Point(13) = { (ri+tw+ch)*Cos(phi), (ri+tw+ch)*Sin(phi), 0, lc };

phi = cratio*phi;
Point(8) = { ri*Cos(phi), ri*Sin(phi), 0, lc };
Point(9) = { (ri+tw)*Cos(phi), (ri+tw)*Sin(phi), 0, lc };
Point(10) = { (ri+tw+ch)*Cos(phi), (ri+tw+ch)*Sin(phi), 0, lc };
Point(11) = { (ri+tw+ch+ow)*Cos(phi), (ri+tw+ch+ow)*Sin(phi), 0, lc };

Circle(1) = {2, 1, 8};
Circle(2) = {8, 1, 6};

Circle(3) = {3, 1, 9};
Circle(4) = {4, 1, 10};
Circle(5) = {5, 1, 11};
Circle(6) = {11, 1, 7};
Line(7) = {2, 3};
Line(8) = {3, 4};
Line(9) = {4, 5};
Line(10) = {10, 9};
Line(11) = {6, 12};
Line(14) = {12,13};
Line(15) = {13,7};

Line(12) = {8, 9};
Line(13) = {10, 11};
Circle(16) = {9, 1, 12};
Circle(17) = {10, 1, 13};
Transfinite Line {1, 3, 4, 5} = 10 Using Progression 1;
Transfinite Line {2, 16, 17, 6} = 5 Using Progression 1;
Transfinite Line {7, 12, 11} = 5 Using Progression 1;
Transfinite Line {8, 10, 14} = 10 Using Progression 1;
Transfinite Line {9, 13, 15} = 4 Using Progression 1;


Line Loop(18) = {7, 3, -12, -1};
Plane Surface(19) = {18};
Line Loop(20) = {8, 4, 10, -3};
Plane Surface(21) = {20};
Line Loop(22) = {9, 5, -13, -4};
Plane Surface(23) = {22};
Line Loop(24) = {2, 11, -16, -12};
Plane Surface(25) = {24};
Line Loop(26) = {10, 16, 14, -17};
Plane Surface(27) = {26};
Line Loop(28) = {13, 6, -15, -17};
Plane Surface(29) = {28};
Transfinite Surface {19} = {2, 3, 9, 8};
Transfinite Surface {25} = {8, 9, 12, 6};
Transfinite Surface {21} = {3, 4, 10, 9};
Transfinite Surface {27} = {9, 10, 13, 12};
Transfinite Surface {23} = {4, 5, 11, 10};
Transfinite Surface {29} = {10, 11, 7, 13};
Recombine Surface {19, 21, 23, 29, 27, 25};
Extrude {0, 0, 4} {
  Surface{19, 21, 23, 29, 27, 25}; Layers{10}; Recombine;
}
Physical Volume(162) = {2};
Physical Volume(163) = {1, 6, 5, 4, 3};
Physical Surface(164) = {21};
Physical Surface(165) = {148, 50};
