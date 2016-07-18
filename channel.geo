// Gmsh project created on Mon Jul 18 12:12:00 2016
lc = DefineNumber[ 0.1, Name "Parameters/lc" ];
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
Line(11) = {6, 7};

Line Loop(12) = {3, -10, -4, -8};
Plane Surface(13) = {12};
Line Loop(14) = {7, 3, -10, -4, 9, 5, 6, -11, -2, -1};
Plane Surface(15) = {14};
Physical Surface(16) = {13};
Physical Surface(17) = {15};
Physical Line(18) = {11, 9, 8, 7};
Physical Line(19) = {1, 2};
Physical Line(20) = {5, 6};
Recombine Surface {15, 13};
