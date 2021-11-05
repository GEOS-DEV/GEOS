// Gmsh project created on Wed Nov 11 14:32:31 2020
ta = 1.0; //0.1;
a = 0.5;
c = 0.005;
l = 1;
h = 0.5;
Point(1) = {0.0,0.0,0.0,ta};
Point(2) = {l, 0, 0.0, ta};
Point(3) = {l, h, 0.0, ta};
Point(4) = {l, l , 0.0, ta};
Point(5) = {0.0, l, 0.0, ta};
Point(6) = {0.0, h + c, 0.0, ta};
Point(7) = {0.0, h - c, 0.0, ta};
Point(8) = {a, h, 0.0, ta};
Point(9) = {0.5, 0, 0.0, ta};
Point(10) = {0.5, l, 0.0, ta};

Line(1) = {1, 9};
Line(2) = {9, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 10};
Line(6) = {10, 5};
Line(7) = {5, 6};
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 1};
Line(11) = {8, 10};
Line(12) = {8, 9};
Line(13) = {8, 3};

Transfinite Line {1, 8, 9, 6} = 4 Using Progression 1;
Transfinite Line {7, 10, 11, 12, 3, 4} = 4 Using Progression 1;
Transfinite Line {2, 5, 13} = 4 Using Progression 1;

Line Loop(14) = {1, -12, 9, 10};
Ruled Surface(15) = {14};
Transfinite Surface {15} = {1, 9, 8, 7};
Recombine Surface {15};

Line Loop(16) = {8, 11, 6, 7};
Ruled Surface(17) = {16};
Transfinite Surface {17} = {8, 10, 5, 6};
Recombine Surface {17};

Line Loop(18) = {2, 3, -13, 12};
Ruled Surface(19) = {18};
Transfinite Surface {19} = {9, 2, 3, 8};
Recombine Surface {19};

Line Loop(20) = {13, 4, 5, -11};
Ruled Surface(21) = {20};
Transfinite Surface {21} = {8, 3, 4, 10};
Recombine Surface {21};

//Physical Surface(23) = {17, 15, 19, 21};
Extrude {0, 0, 0.1} {Surface{15, 17, 19, 21}; Layers{1}; Recombine;}
//Extrude {0, 0, 0.2} {Surface{23}; Layers{1}; Recombine;}
//Compound Volume(5) = {1, 2, 3, 4};
Physical Volume("all") = {1, 2, 3, 4};
Physical Surface("top") = {60, 104};
Physical Surface("bottom") = {30, 74};
Physical Surface("left") = {42, 64};
Physical Surface("right") = {78, 100};
//Physical Line(22) = {6, 7, 8, 10, 11, 1, 2, 3, 4, 5};
//Physical Surface(23) = {17, 15, 19, 21};