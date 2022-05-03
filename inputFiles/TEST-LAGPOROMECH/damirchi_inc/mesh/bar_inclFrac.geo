lc = 1.00; // target mesh size
th = 2.00; // thickness of the PARTIAL extrusion 

Point(1) = { 0, 0, 0, lc};
Point(2) = { 2.17, 0, 0, lc};
Point(3) = {10, 0, 0, lc};
Point(4) = {10, 0, 5.17, lc};
Point(5) = {10, 0, 10.83, lc};
Point(6) = {10, 0, 16, lc};
Point(7) = { 7.83, 0, 16, lc};
Point(8) = { 0, 0, 16, lc};
Point(9) = { 0, 0, 10.83, lc};
Point(10) = {0, 0, 5.17, lc};
Point(11) = {7.83, 0, 10.83, lc};
Point(12) = {2.17, 0, 5.17, lc};

Line(1) = {1 , 2};
Line(2) = {2 , 3};
Line(3) = {1 , 10};
Line(4) = {2 , 12};
Line(5) = {3 , 4};
Line(6) = {10 , 12};
Line(7) = {12 , 4};
Line(8) = {10 , 9};
Line(9) = {12 , 11};
Line(10) = {4 , 5};
Line(11) = {9 , 11};
Line(12) = {11 , 5};
Line(13) = {9 , 8};
Line(14) = {11 , 7};
Line(15) = {5 , 6};
Line(16) = {8 , 7};
Line(17) = {7 , 6};

Curve Loop (1) = {1 , 4 , -6 , -3}; Plane Surface(1) = {1};
Curve Loop (2) = {2 , 5 , -7 , -4}; Plane Surface(2) = {2};
Curve Loop (3) = {7 , 10 , -12 , -9}; Plane Surface(3) = {3};
Curve Loop (4) = {6 , 9 , -11 , -8}; Plane Surface(4) = {4};
Curve Loop (5) = {11 , 14 , -16 , -13}; Plane Surface(5) = {5};
Curve Loop (6) = {12 , 15 , -17 , -14}; Plane Surface(6) = {6};

Extrude {0, 10*th, 0} { 
   Surface{1, 2, 3, 4, 5, 6 };
   Layers{10,1};
 }

Physical Volume("DOMAIN",1) = {1, 2, 3, 4, 5, 6 };

Physical Surface("FRACTURE_0") = { 82 };
Physical Surface("ZMAX") = {122, 144};
Physical Surface("ZMIN") = {26, 48};
Physical Surface("XMIN") = {38, 104, 126};
Physical Surface("XMAX") = {52, 74, 140};
//Physical Surface("DOMAIN_BOUNDAY_SURFACE_4") = {303, 325, 347, 369, 391, 413};
Physical Surface("YMAX") = {39, 61, 83, 105, 127, 149};
Physical Surface("YMIN") = {1, 2, 3, 4, 5, 6};

Mesh 3;
Mesh.MshFileVersion = 2.2;
Mesh.IgnorePeriodicity = 2;
//Save "bar_1_frac_ref.msh";
