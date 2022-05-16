// Gmsh project created on Wed Apr 27 20:30:50 2022
// Add points 
Point(1) = {0, 0, 0, 1};
Point(2) = {45, 0, 0, 1};
Point(3) = {45, 14, 0, 1};
Point(4) = {45, 28, 0, 1};
Point(5) = {0, 28, 0, 1};
Point(6) = {0, 14, 0, 1};
Point(7) = {63, 0, 0, 1};
Point(8) = {63, 14, 0, 0.1};
Point(9) = {63, 28, 0, 0.1};
Point(10) = {72, 14, 0, 0.1};
Point(11) = {72, 28, 0, 0.1};
Point(12) = {119, 14, 0, 1};
Point(13) = {119, 28, 0, 1};
Point(14) = {126, 14, 0, 1};
Point(15) = {126, 28, 0, 1};

// Add lines 
Line(1) = {1, 2}; 
Line(2) = {2, 3}; 
Line(3) = {3, 6}; 
Line(4) = {6, 1}; 
Line(5) = {3, 4}; 
Line(6) = {4, 5}; 
Line(7) = {5, 6}; 
Line(8) = {2, 7}; 
Line(9) = {7, 8}; 
Line(10) = {8, 3}; 
Line(11) = {8, 9}; 
Line(12) = {9, 4}; 
Line(13) = {8, 10}; 
Line(14) = {10, 11}; 
Line(15) = {11, 9}; 
Line(16) = {10, 12}; 
Line(17) = {12, 13}; 
Line(18) = {13, 11}; 
Line(19) = {12, 14}; 
Line(20) = {14, 15}; 
Line(21) = {15, 13}; 

Transfinite Line {-8, 10, 12} = 101 Using Progression 1.02;
Transfinite Line {13, -15} = 101 Using Progression 1;
Transfinite Line {16, -18} = 101 Using Progression 1.02;
Transfinite Line {1, 3, 6, 19, 21} = 21 Using Progression 1;
Transfinite Line {-2, 4, 5, -7, -9, 11, 14, 17, 20} = 81 Using Progression 1.02;

// Add line loops and surfaces 
Line Loop(1) = {1, 2, 3, 4};
Ruled Surface(1) = {1};
Transfinite Surface {1} = {1, 2, 3, 6};
Recombine Surface {1};

Line Loop(2) = {8, 9, 10, -2};
Ruled Surface(2) = {2};
Transfinite Surface {2} = {2, 7, 8, 3};
Recombine Surface {2};

Line Loop(3) = {-3, 5, 6, 7};
Ruled Surface(3) = {3};
Transfinite Surface {3} = {6, 3, 4, 5};
Recombine Surface {3};

Line Loop(4) = {-10, 11, 12, -5};
Ruled Surface(4) = {4};
Transfinite Surface {4} = {3, 8, 9, 4};
Recombine Surface {4};

Line Loop(5) = {13, 14, 15, -11};
Ruled Surface(5) = {5};
Transfinite Surface {5} = {8, 10, 11, 9};
Recombine Surface {5};

Line Loop(6) = {16, 17, 18, -14};
Ruled Surface(6) = {6};
Transfinite Surface {6} = {10, 12, 13, 11};
Recombine Surface {6};

Line Loop(7) = {19, 20, 21, -17};
Ruled Surface(7) = {7};
Transfinite Surface {7} = {12, 14, 15, 13};
Recombine Surface {7};

// Extrude to volumes 
Extrude {0, 0, 14} {Surface{1, 2, 3, 4, 5, 6, 7}; Layers{2}; Recombine;}

// Add physical sets 
Physical Volume("all") = {1, 2, 3, 4, 5, 6, 7};