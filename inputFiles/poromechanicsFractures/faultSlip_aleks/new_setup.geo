W = 4500.0;
H = 4500.0;
a = 75.0;
D = H / 2;
Y1 = D - H;
Y2 = D;
b = 150;
Lplus = b + 100;
phi = 90 * Pi / 180;
lc = 150;
mult = 1.0;
mult1 = 0.1;//0.005;

Point(1) = {-W/2, Y1, 0, lc};
Point(2) = {W/2, Y1, 0, lc};
Point(3) = {W/2, -a, 0, mult * lc};
Point(4) = {W/2, b, 0, mult * lc};
Point(5) = {W/2, Y2, 0, lc};
Point(6) = {-W/2, Y2, 0, lc};
Point(7) = {-W/2, a, 0, mult * lc};
Point(8) = {-W/2, -b, 0, mult * lc};

Point(15) = {-b / Tan(phi), -b, 0, mult1 * lc};
Point(16) = {-a / Tan(phi), -a, 0, mult1 * lc};
Point(17) = {a / Tan(phi), a, 0, mult1 * lc};
Point(18) = {b / Tan(phi), b, 0, mult1 * lc};
Point(19) = {Lplus / Tan(phi), Lplus, 0, mult1 * lc};
Point(20) = {-Lplus / Tan(phi), -Lplus, 0, mult1 * lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,1};

Line(9) = {3, 16};
Line(10) = {4, 18};
Line(11) = {8, 15};
Line(12) = {7, 17};
Line(13) = {18, 17};
Line(14) = {17, 16};
Line(15) = {16, 15};
Line(16) = {15, 20};
Line(17) = {18, 19};

Line Loop(1) = {1, 2, 9, 15, -11, 8};
Line Loop(2) = {3, 10, 13, 14, -9};
Line Loop(3) = {-12,7,11,-15,-14};
Line Loop(4) = {4, 5, 6, 12, -13, -10};
Line Loop(5) = {16, -16};
Line Loop(6) = {17, -17};

Plane Surface(1) = {1, -5};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4, -6};

out[] = Extrude {0, 0, 500} {
	 Surface{:};
	 Layers{1};
	 Recombine;
	};

LEFT = 991;
RIGHT = 992;
BOT = 993;
TOP = 994;
ZM = 995;
ZP = 996;
RES = 99991;
//OUTER = 99992;
OVER = 99992;
UNDER = 99993;
FRAC = 9991;
FRAC_BOUND_STICK = 1;
FRAC_BOUND_FREE = 2;

Physical Volume(RES) = {out[11], out[18]};
Physical Volume(UNDER) = { out[1] };
Physical Volume(OVER) = { out[25] };
Physical Surface(ZM) = {1,2,3,4};
Physical Surface(ZP) = {out[0], out[10], out[17], out[24]};
Physical Surface(BOT) = {out[2]};
Physical Surface(RIGHT) = {out[3], out[12], out[26]};
Physical Surface(TOP) = {out[27]};
Physical Surface(LEFT) = {out[7], out[20], out[28]};
Physical Surface(FRAC) = {out[32], out[14], out[15], out[22], out[9]};//, out[4], out[12]};
//Physical Curve(FRAC_BOUND_STICK) = {53, 149};
//Physical Curve(FRAC_BOUND_FREE) = {121, 17, 13, 14, 15, 16, 25, 22, 64, 63}; //,38,66,78,10,13,6,43};

Mesh 3;
Coherence Mesh;
Mesh.MshFileVersion = 2.1;                                                     numerable IDisposable GetPreamble RuntimeFieldHandle nStdHandle GetStdHandle SafeHandle findFileHandle SafeFileHandle moduleHandle RuntimeTypeHandle CloseHandle G