RO = 1;                 // outer radius 
RI = 0.4;               // inner radius
AL = 2*Pi/4;            // arclength
D = 0.02*2*Pi*RI;       // mesh density
H = 0.5;                 // Extension in Z direction
DR = 0.1;             // Radius of the domain of observation

CONTDOMS = 0;           // Start of range for control domains
OBSDOMS = 10;          // Start of range for observation domains
VOLDOMS = 20;          // Start of range for physical volumes

Point(1) = {0, 0, 0, D};  // Center point
Point(2) = {0, 0, H, D};  // Center point up
Point(3) = {0, 0, H/2, D};  // Center point middle

// *************
// STARTING FACE
// *************

ppp = 3;
Point(ppp+1) = {RI, 0, 0, D};  // inn bot
Point(ppp+2) = {RO, 0, 0, D};  // out bot
Point(ppp+3) = {RO, 0, H/2, D};  // out mid
Point(ppp+4) = {RO, 0, H, D};  // out top
Point(ppp+5) = {RI+DR, 0, H, D};  // mid top
Point(ppp+6) = {RI, 0, H, D};  // inn top
Point(ppp+7) = {RI, 0, H/2, D};  // inn mid

Line(1) = {ppp+1, ppp+2};
Line(2) = {ppp+2, ppp+3};
Line(3) = {ppp+3, ppp+4};
Line(4) = {ppp+4, ppp+5};
Line(5) = {ppp+5, ppp+6};
Line(6) = {ppp+6, ppp+7};
Line(7) = {ppp+7, ppp+1};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

lpp = 7;
spp = 1;
ppp = ppp+7;

// *************
// FIRST SEGMENT
// *************

segid = 1;

Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-6}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-5}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-4}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-3}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-2}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-1}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp}; }}

Circle(lpp+1) = {ppp-6, 1, ppp+1};
Circle(lpp+2) = {ppp-5, 1, ppp+2};
Circle(lpp+3) = {ppp-4, 3, ppp+3};
Circle(lpp+4) = {ppp-3, 2, ppp+4};
Circle(lpp+5) = {ppp-2, 2, ppp+5};
Circle(lpp+6) = {ppp-1, 2, ppp+6};
Circle(lpp+7) = {ppp, 3, ppp+7};

lpp = lpp + 7;

Line(lpp+1) = {ppp+1, ppp+2};
Line(lpp+2) = {ppp+2, ppp+3};
Line(lpp+3) = {ppp+3, ppp+4};
Line(lpp+4) = {ppp+4, ppp+5};
Line(lpp+5) = {ppp+5, ppp+6};
Line(lpp+6) = {ppp+6, ppp+7};
Line(lpp+7) = {ppp+7, ppp+1};

Printf('lpp: %g', lpp);

Line Loop(spp+1) = {lpp-7-6, lpp-5, -lpp-1, -lpp+6};  // bottom
Plane Surface(spp+1) = {spp+1};
Line Loop(spp+2) = {lpp-7-5, lpp-4, -lpp-2, -lpp+5};  // outer lower
Surface(spp+2) = {spp+2};
Line Loop(spp+3) = {lpp-7-4, lpp-3, -lpp-3, -lpp+4};  // outer upper
Surface(spp+3) = {spp+3};
Line Loop(spp+4) = {lpp-7-3, lpp-2, -lpp-4, -lpp+3};  // ceiling outer
Plane Surface(spp+4) = {spp+4};
Line Loop(spp+5) = {lpp-7-2, lpp-1, -lpp-5, -lpp+2};  // ceiling inner
Plane Surface(spp+5) = {spp+5};
Line Loop(spp+6) = {lpp-7-1, lpp  , -lpp-6, -lpp+1};  // inner upper
Surface(spp+6) = {spp+6};
Line Loop(spp+7) = {lpp-7  , lpp-6, -lpp-7, -lpp  };  // inner lower
Surface(spp+7) = {spp+7};

Line Loop(spp+8) = {lpp+1, lpp+2, lpp+3, lpp+4, lpp+5, lpp+6, lpp+7};
Plane Surface(spp+8) = {spp+8};

Physical Surface(0) = {1};              // A DUMMY FOR FENICS
Physical Surface(CONTDOMS+segid) = {spp+1};
Physical Surface(OBSDOMS+segid) = {spp+5};

Surface Loop(segid) = {1, spp+1, spp+2, spp+3, spp+4, spp+5, spp+6, spp+7, spp+8};
Volume(segid) = {segid};
Physical Volume(VOLDOMS+segid) = {segid};

lpp = lpp+7;
spp = spp+8;
ppp = ppp+7;

// *************
// SECOND SEGMENT
// *************

segid = 2;

Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-6}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-5}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-4}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-3}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-2}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-1}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp}; }}

Circle(lpp+1) = {ppp-6, 1, ppp+1};
Circle(lpp+2) = {ppp-5, 1, ppp+2};
Circle(lpp+3) = {ppp-4, 3, ppp+3};
Circle(lpp+4) = {ppp-3, 2, ppp+4};
Circle(lpp+5) = {ppp-2, 2, ppp+5};
Circle(lpp+6) = {ppp-1, 2, ppp+6};
Circle(lpp+7) = {ppp, 3, ppp+7};

lpp = lpp + 7;

Line(lpp+1) = {ppp+1, ppp+2};
Line(lpp+2) = {ppp+2, ppp+3};
Line(lpp+3) = {ppp+3, ppp+4};
Line(lpp+4) = {ppp+4, ppp+5};
Line(lpp+5) = {ppp+5, ppp+6};
Line(lpp+6) = {ppp+6, ppp+7};
Line(lpp+7) = {ppp+7, ppp+1};

Printf('lpp: %g', lpp);

Line Loop(spp+1) = {lpp-7-6, lpp-5, -lpp-1, -lpp+6};  // bottom
Plane Surface(spp+1) = {spp+1};
Line Loop(spp+2) = {lpp-7-5, lpp-4, -lpp-2, -lpp+5};  // outer lower
Surface(spp+2) = {spp+2};
Line Loop(spp+3) = {lpp-7-4, lpp-3, -lpp-3, -lpp+4};  // outer upper
Surface(spp+3) = {spp+3};
Line Loop(spp+4) = {lpp-7-3, lpp-2, -lpp-4, -lpp+3};  // ceiling outer
Plane Surface(spp+4) = {spp+4};
Line Loop(spp+5) = {lpp-7-2, lpp-1, -lpp-5, -lpp+2};  // ceiling inner
Plane Surface(spp+5) = {spp+5};
Line Loop(spp+6) = {lpp-7-1, lpp  , -lpp-6, -lpp+1};  // inner upper
Surface(spp+6) = {spp+6};
Line Loop(spp+7) = {lpp-7  , lpp-6, -lpp-7, -lpp  };  // inner lower
Surface(spp+7) = {spp+7};

Line Loop(spp+8) = {lpp+1, lpp+2, lpp+3, lpp+4, lpp+5, lpp+6, lpp+7};
Plane Surface(spp+8) = {spp+8};

Physical Surface(CONTDOMS+segid) = {spp+1};
Physical Surface(OBSDOMS+segid) = {spp+5};

Surface Loop(segid) = {spp, spp+1, spp+2, spp+3, spp+4, spp+5, spp+6, spp+7, spp+8};
Volume(segid) = {segid};
Physical Volume(VOLDOMS+segid) = {segid};

lpp = lpp+7;
spp = spp+8;
ppp = ppp+7;

// *************
// THIRD SEGMENT
// *************

segid = 3;

Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-6}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-5}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-4}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-3}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-2}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp-1}; }}
Rotate { {0,0,1}, {0, 0, 0},  AL } { Duplicata { Point{ppp}; }}

Circle(lpp+1) = {ppp-6, 1, ppp+1};
Circle(lpp+2) = {ppp-5, 1, ppp+2};
Circle(lpp+3) = {ppp-4, 3, ppp+3};
Circle(lpp+4) = {ppp-3, 2, ppp+4};
Circle(lpp+5) = {ppp-2, 2, ppp+5};
Circle(lpp+6) = {ppp-1, 2, ppp+6};
Circle(lpp+7) = {ppp, 3, ppp+7};

lpp = lpp + 7;

Line(lpp+1) = {ppp+1, ppp+2};
Line(lpp+2) = {ppp+2, ppp+3};
Line(lpp+3) = {ppp+3, ppp+4};
Line(lpp+4) = {ppp+4, ppp+5};
Line(lpp+5) = {ppp+5, ppp+6};
Line(lpp+6) = {ppp+6, ppp+7};
Line(lpp+7) = {ppp+7, ppp+1};

Printf('lpp: %g', lpp);

Line Loop(spp+1) = {lpp-7-6, lpp-5, -lpp-1, -lpp+6};  // bottom
Plane Surface(spp+1) = {spp+1};
Line Loop(spp+2) = {lpp-7-5, lpp-4, -lpp-2, -lpp+5};  // outer lower
Surface(spp+2) = {spp+2};
Line Loop(spp+3) = {lpp-7-4, lpp-3, -lpp-3, -lpp+4};  // outer upper
Surface(spp+3) = {spp+3};
Line Loop(spp+4) = {lpp-7-3, lpp-2, -lpp-4, -lpp+3};  // ceiling outer
Plane Surface(spp+4) = {spp+4};
Line Loop(spp+5) = {lpp-7-2, lpp-1, -lpp-5, -lpp+2};  // ceiling inner
Plane Surface(spp+5) = {spp+5};
Line Loop(spp+6) = {lpp-7-1, lpp  , -lpp-6, -lpp+1};  // inner upper
Surface(spp+6) = {spp+6};
Line Loop(spp+7) = {lpp-7  , lpp-6, -lpp-7, -lpp  };  // inner lower
Surface(spp+7) = {spp+7};

Line Loop(spp+8) = {lpp+1, lpp+2, lpp+3, lpp+4, lpp+5, lpp+6, lpp+7};
Plane Surface(spp+8) = {spp+8};

Physical Surface(CONTDOMS+segid) = {spp+1};
Physical Surface(OBSDOMS+segid) = {spp+5};

Surface Loop(segid) = {spp, spp+1, spp+2, spp+3, spp+4, spp+5, spp+6, spp+7, spp+8};
Volume(segid) = {segid};
Physical Volume(VOLDOMS+segid) = {segid};

lpp = lpp+7;
spp = spp+8;
ppp = ppp+7;

// *************
// FOURTH SEGMENT
// *************

segid = 4;

Printf('ppp: %g', ppp);

Circle(lpp+1) = {ppp-6, 1, 3+1};
Circle(lpp+2) = {ppp-5, 1, 3+2};
Circle(lpp+3) = {ppp-4, 3, 3+3};
Circle(lpp+4) = {ppp-3, 2, 3+4};
Circle(lpp+5) = {ppp-2, 2, 3+5};
Circle(lpp+6) = {ppp-1, 2, 3+6};
Circle(lpp+7) = {ppp  , 3, 3+7};

lpp = lpp+7;

Printf('spp: %g', spp);
Printf('lpp: %g', lpp);

Line Loop(spp+1) = {lpp-7-6, lpp-5, -1, -lpp+6};  // bottom
Plane Surface(spp+1) = {spp+1};
Line Loop(spp+2) = {lpp-7-5, lpp-4, -2, -lpp+5};  // outer lower
Surface(spp+2) = {spp+2};
Line Loop(spp+3) = {lpp-7-4, lpp-3, -3, -lpp+4};  // outer upper
Surface(spp+3) = {spp+3};
Line Loop(spp+4) = {lpp-7-3, lpp-2, -4, -lpp+3};  // ceiling outer
Plane Surface(spp+4) = {spp+4};
Line Loop(spp+5) = {lpp-7-2, lpp-1, -5, -lpp+2};  // ceiling inner
Plane Surface(spp+5) = {spp+5};
Line Loop(spp+6) = {lpp-7-1, lpp  , -6, -lpp+1};  // inner upper
Surface(spp+6) = {spp+6};
Line Loop(spp+7) = {lpp-7  , lpp-6, -7, -lpp  };  // inner lower
Surface(spp+7) = {spp+7};

Physical Surface(CONTDOMS+segid) = {spp+1};
Physical Surface(OBSDOMS+segid) = {spp+5};

Surface Loop(segid) = {spp, spp+1, spp+2, spp+3, spp+4, spp+5, spp+6, spp+7, 1};
Volume(segid) = {segid};
Physical Volume(VOLDOMS+segid) = {segid};

Field[1] = Box;
Field[1].VIn = D/2;
Field[1].VOut = D;
Field[1].XMax = 1;
Field[1].XMin = -1;
Field[1].YMax = 1;
Field[1].YMin = -1;
Field[1].ZMax = 0.5*D;
Field[1].ZMin = -0.01;

// Field[2] = Box;
// Field[2].VIn = D/2;
// Field[2].VOut = D;
// Field[2].XMax = RI+DR;
// Field[2].XMin = -RI-DR;
// Field[2].YMax = RI+DR;
// Field[2].YMin = -RI-DR;
// Field[2].ZMax = H;
// Field[2].ZMin = H-D;

Field[2] = Attractor;
Field[2].NNodesByEdge = 50;
Field[2].EdgesList = {11,25,39,53,13,27,41,55};
Field[2].NodesList = {7,14,21,28,9,16,23,30};

Field[3] = Threshold;
Field[3].IField = 2;
Field[3].LcMin = D/3;
Field[3].LcMax = D;
Field[3].DistMin = D/3;
Field[3].DistMax = D/2;

Field[4] = Min;
Field[4].FieldsList = {1,3};
Background Field = 4;
