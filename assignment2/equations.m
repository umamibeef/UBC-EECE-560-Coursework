% Source voltage

syms vs delta_t t

vs = 230*cos(377*t)

% Circuit values

syms R1 R2 L1 L2 C1 C2

R1 = 3
L1 = 0.35
C1 = 10e-9
R2 = 50
L2 = 0.1
C2 = 600e-9
delta_t = 10e-6

% Discretized Resistances:
R23 = (2*L1)/delta_t
R50 = (2*L2)/delta_t
R30 = delta_t/(2*C1)
R40 = delta_t/(2*C2)

% History Equations:

syms eh23 eh50 eh30 eh40
syms h23 h50 h30 h40
syms v2 v3 v4 v5
syms v2_next v3_next v4_next v5_next
syms i23 i50 i30 i40
syms i23_next i50_next i30_next i40_next

eh23 = -(v2 - v3) - ((2*L1)/delta_t)*i23
eh50 = -v5 - ((2*L2)/delta_t)*i50
eh30 = v3 + (delta_t/(2*C1))*i30
eh40 = v4 + (delta_t/(2*C2))*i40

h23 = eh23/R23
h50 = eh50/R50
h30 = eh30/R30
h40 = eh40/R40

i23_next = (v2_next - v3_next)/R23 - h23
i30_next = v3_next/R30 - h30
i40_next = v4_next/R40 - h40
i50_next = v5_next/R50 - h50

% Switch Closed (v3 = v4):

syms g11 g12 g13 g14 g15
syms g21 g22 g23 g24 g25
syms g31 g32 g33 g34 g55
syms g41 g42 g43 g44 g45
syms g51 g52 g53 g54 g55
syms h1 h2 h3 h4 h5

h1 = 0
h2 = h23
h3 = -h23 + h30 + h40
h5 = h50

hA = [h2 ; h3 ; h5]
hB = h1
H = [hA ; hB]

g11 = 1/R1
g12 = -1/R1
g13 = 0
g15 = 0

g21 = -1/R1
g22 = 1/R1 + 1/R23
g23 = -1/R23
g25 = 0

g31 = 0
g32 = -1/R23
g33 = 1/R23 + 1/R30 + 1/R40 + 1/R2
g35 = -1/R2

g51 = 0
g52 = 0
g53 = -1/R2
g55 = 1/R50 + 1/R2

GAA = [g22 g23 g25 ; g32 g33 g35 ; g52 g53 g55]
GAB = [g21 ; g31 ; g51]
GBA = [g12 g13 g15]
GBB = [g11]
G = [GAA GAB ; GBA GBB]

% Nodal voltage solution vector:

v_next = inv(GAA)*hA - inv(GAA)*GAB*vs
vpa(subs(v_next))

% Switch Open (v3 != v4):

syms g11 g12 g13 g14 g15
syms g21 g22 g23 g24 g25
syms g31 g32 g33 g34 g55
syms g41 g42 g43 g44 g45
syms g51 g52 g53 g54 g55
syms h1 h2 h3 h4 h5

h1 = 0
h2 = h23
h3 = -h23 + h30
h4 = h40
h5 = h50

hA = [h2 ; h3 ; h4 ; h5]
hB = h1
H = [hA ; hB]

g11 = 1/R1
g12 = -1/R1
g13 = 0
g14 = 0
g15 = 0

g21 = -1/R1
g22 = 1/R1 + 1/R23
g23 = -1/R23
g24 = 0
g25 = 0

g31 = 0
g32 = -1/R23
g33 = 1/R23 + 1/R30
g34 = 0
g35 = 0

g41 = 0
g42 = 0
g43 = 0
g44 = 1/R40 + 1/R2
g45 = -1/R2

g51 = 0
g52 = 0
g53 = 0
g54 = -1/R2
g55 = 1/R50 + 1/R2

GAA = [g22 g23 g24 g25 ; g32 g33 g34 g35 ; g42 g43 g44 g45 ; g52 g53 g54 g55]
GAB = [g21 ; g31 ; g41 ; g51]
GBA = [g12 g13 g14 g15]
GBB = [g11]
G = [GAA GAB ; GBA GBB]

% Nodal voltage solution vector:

v_next = inv(GAA)*hA - inv(GAA)*GAB*vs
vpa(subs(v_next))

% Resonant frequency approximation

% f1 = 1/(2*pi*sqrt(L1*C1)) = 2.69 kHz => t1 = 372 us
% f2 = 1/(2*pi*sqrt(L2*C2)) = 650 Hz
% f3 = 1/(2*pi*sqrt((L1+L2)*(C1+C2))) = 303 Hz
% choose 10 us time step
