% Source voltage

syms vs delta_t t

vs = 10

% Circuit values

syms r_sw r l

% Discretized Resistances:
r30_trap = 2*l/delta_t
r30_back = l/delta_t

% History Equations:

syms eh30_trap eh30_back
syms h30_trap h30_back
syms v3
syms v3_next_trap
syms v3_next_back
syms flux flux_knee
syms i30
syms i30_trap_next i30_back_next

eh30_trap = -v3 - (2/delta_t)*flux + (2/delta_t)*flux_knee
eh30_back = -(1/delta_t)*flux + (1/delta_t)*flux_knee

h30_trap = eh30_trap/r30_trap
h30_back = eh30_back/r30_back

i30_trap_next = v3_next_trap/r30_trap - h30_trap
i30_back_next = v3_next_back/r30_back - h30_back

syms g11 g12 g13
syms g21 g22 g23
syms g31 g32 g33
syms h1 h2 h3

% Trapezoidal solution

h1 = 0
h2 = 0
h3 = h30_trap

hA = [h2 ; h3]
hB = h1
H = [hA ; hB]

g11 = 1/r_sw
g12 = -1/r_sw
g13 = 0

g21 = -1/r_sw
g22 = 1/r_sw + 1/r
g23 = -1/r

g31 = 0
g32 = -1/r
g33 = 1/r + 1/r30_trap

GAA = [g22 g23 ; g32 g33]
GAB = [g21 ; g31]
GBA = [g12 g13]
GBB = [g11]
G = [GAA GAB ; GBA GBB]

% Nodal voltage solution vector:

v_next_trap = inv(GAA)*hA - inv(GAA)*GAB*vs

% Backward Euler solution

h1 = 0
h2 = 0
h3 = h30_back

hA = [h2 ; h3]
hB = h1
H = [hA ; hB]

g11 = 1/r_sw
g12 = -1/r_sw
g13 = 0

g21 = -1/r_sw
g22 = 1/r_sw + 1/r
g23 = -1/r

g31 = 0
g32 = -1/r
g33 = 1/r + 1/r30_back

GAA = [g22 g23 ; g32 g33]
GAB = [g21 ; g31]
GBA = [g12 g13]
GBB = [g11]
G = [GAA GAB ; GBA GBB]

% Nodal voltage solution vector:

v_next_back = inv(GAA)*hA - inv(GAA)*GAB*vs

v_next_trap
v_next_back
i30_trap_next
i30_back_next