% Source voltage

syms vs delta_t t

vs = 230*cos(377*t)

% Circuit values

syms R1 R_line L_line C_line line_length Zc Zo a Rbrk

R1 = 3

R_line = 0.5 % Ohm/km
L_line = 1e-3 % mH/km
C_line = 12e-9 % nF/km
line_length = 100 % km

delta_t = 10e-6

% Velocity of propagation
a = sqrt(1/(L_line*C_line)) % km/s

% Transmission line history time delta
tau = line_length/a % s

% Characteristic line impedance
Zc = sqrt(L_line/C_line)

% Bergeron equivalent impedance
Zo = Zc + (R_line*line_length/4)
Hb = (Zo - (R_line*line_length/4))/(Zo + (R_line*line_length/4))

% History Equations:
syms h23 h30
syms v2 v3
syms v2_tmintau v3_tmintau
syms v2_next v3_next
syms i20 i30
syms i20_tmintau i30_tmintau
syms i20_next i30_next

% History current sources
h20 = ((1 + Hb)/2)*((-1/Zo)*v3_tmintau - Hb*i30_tmintau) + ((1 - Hb)/2)*((-1/Zo)*v2_tmintau - Hb*i20_tmintau)
h30 = ((1 + Hb)/2)*((-1/Zo)*v2_tmintau - Hb*i20_tmintau) + ((1 - Hb)/2)*((-1/Zo)*v3_tmintau - Hb*i30_tmintau)

i20_next = v2_next/Zo + h20
i30_next = h30

syms g11 g12 g13
syms g21 g22 g23
syms g31 g32 g33
syms h1 h2 h3

h1 = 0
h2 = h20

hA = h2
hB = h1
H = [hA ; hB]

g11 = 1/R1
g12 = -1/R1

g21 = -1/R1
g22 = 1/R1 + 1/Zo

GAA = [g22]
GAB = [g21]
GBA = [g12]
GBB = [g11]
G = [GAA GAB ; GBA GBB]

% Nodal voltage solution vector:
v_next = inv(GAA)*hA - inv(GAA)*GAB*vs
vpa(subs(v_next))
