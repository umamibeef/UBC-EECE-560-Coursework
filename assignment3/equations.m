% Source voltage

syms vs delta_t t

vs = 230*cos(377*t)

% Circuit values

syms R1 R_line L1 L_line C1 C_line line_length Zc Zo a Rbrk frequency

R1 = 3
L1 = 350e-3
C1 = 10e-9

R_line = 0.5 % Ohm/km
L_line = 1e-3 % mH/km
C_line = 12e-9 % nF/km
line_length = 100 % km

frequency = 60 % Hz

Xl_line = 2*pi*frequency*L_line % Ohm/km
Xc_line = 1/(2*pi*frequency*C_line) % Ohm/km

% delta_t = 10e-6

% Velocity of propagation
a = sqrt(1/(L_line*C_line)) % km/s

% Transmission line history time delta
tau = line_length/a % s

% Characteristic line impedance
Zc = sqrt(L_line/C_line)

% Bergeron equivalent impedance
Zo = Zc + (R_line*line_length/4)
Hb = (Zo - (R_line*line_length/4))/(Zo + (R_line*line_length/4))

% Discretized Resistances:
R23 = (2*L1)/delta_t
R30 = delta_t/(2*C1)

% History Equations:
syms eh23 eh30
syms h23 h30 h40 h50
syms v2 v3 v4 v5
syms v4_tmintau v5_tmintau
syms v2_next v3_next v4_next v5_next
syms i23 i50 i30 i40
syms i40_tmintau i50_tmintau
syms i23_next i50_next i30_next i40_next

% History voltage source for L1 and C1
eh23 = -(v2 - v3) - ((2*L1)/delta_t)*i23
eh30 = v3 + (delta_t/(2*C1))*i30

% History current sources
h23 = eh23/R23
h30 = eh30/R30
h40 = ((1 + Hb)/2)*((-1/Zo)*v5_tmintau - Hb*i50_tmintau) + ((1 - Hb)/2)*((-1/Zo)*v4_tmintau - Hb*i40_tmintau)
h50 = ((1 + Hb)/2)*((-1/Zo)*v4_tmintau - Hb*i40_tmintau) + ((1 - Hb)/2)*((-1/Zo)*v5_tmintau - Hb*i50_tmintau)

i23_next = (v2_next - v3_next)/R23 - h23
i30_next = v3_next/R30 - h30
i40_next = v4_next/Zo + h40
i50_next = h50

syms g11 g12 g13 g14
syms g21 g22 g23 g24
syms g31 g32 g33 g34
syms g41 g42 g43 g44
syms h1 h2 h3 h4

h1 = 0
h2 = h23
h3 = h30 - h23
h4 = -h40

hA = [h2 ; h3 ; h4]
hB = h1
H = [hA ; hB]

g11 = 1/R1
g12 = -1/R1
g13 = 0
g14 = 0

g21 = -1/R1
g22 = 1/R1 + 1/R23
g23 = -1/R23
g24 = 0

g31 = 0
g32 = -1/R23
g33 = 1/R23 + 1/R30 + 1/Rbrk
g34 = -1/Rbrk

g41 = 0
g42 = 0
g43 = -1/Rbrk
g44 = 1/Rbrk + 1/Zo

GAA = [g22 g23 g24; g32 g33 g34; g42 g43 g44]
GAB = [g21 ; g31 ; g41]
GBA = [g12 g13 g14]
GBB = [g11]
G = [GAA GAB ; GBA GBB]

% Nodal voltage solution vector:
v_next = inv(GAA)*hA - inv(GAA)*GAB*vs
vpa(subs(v_next))
