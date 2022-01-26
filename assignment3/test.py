import argparse
import csv
from math import cos, sqrt
import matplotlib
import numpy as np
 
 # Matplotlib export settings
matplotlib.use("pgf")
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.size": 10 ,
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False     # don't setup fonts from rc parameters
})

# Interpolation function
# t0 < t < t1
def interpolate(t, t0, y0, t1, y1):
    # Linear interpolation between points (t0, y0) and (t1, y1)
    return y0 + (t - t0)*(y1 - y0)/(t1 - t0)

# Main function
def main(args):
    # Lists for results
    pscad_t_vals = []
    pscad_v1_vals = []
    pscad_v2_vals = []
    pscad_v3_vals = []
    pscad_ishort_vals = []

    # tuple indices
    V1_INDEX = 0
    V2_INDEX = 1
    V3_INDEX = 2
    I20_INDEX = 3
    I30_INDEX = 4
    # (v1, v2, v3, i20, i30)
    vals = dict()

    # Read in PSCAD .CSV data
    print('Opening PSCAD CSV file...')
    with open('pscad_test_dat.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        pscad_breaker_open_time = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                # print('Column names are: ' + ', '.join(row))
                line_count += 1
            else:
                pscad_t_vals.append(row[0])
                pscad_v1_vals.append(row[1])
                pscad_v2_vals.append(row[2])
                pscad_v3_vals.append(row[3])
                pscad_ishort_vals.append(row[4])
                line_count += 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' lines.')

    # Time delta used for calculations
    delta_t = 10e-6
    # Simulation time
    simulation_time = 200e-3

    # Propagation constants
    Lline = 1e-3 # mH/km
    Cline = 12e-9 # nF/km
    line_length = 100 # km
    # Velocity of propagation
    a = sqrt(1/(Lline*Cline)) # km/s
    # Transmission line history time delta
    tau = line_length/a # s
    print("a = " + str(a) + " km/s")
    print("tau = " + str(tau) + " s")

    # Zero initial conditions
    v1 = 0.0
    v2 = 0.0
    v2_tmintau = 0.0
    v3 = 0.0
    v3_tmintau = 0.0
    i20 = 0.0
    i20_tmintau = 0.0
    i30 = 0.0
    i30_tmintau = 0.0
    t_vals = []
    vals = []
    t_vals.append(0)
    vals.append((230, v2, v3, i20, i30))

    # Start at t = delta_t
    t = delta_t

    # Run simulation
    while t < simulation_time:

        # Iterpolate required values at t - tau
        # We are guaranteed a t entry for every delta_t
        # Make sure we have appropriate data to extrapolate from
        if len(t_vals) > 1 and (t - tau) > 0:
            # Find the smallest value that t - tau is above
            t0 = int(round((t - tau)/delta_t)) * delta_t
            # Next value to interpolate with is at t0 + delta_t
            t1 = t0 + delta_t
            # print("t: %g" % t)
            # print("t0: %g" % t0)
            # print("t1: %g" % t1)
            # print("t-tau: %g\n" % (t-tau))
            # With upper and lower values, interpolate a new one in between for t - tau
            t0_index = int(round(t0/delta_t))
            t1_index = int(round(t1/delta_t))
            v2_tmintau = interpolate(t - tau, t0, vals[t0_index][V2_INDEX], t1, vals[t1_index][V2_INDEX])
            v3_tmintau = interpolate(t - tau, t0, vals[t0_index][V3_INDEX], t1, vals[t1_index][V3_INDEX])
            i20_tmintau = interpolate(t - tau, t0, vals[t0_index][I20_INDEX], t1, vals[t1_index][I20_INDEX])
            i30_tmintau = interpolate(t - tau, t0, vals[t0_index][I30_INDEX], t1, vals[t1_index][I30_INDEX])
        else:
            # Can't extrapolate, use 0
            t0 = 0
            t1 = 0
            t0_index = int(round(t0/delta_t))
            t1_index = int(round(t1/delta_t))
            v2_tmintau = 0
            v3_tmintau = 0
            i20_tmintau = 0
            i30_tmintau = 0

        # Voltage calculations
        # These voltage calculations are generated in MATLAB
        # by v_next[] = inv(GAA)*hA - inv(GAA)*GAB*vs
        v2_next = 227.73157000187032039661971793976*cos(377.0*t) - 2.6247319811815387817513823036834*i30_tmintau - 0.0003930315968080709353401758969589*v2_tmintau - 0.0094697075254947917026040267266773*v3_tmintau - 0.10893711331418238288249439638031*i20_tmintau
        v3_next = 0
        # Next currents
        i20_next = 0.0033203272286915516454473282287533*v2_next - 0.88362562945898092608117343658749*i30_tmintau - 0.036674077968279601044949287851997*i20_tmintau - 0.00013231552578197494803240204626591*v2_tmintau - 0.003188011702909576512471947744615*v3_tmintau
        i30_next = - 0.88362562945898092608117343658749*i20_tmintau - 0.036674077968279601044949287851997*i30_tmintau - 0.003188011702909576512471947744615*v2_tmintau - 0.00013231552578197494803240204626591*v3_tmintau

        # Next iteration
        i20 = i20_next  
        i30 = i30_next 
        v2 = v2_next
        v3 = v3_next

        # Append the values for the current iteration
        # (v1, v2, v3, v34_vals, v4, v5, i23, i30, i40, i50)
        t_vals.append(t)
        vals.append((230*cos(377*t), v2, v3, i20, i30))

        # Go to next time step
        t = t + delta_t

    # Transpose data
    transposed_vals = zip(*vals)

    # Plots for publication

    # Single
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V1_INDEX], label='$Simulation\\ V1$')
    ax.plot(pscad_t_vals, pscad_v1_vals, linestyle='dashed', label='$PSCAD\\ V_1$')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)')
    ax.legend(loc='best', fancybox=True, shadow=True)
    ax.grid()
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('test_v1.png')

    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V2_INDEX], label='$Simulation\\ V2$')
    ax.plot(pscad_t_vals, pscad_v2_vals, linestyle='dashed', label='$PSCAD\\ V_2$')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)')
    ax.legend(loc='best', fancybox=True, shadow=True)
    ax.grid()
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('test_v2.png')

    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V3_INDEX], label='$Simulation\\ V3$')
    ax.plot(pscad_t_vals, pscad_v3_vals, linestyle='dashed', label='$PSCAD\\ V_3$')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)')
    ax.legend(loc='best', fancybox=True, shadow=True)
    ax.grid()
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('test_v3.png')

    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[I30_INDEX], label='$Simulation\\ I_{30}$')
    ax.plot(pscad_t_vals, pscad_ishort_vals, linestyle='dashed', label='$PSCAD\\ I_{30}$')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='Simulation vs. PSCAD Short Circuit Current')
    ax.legend(loc='best', fancybox=True, shadow=True)
    ax.grid()
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('test_ishort.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 3 solution generator')

    args = parser.parse_args()

    main(args)