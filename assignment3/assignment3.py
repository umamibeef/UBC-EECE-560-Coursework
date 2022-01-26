import argparse
import csv
import math
import matplotlib
import numpy as np
 
 # Matplotlib export settings
matplotlib.use("pgf")
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.size": 10,
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False     # don't setup fonts from rc parameters
})

# Interpolation function
# t0 < t < t1
def interpolate(t, t0, y0, t1, y1):
    # Linear interpolation between points (t0, y0) and (t1, y1)
    return y0 + (t - t0)*(y1 - y0)/(t1 - t0)
    # return y1

# Main function
def main(args):
    # Lists for results
    pscad_assignment2_t_vals = []
    pscad_assignment2_v1_vals = []
    pscad_assignment2_v2_vals = []
    pscad_assignment2_v3_vals = []
    pscad_assignment2_v34_vals = []
    pscad_assignment2_v4_vals = []
    pscad_assignment2_v5_vals = []
    pscad_assignment2_ishort_vals = []

    pscad_assignment3_t_vals = []
    pscad_assignment3_v1_vals = []
    pscad_assignment3_v2_vals = []
    pscad_assignment3_v3_vals = []
    pscad_assignment3_v34_vals = []
    pscad_assignment3_v4_vals = []
    pscad_assignment3_v5_vals = []
    pscad_assignment3_ishort_vals = []

    # tuple indices
    V1_INDEX = 0
    V2_INDEX = 1
    V3_INDEX = 2
    V34_INDEX = 3
    V4_INDEX = 4
    V5_INDEX = 5
    I23_INDEX = 6
    I30_INDEX = 7
    I40_INDEX = 8
    I50_INDEX = 9
    ISHORT_INDEX = 10
    # (v1, v2, v3, v34_vals, v4, v5, i23, i30, i40, i50, ishort)
    vals = dict()

    # Read in PSCAD .CSV data
    print('*** Opening Assignment 2 PSCAD CSV file...')
    with open('pscad_assignment2_dat.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        pscad_assignment2_breaker_open_time = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                # print('Column names are: ' + ', '.join(row))
                line_count += 1
            else:
                pscad_assignment2_t_vals.append(row[0])
                pscad_assignment2_v1_vals.append(row[1])
                pscad_assignment2_v2_vals.append(row[2])
                pscad_assignment2_v3_vals.append(row[3])
                pscad_assignment2_v34_vals.append(row[4])
                pscad_assignment2_v4_vals.append(row[5])
                pscad_assignment2_v5_vals.append(row[6])
                pscad_assignment2_ishort_vals.append(row[7])
                # Check if break opned
                if float(row[8]) > 0 and pscad_assignment2_breaker_open_time == 0:
                    pscad_assignment2_breaker_open_time = float(row[0])
                    print('PSCAD breaker opened at %g' % float(pscad_assignment2_breaker_open_time))
                line_count += 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' lines.')
    print('*** Opening Assignment 3 PSCAD CSV file...')
    with open('pscad_assignment3_dat.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        pscad_assignment3_breaker_open_time = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                # print('Column names are: ' + ', '.join(row))
                line_count += 1
            else:
                pscad_assignment3_t_vals.append(row[0])
                pscad_assignment3_v1_vals.append(row[1])
                pscad_assignment3_v2_vals.append(row[2])
                pscad_assignment3_v3_vals.append(row[3])
                pscad_assignment3_v34_vals.append(row[4])
                pscad_assignment3_v4_vals.append(row[5])
                pscad_assignment3_ishort_vals.append(row[6])
                # Check if break opned
                if float(row[7]) > 0 and pscad_assignment3_breaker_open_time == 0:
                    pscad_assignment3_breaker_open_time = float(row[0])
                    print('PSCAD breaker opened at %g' % float(pscad_assignment3_breaker_open_time))
                line_count += 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' lines.')

    # Time delta used for calculations
    delta_t = 10e-6
    # Time for breaker switch to open
    breaker_time = 66.67e-3
    # Zero crossing flag
    zero_crossing_flag = False
    # Breaker open flag
    breaker_open_flag = False
    # Time breaker opened
    breaker_open_time = 0
    # Simulation time
    simulation_time = 200e-3

    # Propagation constants
    Lline = 1e-3 # mH/km
    Cline = 12e-9 # nF/km
    line_length = 100 # km
    # Velocity of propagation
    a = math.sqrt(1/(Lline*Cline)) # km/s
    # Transmission line history time delta
    tau = line_length/a # s
    print("a = " + str(a) + " km/s")
    print("tau = " + str(tau) + " s")

    # Zero initial conditions
    v1 = 0.0
    v2 = 0.0
    v3 = 0.0
    v34 = 0.0
    v4 = 0.0
    v4_tmintau = 0.0
    v5 = 0.0
    v5_tmintau = 0.0
    i23 = 0.0
    i30 = 0.0
    i40 = 0.0
    i40_tmintau = 0.0
    i50 = 0.0
    i50_tmintau = 0.0
    ishort = 0.0
    t_vals = []
    vals = []
    t_vals.append(0)
    vals.append((230, v2, v3, v3-v4, v4, v5, i23, i30, i40, i50, ishort))

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
            t0_index = int(round(t0/delta_t))
            t1_index = int(round(t1/delta_t))
            v4_tmintau = interpolate(t - tau, t0, vals[t0_index][V4_INDEX], t1, vals[t1_index][V4_INDEX])
            v5_tmintau = 0.0 # always 0
            i40_tmintau = interpolate(t - tau, t0, vals[t0_index][I40_INDEX], t1, vals[t1_index][I40_INDEX])
            i50_tmintau = interpolate(t - tau, t0, vals[t0_index][I50_INDEX], t1, vals[t1_index][I50_INDEX])
        else:
            # Can't extrapolate, use 0
            t0 = 0.0
            t1 = 0.0
            t0_index = int(round(t0/delta_t))
            t1_index = int(round(t1/delta_t))
            v4_tmintau = 0.0
            v5_tmintau = 0.0
            i40_tmintau = 0.0
            i50_tmintau = 0.0

        # Check if breaker is opened or closed and set breaker resistance accordingly
        if breaker_open_flag == False:
            # Breaker has very small resistance (closed)
            Rbrk = 1e-9
        else:
            # Break has very large resistance (open)
            Rbrk = 1e9

        # Voltage calculations
        # These voltage calculations are generated in MATLAB
        # by v_next[] = inv(GAA)*hA - inv(GAA)*GAB*vs
        v2_next = (1.1358149378044935437410460720169e+37*(0.036674077968279601044949287851997*i40_tmintau + 0.88362562945898092608117343658749*i50_tmintau + 0.00013231552578197494803240204626591*v4_tmintau + 0.003188011702909576512471947744615*v5_tmintau))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) - (59109745109237760000.0*(89959801610818431.0*Rbrk + 71754549774979273936.0)*(i23 + 0.000014285714285714285714285714285714*v2 - 0.000014285714285714285714285714285714*v3))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) + (19703248369745920000.0*(1914038332145073.0*Rbrk + 576460752303423488.0)*(i23 + i30 + 0.000014285714285714285714285714285714*v2 + 0.0019857142857142857142857142857143*v3))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) + (4531747125041561600000.0*math.cos(377.0*t)*(89959801610818431.0*Rbrk + 71754549774979273936.0))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39)
        v3_next = (2.650348436370931702638844437332e+41*(0.036674077968279601044949287851997*i40_tmintau + 0.88362562945898092608117343658749*i50_tmintau + 0.00013231552578197494803240204626591*v4_tmintau + 0.003188011702909576512471947744615*v5_tmintau))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) - (19703248369745920000.0*(1914038332145073.0*Rbrk + 576460752303423488.0)*(i23 + 0.000014285714285714285714285714285714*v2 - 0.000014285714285714285714285714285714*v3))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) + (1510582375013853866666.6666666667*math.cos(377.0*t)*(1914038332145073.0*Rbrk + 576460752303423488.0))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) + (459762165209107818750000.0*(1914038332145073.0*Rbrk + 576460752303423488.0)*(i23 + i30 + 0.000014285714285714285714285714285714*v2 + 0.0019857142857142857142857142857143*v3))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39)
        v4_next = (8.7079145231677838353480198854629e+38*math.cos(377.0*t))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) - (1.1358149378044935437410460720169e+37*(i23 + 0.000014285714285714285714285714285714*v2 - 0.000014285714285714285714285714285714*v3))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) + (2.650348436370931702638844437332e+41*(i23 + i30 + 0.000014285714285714285714285714285714*v2 + 0.0019857142857142857142857142857143*v3))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39) + (576460752303423488.0*(926092079874797609969.0*Rbrk + 459762165209107818750000.0)*(0.036674077968279601044949287851997*i40_tmintau + 0.88362562945898092608117343658749*i50_tmintau + 0.00013231552578197494803240204626591*v4_tmintau + 0.003188011702909576512471947744615*v5_tmintau))/(1.772575739976319342526327845379e+36*Rbrk + 1.4138581449471162223771073543884e+39)
        v5_next = 0.0 # v5 always zero!

        # Next currents
        i23_next = i23 + 0.000014285714285714285714285714285714*v2 - 0.000014285714285714285714285714285714*v3 + 0.000014285714285714285714285714285714*v2_next - 0.000014285714285714285714285714285714*v3_next
        i30_next = 0.002*v3_next - 0.002*v3 - 1.0*i30
        i40_next = 0.0033203272286915516454473282287533*v4_next - 0.88362562945898092608117343658749*i50_tmintau - 0.036674077968279601044949287851997*i40_tmintau - 0.00013231552578197494803240204626591*v4_tmintau - 0.003188011702909576512471947744615*v5_tmintau
        i50_next = - 0.88362562945898092608117343658749*i40_tmintau - 0.036674077968279601044949287851997*i50_tmintau - 0.003188011702909576512471947744615*v4_tmintau - 0.00013231552578197494803240204626591*v5_tmintau

        # Check if breaker has opened
        if (t >= breaker_time) and (not breaker_open_flag):
            # Check for zero crossing
            i_break_prev = (v4 - v3)/Rbrk
            i_break_next = (v4_next - v3_next)/Rbrk
            if ((i_break_next * i_break_prev) < 0):
                # Sign has changed, set our flag to true
                zero_crossing_flag = True
                # Breaker has now officially opened
                breaker_open_flag = True
                # Record time breaker has opened
                breaker_open_time = t
                print("Snap! Breaker opened @ t = %f" % breaker_open_time)

        # Next iteration
        i23 = i23_next  
        i30 = i30_next 
        i40 = i40_next 
        i50 = i50_next
        ishort = -i50_next # ishort's direction is opposite of i50
        v2 = v2_next
        v3 = v3_next
        v4 = v4_next
        v5 = v5_next

        # Append the values for the current iteration
        # (v1, v2, v3, v34_vals, v4, v5, i23, i30, i40, i50, ishort)
        t_vals.append(t)
        vals.append((230*math.cos(377*t), v2, v3, v3-v4, v4, v5, i23, i30, i40, i50, ishort))

        # Go to next time step
        t = t + delta_t

    # Transpose data
    transposed_vals = zip(*vals)

    # Plots for publication
    legend_font_size = 6

    # Superimposed
    # Nodal Voltage Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V1_INDEX], label='$Simulation\\ V_{one}$')
    ax.plot(t_vals, transposed_vals[V3_INDEX], label='$Simulation\\ V_{two}$')
    ax.plot(t_vals, transposed_vals[V4_INDEX], label='$Simulation\\ V_{three}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v1_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{one}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v3_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{two}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v4_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation vs. PSCAD Nodal Voltages')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(list(ax.get_xticks()) + [breaker_open_time])
    ax.set_ylim(-400,400)
    ax.set_xlim(0, 0.2)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('nodal_voltage_comparison_plots.pgf')
    fig.savefig('nodal_voltage_comparison_plots.png')

    # Superimposed
    # Nodal Voltage Comparison Plots (Zoomed)
    # PSCAD & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v1_vals, label='$PSCAD_3\\ V_{one}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v3_vals, label='$PSCAD_3\\ V_{two}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v4_vals, label='$PSCAD_3\\ V_{three}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_v1_vals, linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ V_{one}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_v3_vals, linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ V_{two}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_v4_vals, linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ V_{three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='$PSCAD_2$ vs. $PSCAD_3$ Nodal Voltages')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(list(ax.get_xticks()) + [breaker_open_time])
    ax.set_ylim(-400,400)
    ax.set_xlim(0, 0.2)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('pscad_nodal_voltage_comparison_plots.pgf')
    fig.savefig('pscad_nodal_voltage_comparison_plots.png')

    # Nodal Voltage Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V1_INDEX], label='$Simulation\\ V_{one}$')
    ax.plot(t_vals, transposed_vals[V3_INDEX], label='$Simulation\\ V_{two}$')
    ax.plot(t_vals, transposed_vals[V4_INDEX], label='$Simulation\\ V_{three}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v1_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{one}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v3_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{two}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v4_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation vs. PSCAD Nodal Voltages')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-400,400)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('nodal_voltage_comparison_plots_zoom_1.pgf')
    fig.savefig('nodal_voltage_comparison_plots_zoom_1.png')

    # Nodal Voltage Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V1_INDEX], label='$Simulation\\ V_{one}$')
    ax.plot(t_vals, transposed_vals[V3_INDEX], label='$Simulation\\ V_{two}$')
    ax.plot(t_vals, transposed_vals[V4_INDEX], label='$Simulation\\ V_{three}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v1_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{one}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v3_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{two}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v4_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation vs. PSCAD Nodal Voltages')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.0, 0.01, step=0.001).tolist())
    ax.set_ylim(-400,400)
    ax.set_xlim(0.0, 0.01)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('nodal_voltage_comparison_plots_zoom_2.pgf')
    fig.savefig('nodal_voltage_comparison_plots_zoom_2.png')

    # Short Circuit Current Comparison Plots
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[ISHORT_INDEX], label='$Simulation\\ I_{four, ground}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_ishort_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ I_{four, ground}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='Simulation vs. PSCAD Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(list(ax.get_xticks()) + [breaker_open_time])
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0, 0.2)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('short_circuit_current_comparison_plots.pgf')
    fig.savefig('short_circuit_current_comparison_plots.png')

    # Short Circuit Current Comparison Plots
    # PSCAD & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_ishort_vals, linestyle='solid', label='$PSCAD_3\\ I_{four, ground}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_ishort_vals, c='r', linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ I_{four, ground}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='$PSCAD_2$ vs. $PSCAD_3$ Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(list(ax.get_xticks()) + [breaker_open_time])
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0, 0.2)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('pscad_short_circuit_current_comparison_plots.pgf')
    fig.savefig('pscad_short_circuit_current_comparison_plots.png')

    # Short Circuit Current Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[ISHORT_INDEX], label='$Simulation\\ I_{four, ground}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_ishort_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ I_{four, ground}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='Simulation vs. PSCAD Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('short_circuit_current_comparison_plots_zoom_1.pgf')
    fig.savefig('short_circuit_current_comparison_plots_zoom_1.png')

    # Short Circuit Current Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_ishort_vals, linestyle='solid', label='$PSCAD_3\\ I_{four, ground}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_ishort_vals, c='r', linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ I_{four, ground}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='$PSCAD_2$ vs. $PSCAD_3$ Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('pscad_short_circuit_current_comparison_plots_zoom_1.pgf')
    fig.savefig('pscad_short_circuit_current_comparison_plots_zoom_1.png')

    # Short Circuit Current Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_ishort_vals, linestyle='solid', label='$PSCAD_3\\ I_{four, ground}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_ishort_vals, c='r', linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ I_{four, ground}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='$PSCAD_2$ vs. $PSCAD_3$ Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.0, 0.01, step=0.001).tolist())
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0.0, 0.01)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('pscad_short_circuit_current_comparison_plots_zoom_2.pgf')
    fig.savefig('pscad_short_circuit_current_comparison_plots_zoom_2.png')

    # Short Circuit Current Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[ISHORT_INDEX], label='$Simulation\\ I_{four, ground}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_ishort_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ I_{four, ground}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='Simulation vs. PSCAD Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.0, 0.01, step=0.001).tolist())
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0.0, 0.01)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('short_circuit_current_comparison_plots_zoom_2.pgf')
    fig.savefig('short_circuit_current_comparison_plots_zoom_2.png')

    # Breaker Voltage Comparison Plots
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V34_INDEX], label='$Simulation\\ V_{two, three}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v34_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{two, three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation vs. PSCAD Breaker Voltage')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(list(ax.get_xticks()) + [breaker_open_time])
    ax.set_ylim(-400,400)
    ax.set_xlim(0, 0.2)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('breaker_voltage_comparison_plots.pgf')
    fig.savefig('breaker_voltage_comparison_plots.png')

    # Breaker Voltage Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, transposed_vals[V34_INDEX], label='$Simulation\\ V_{two, three}$')
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v34_vals, linestyle='--', dashes=(5, 1), label='$PSCAD\\ V_{two, three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation vs. PSCAD Breaker Voltage')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-400,400)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('breaker_voltage_comparison_plots_zoom.pgf')
    fig.savefig('breaker_voltage_comparison_plots_zoom.png')

    # Breaker Voltage Comparison Plots
    # PSCAD & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v34_vals, linestyle='solid', label='$PSCAD_3\\ V_{two, three}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_v34_vals, linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ V_{two, three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='$PSCAD_2$ vs. $PSCAD_3$ Breaker Voltage')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(list(ax.get_xticks()) + [breaker_open_time])
    ax.set_ylim(-400,400)
    ax.set_xlim(0, 0.2)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('pscad_breaker_voltage_comparison_plots.pgf')
    fig.savefig('pscad_breaker_voltage_comparison_plots.png')

    # Breaker Voltage Comparison Plots (Zoomed)
    # PSCAD & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(pscad_assignment3_t_vals, pscad_assignment3_v34_vals, linestyle='solid', label='$PSCAD_3\\ V_{two, three}$')
    ax.plot(pscad_assignment2_t_vals, pscad_assignment2_v34_vals, linestyle='--', dashes=(5, 1), label='$PSCAD_2\\ V_{two, three}$')
    ax.axvline(x=breaker_open_time, linestyle='solid', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='$PSCAD_2$ vs. $PSCAD_3$ Breaker Voltage')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-400,400)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('pscad_breaker_voltage_comparison_plots_zoom.pgf')
    fig.savefig('pscad_breaker_voltage_comparison_plots_zoom.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 3 solution generator')

    args = parser.parse_args()

    main(args)