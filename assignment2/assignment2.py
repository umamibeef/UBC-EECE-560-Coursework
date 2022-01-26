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
    "font.size": 10 ,
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False     # don't setup fonts from rc parameters
})

# Main function
def main(args):
    # Lists for results
    pscad_t_vals = []
    pscad_v1_vals = []
    pscad_v2_vals = []
    pscad_v3_vals = []
    pscad_v34_vals = []
    pscad_v4_vals = []
    pscad_v5_vals = []
    pscad_ishort_vals = []

    t_vals = []
    v1_vals = []
    v2_vals = []
    v3_vals = []
    v34_vals = []
    v4_vals = []
    v5_vals = []

    i23_vals = []
    i30_vals = []
    i40_vals = []
    i50_vals = []

    # Read in PSCAD .CSV data
    print('Opening PSCAD CSV file...')
    with open('PSCAD_data.csv') as csv_file:
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
                pscad_v34_vals.append(row[4])
                pscad_v4_vals.append(row[5])
                pscad_v5_vals.append(row[6])
                pscad_ishort_vals.append(row[7])
                # Check if break opned
                if float(row[8]) > 0 and pscad_breaker_open_time == 0:
                    pscad_breaker_open_time = float(row[0])
                    print('PSCAD breaker opened at %g' % float(pscad_breaker_open_time))
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

    # Component values
    R1 = 3
    L1 = 0.35
    C1 = 10e-9
    R2 = 50
    L2 = 0.1
    C2 = 600e-9

    # Start at t = delta_t
    t = delta_t
    # Zero initial conditions
    v2 = 0.0
    v3 = 0.0
    v4 = 0.0
    v5 = 0.0
    i23 = 0.0
    i30 = 0.0
    i40 = 0.0
    i50 = 0.0

    # Run simulation
    while t < simulation_time:

        # Append the values for the current iteration
        t_vals.append(t)
        v1_vals.append(230*math.cos(377*t))
        v2_vals.append(v2)
        v3_vals.append(v3)
        v34_vals.append(v3-v4)
        v4_vals.append(v4)
        v5_vals.append(v5)
        i23_vals.append(i23)
        i30_vals.append(i30)
        i40_vals.append(i40)
        i50_vals.append(i50)

        # Voltage calculations
        # These voltage calculations are generated in MATLAB
        # by v_next[] = inv(GAA)*hA - inv(GAA)*GAB*vs
        if breaker_open_flag == False:
            # Closed switch solution (breaker closed)
            v2_next = (230.0*math.cos(377.0*t)*(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 4.0*C1*L1*R2*delta_t + 4.0*C2*L1*R2*delta_t))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) - (1.0*R1*delta_t**3*(v5 + (2.0*L2*i50)/delta_t))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) + (R1*delta_t**2*(2.0*L2 + R2*delta_t)*((2.0*C1*(v3 + (0.5*delta_t*i30)/C1))/delta_t + (2.0*C2*(v4 + (0.5*delta_t*i40)/C2))/delta_t + (0.5*delta_t*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/L1))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) - (0.5*R1*delta_t*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t)*(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 4.0*C1*L1*R2*delta_t + 4.0*C2*L1*R2*delta_t))/(L1*(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t))
            v3_next = (230.0*delta_t**2*math.cos(377.0*t)*(2.0*L2 + R2*delta_t))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) - (1.0*delta_t**2*(2.0*L1 + R1*delta_t)*(v5 + (2.0*L2*i50)/delta_t))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) + (delta_t*(2.0*L1 + R1*delta_t)*(2.0*L2 + R2*delta_t)*((2.0*C1*(v3 + (0.5*delta_t*i30)/C1))/delta_t + (2.0*C2*(v4 + (0.5*delta_t*i40)/C2))/delta_t + (0.5*delta_t*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/L1))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) - (0.5*R1*delta_t**3*(2.0*L2 + R2*delta_t)*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/(L1*(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t))
            v4_next = v3_next
            v5_next = (460.0*L2*delta_t**2*math.cos(377.0*t))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) - (1.0*delta_t*(v5 + (2.0*L2*i50)/delta_t)*(2.0*L1*delta_t + R1*delta_t**2 + R2*delta_t**2 + 4.0*C1*L1*R2 + 4.0*C2*L1*R2 + 2.0*C1*R1*R2*delta_t + 2.0*C2*R1*R2*delta_t))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) + (2.0*L2*delta_t*(2.0*L1 + R1*delta_t)*((2.0*C1*(v3 + (0.5*delta_t*i30)/C1))/delta_t + (2.0*C2*(v4 + (0.5*delta_t*i40)/C2))/delta_t + (0.5*delta_t*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/L1))/(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t) - (1.0*L2*R1*delta_t**3*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/(L1*(2.0*L1*delta_t**2 + 2.0*L2*delta_t**2 + R1*delta_t**3 + R2*delta_t**3 + 8.0*C1*L1*L2 + 8.0*C2*L1*L2 + 2.0*C1*R1*R2*delta_t**2 + 2.0*C2*R1*R2*delta_t**2 + 4.0*C1*L1*R2*delta_t + 4.0*C1*L2*R1*delta_t + 4.0*C2*L1*R2*delta_t + 4.0*C2*L2*R1*delta_t))
        else:
            # Open switch solution (breaker opened)
            v2_next = (230.0*math.cos(377.0*t)*(delta_t**2 + 4.0*C1*L1))/(delta_t**2 + 2.0*C1*R1*delta_t + 4.0*C1*L1) + (R1*delta_t**2*((2.0*C1*(v3 + (0.5*delta_t*i30)/C1))/delta_t + (0.5*delta_t*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/L1))/(delta_t**2 + 2.0*C1*R1*delta_t + 4.0*C1*L1) - (0.5*R1*delta_t*(delta_t**2 + 4.0*C1*L1)*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/(L1*(delta_t**2 + 2.0*C1*R1*delta_t + 4.0*C1*L1))
            v3_next = (230.0*delta_t**2*math.cos(377.0*t))/(delta_t**2 + 2.0*C1*R1*delta_t + 4.0*C1*L1) + (delta_t*(2.0*L1 + R1*delta_t)*((2.0*C1*(v3 + (0.5*delta_t*i30)/C1))/delta_t + (0.5*delta_t*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/L1))/(delta_t**2 + 2.0*C1*R1*delta_t + 4.0*C1*L1) - (0.5*R1*delta_t**3*(v2 - 1.0*v3 + (2.0*L1*i23)/delta_t))/(L1*(delta_t**2 + 2.0*C1*R1*delta_t + 4.0*C1*L1))
            v4_next = (2.0*C2*(2.0*L2 + R2*delta_t)*(v4 + (0.5*delta_t*i40)/C2))/(delta_t**2 + 2.0*C2*R2*delta_t + 4.0*C2*L2) - (1.0*delta_t**2*(v5 + (2.0*L2*i50)/delta_t))/(delta_t**2 + 2.0*C2*R2*delta_t + 4.0*C2*L2)
            v5_next = (4.0*C2*L2*(v4 + (0.5*delta_t*i40)/C2))/(delta_t**2 + 2.0*C2*R2*delta_t + 4.0*C2*L2) - (1.0*delta_t*(v5 + (2.0*L2*i50)/delta_t)*(delta_t + 2.0*C2*R2))/(delta_t**2 + 2.0*C2*R2*delta_t + 4.0*C2*L2)
        
        # Next currents
        i23_next = (delta_t*(v2_next - v3_next))/(2*L1) + (delta_t*(v2 - v3 + (2*L1*i23)/delta_t))/(2*L1)
        i30_next = (2*C1*v3_next)/delta_t - (2*C1*(v3 + (delta_t*i30)/(2*C1)))/delta_t
        i40_next = (2*C2*v4_next)/delta_t - (2*C2*(v4 + (delta_t*i40)/(2*C2)))/delta_t
        i50_next = (delta_t*(v5 + (2*L2*i50)/delta_t))/(2*L2) + (delta_t*v5_next)/(2*L2)
        
        # Check if breaker has opened
        if (t >= breaker_time) and (not breaker_open_flag):
            # Check for zero crossing
            i_left_prev = i23 - i30
            i_left_next = i23_next - i30_next
            if ((i_left_next * i_left_prev) < 0):
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
        v2 = v2_next
        v3 = v3_next
        v4 = v4_next
        v5 = v5_next

        # Go to next time step
        t = t + delta_t

    # Plots for publication

    legend_font_size = 6
    # Superimposed
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, v1_vals, label='$Simulation\\ V_1$')
    ax.plot(t_vals, v3_vals, label='$Simulation\\ V_2$')
    ax.plot(t_vals, v4_vals, label='$Simulation\\ V_3$')
    ax.plot(pscad_t_vals, pscad_v1_vals, linestyle='dashed', label='$PSCAD\\ V_1$')
    ax.plot(pscad_t_vals, pscad_v3_vals, linestyle='dashed', label='$PSCAD\\ V_2$')
    ax.plot(pscad_t_vals, pscad_v4_vals, linestyle='dashed', label='$PSCAD\\ V_3$')
    ax.axvline(x=breaker_open_time, linestyle='dashed', c='r')
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

    # 1. Nodal Voltage Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, v1_vals, label='$Simulation\\ V_1$')
    ax.plot(t_vals, v3_vals, label='$Simulation\\ V_2$')
    ax.plot(t_vals, v4_vals, label='$Simulation\\ V_3$')
    ax.plot(pscad_t_vals, pscad_v1_vals, linestyle='dashed', label='$PSCAD\\ V_1$')
    ax.plot(pscad_t_vals, pscad_v3_vals, linestyle='dashed', label='$PSCAD\\ V_2$')
    ax.plot(pscad_t_vals, pscad_v4_vals, linestyle='dashed', label='$PSCAD\\ V_3$')
    ax.axvline(x=breaker_open_time, linestyle='dashed', c='r')
    ax.set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation vs. PSCAD Nodal Voltages')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-400,400)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('nodal_voltage_comparison_plots_zoom.pgf')
    fig.savefig('nodal_voltage_comparison_plots_zoom.png')

    # 2. Short Circuit Current Comparison Plots
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, i50_vals, label='$Simulation\\ I_{40}$')
    ax.plot(pscad_t_vals, pscad_ishort_vals, linestyle='dashed', label='$PSCAD\\ I_{40}$')
    ax.axvline(x=breaker_open_time, linestyle='dashed', c='r')
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

    # 2. Short Circuit Current Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, i50_vals, label='$Simulation\\ I_{40}$')
    ax.plot(pscad_t_vals, pscad_ishort_vals, linestyle='dashed', label='$PSCAD\\ I_{40}$')
    ax.axvline(x=breaker_open_time, linestyle='dashed', c='r')
    ax.set(xlabel='time (s)', ylabel='current (kA)', title='Simulation vs. PSCAD Short Circuit Current')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax.grid()
    ax.set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    ax.set_ylim(-1.5,1.5)
    ax.set_xlim(0.06, 0.09)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('short_circuit_current_comparison_plots_zoom.pgf')
    fig.savefig('short_circuit_current_comparison_plots_zoom.png')

    # 3. Breaker Voltage Comparison Plots
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, v34_vals, label='$Simulation\\ V_{23}$')
    ax.plot(pscad_t_vals, pscad_v34_vals, linestyle='dashed', label='$PSCAD\\ V_{23}$')
    ax.axvline(x=breaker_open_time, linestyle='dashed', c='r')
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

    # 3. Breaker Voltage Comparison Plots (Zoomed)
    # Simulation & PSCAD
    fig, ax = plt.subplots(1)
    ax.plot(t_vals, v34_vals, label='$Simulation\\ V_{23}$')
    ax.plot(pscad_t_vals, pscad_v34_vals, linestyle='dashed', label='$PSCAD\\ V_{23}$')
    ax.axvline(x=breaker_open_time, linestyle='dashed', c='r')
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

    # Side-by-side

    # # 1. Nodal Voltage Comparison Plots
    # # Simulation & PSCAD
    # fig, ax = plt.subplots(2)
    # ax[0].plot(t_vals, v1_vals, label='$Simulation\\ V_1$')
    # ax[0].plot(t_vals, v3_vals, label='$Simulation\\ V_2$')
    # ax[0].plot(t_vals, v4_vals, label='$Simulation\\ V_3$')
    # ax[0].axvline(x=breaker_open_time, linestyle='dashed', c='r')
    # ax[0].set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation Nodal Voltages')
    # ax[0].legend(loc='best', fancybox=True, shadow=True)
    # ax[0].grid()
    # ax[0].set_xticks(list(ax[0].get_xticks()) + [breaker_open_time])
    # ax[0].set_ylim(-400,400)
    # ax[0].set_xlim(0, 0.2)
    # # PSCAD
    # ax[1].plot(pscad_t_vals, pscad_v1_vals, label='$PSCAD\\ V_1$')
    # ax[1].plot(pscad_t_vals, pscad_v3_vals, label='$PSCAD\\ V_2$')
    # ax[1].plot(pscad_t_vals, pscad_v4_vals, label='$PSCAD\\ V_3$')
    # ax[1].axvline(x=pscad_breaker_open_time, linestyle='dashed', c='r')
    # ax[1].set(xlabel='time (s)', ylabel='voltage (kV)', title='PSCAD Nodal Voltages')
    # ax[1].legend(loc='best', fancybox=True, shadow=True)
    # ax[1].grid()
    # ax[1].set_xticks(list(ax[1].get_xticks()) + [pscad_breaker_open_time])
    # ax[1].set_ylim(-400,400)
    # ax[1].set_xlim(0, 0.2)
    # fig.set_size_inches(6.5,8)
    # fig.tight_layout()
    # plt.show()
    # fig.savefig('nodal_voltage_comparison_plots.pgf')
    # fig.savefig('nodal_voltage_comparison_plots.png')

    # # 1. Nodal Voltage Comparison Plots (Zoomed)
    # # Simulation & PSCAD
    # fig, ax = plt.subplots(2)
    # ax[0].plot(t_vals, v1_vals, label='$Simulation\\ V_1$')
    # ax[0].plot(t_vals, v3_vals, label='$Simulation\\ V_2$')
    # ax[0].plot(t_vals, v4_vals, label='$Simulation\\ V_3$')
    # ax[0].axvline(x=breaker_open_time, linestyle='dashed', c='r')
    # ax[0].set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation Nodal Voltages')
    # ax[0].legend(loc='best', fancybox=True, shadow=True)
    # ax[0].grid()
    # ax[0].set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    # ax[0].set_ylim(-400,400)
    # ax[0].set_xlim(0.06, 0.09)
    # # PSCAD
    # ax[1].plot(pscad_t_vals, pscad_v1_vals, label='$PSCAD\\ V_1$')
    # ax[1].plot(pscad_t_vals, pscad_v3_vals, label='$PSCAD\\ V_2$')
    # ax[1].plot(pscad_t_vals, pscad_v4_vals, label='$PSCAD\\ V_3$')
    # ax[1].axvline(x=pscad_breaker_open_time, linestyle='dashed', c='r')
    # ax[1].set(xlabel='time (s)', ylabel='voltage (kV)', title='PSCAD Nodal Voltages')
    # ax[1].legend(loc='best', fancybox=True, shadow=True)
    # ax[1].grid()
    # ax[1].set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    # ax[1].set_ylim(-400,400)
    # ax[1].set_xlim(0.06, 0.09)
    # fig.set_size_inches(6.5,8)
    # fig.tight_layout()
    # plt.show()
    # fig.savefig('nodal_voltage_comparison_plots_zoom.pgf')
    # fig.savefig('nodal_voltage_comparison_plots_zoom.png')

    # # 2. Short Circuit Current Comparison Plots
    # # Simulation & PSCAD
    # fig, ax = plt.subplots(2)
    # ax[0].plot(t_vals, i50_vals, label='$Simulation\\ I_{40}$')
    # ax[0].axvline(x=breaker_open_time, linestyle='dashed', c='r')
    # ax[0].set(xlabel='time (s)', ylabel='current (kA)', title='Simulation Short Circuit Current')
    # ax[0].legend(loc='best', fancybox=True, shadow=True)
    # ax[0].grid()
    # ax[0].set_xticks(list(ax[0].get_xticks()) + [breaker_open_time])
    # ax[0].set_ylim(-1.5,1.5)
    # ax[0].set_xlim(0, 0.2)
    # # PSCAD
    # ax[1].plot(pscad_t_vals, pscad_ishort_vals, label='$PSCAD\\ I_{40}$')
    # ax[1].axvline(x=pscad_breaker_open_time, linestyle='dashed', c='r')
    # ax[1].set(xlabel='time (s)', ylabel='current (kA)', title='PSCAD Short Circuit Current')
    # ax[1].legend(loc='best', fancybox=True, shadow=True)
    # ax[1].grid()
    # ax[1].set_xticks(list(ax[1].get_xticks()) + [pscad_breaker_open_time])
    # ax[1].set_ylim(-1.5,1.5)
    # ax[1].set_xlim(0, 0.2)
    # fig.set_size_inches(6.5,8)
    # fig.tight_layout()
    # plt.show()
    # fig.savefig('short_circuit_current_comparison_plots.pgf')
    # fig.savefig('short_circuit_current_comparison_plots.png')

    # # 2. Short Circuit Current Comparison Plots (Zoomed)
    # # Simulation & PSCAD
    # fig, ax = plt.subplots(2)
    # ax[0].plot(t_vals, i50_vals, label='$Simulation\\ I_{40}$')
    # ax[0].axvline(x=breaker_open_time, linestyle='dashed', c='r')
    # ax[0].set(xlabel='time (s)', ylabel='current (kA)', title='Simulation Short Circuit Current')
    # ax[0].legend(loc='best', fancybox=True, shadow=True)
    # ax[0].grid()
    # ax[0].set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    # ax[0].set_ylim(-1.5,1.5)
    # ax[0].set_xlim(0.06, 0.09)
    # # PSCAD
    # ax[1].plot(pscad_t_vals, pscad_ishort_vals, label='$PSCAD\\ I_{40}$')
    # ax[1].axvline(x=pscad_breaker_open_time, linestyle='dashed', c='r')
    # ax[1].set(xlabel='time (s)', ylabel='current (kA)', title='PSCAD Short Circuit Current')
    # ax[1].legend(loc='best', fancybox=True, shadow=True)
    # ax[1].grid()
    # ax[1].set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    # ax[1].set_ylim(-1.5,1.5)
    # ax[1].set_xlim(0.06, 0.09)
    # fig.set_size_inches(6.5,8)
    # fig.tight_layout()
    # plt.show()
    # fig.savefig('short_circuit_current_comparison_plots_zoom.pgf')
    # fig.savefig('short_circuit_current_comparison_plots_zoom.png')

    # # 3. Breaker Voltage Comparison Plots
    # # Simulation & PSCAD
    # fig, ax = plt.subplots(2)
    # ax[0].plot(t_vals, v34_vals, label='$Simulation\\ V_{23}$')
    # ax[0].axvline(x=breaker_open_time, linestyle='dashed', c='r')
    # ax[0].set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation Breaker Voltage')
    # ax[0].legend(loc='best', fancybox=True, shadow=True)
    # ax[0].grid()
    # ax[0].set_xticks(list(ax[0].get_xticks()) + [breaker_open_time])
    # ax[0].set_ylim(-400,400)
    # ax[0].set_xlim(0, 0.2)
    # # PSCAD
    # ax[1].plot(pscad_t_vals, pscad_v34_vals, label='$PSCAD\\ V_{23}$')
    # ax[1].axvline(x=pscad_breaker_open_time, linestyle='dashed', c='r')
    # ax[1].set(xlabel='time (s)', ylabel='voltage (kV)', title='PSCAD Breaker Voltage')
    # ax[1].legend(loc='best', fancybox=True, shadow=True)
    # ax[1].grid()
    # ax[1].set_xticks(list(ax[1].get_xticks()) + [pscad_breaker_open_time])
    # ax[1].set_ylim(-400,400)
    # ax[1].set_xlim(0, 0.2)
    # fig.set_size_inches(6.5,8)
    # fig.tight_layout()
    # plt.show()
    # fig.savefig('breaker_voltage_comparison_plots.pgf')
    # fig.savefig('breaker_voltage_comparison_plots.png')

    # # 3. Breaker Voltage Comparison Plots (Zoomed)
    # # Simulation & PSCAD
    # fig, ax = plt.subplots(2)
    # ax[0].plot(t_vals, v34_vals, label='$Simulation\\ V_{23}$')
    # ax[0].axvline(x=breaker_open_time, linestyle='dashed', c='r')
    # ax[0].set(xlabel='time (s)', ylabel='voltage (kV)', title='Simulation Breaker Voltage')
    # ax[0].legend(loc='best', fancybox=True, shadow=True)
    # ax[0].grid()
    # ax[0].set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    # ax[0].set_ylim(-400,400)
    # ax[0].set_xlim(0.06, 0.09)
    # # PSCAD
    # ax[1].plot(pscad_t_vals, pscad_v34_vals, label='$PSCAD\\ V_{23}$')
    # ax[1].axvline(x=pscad_breaker_open_time, linestyle='dashed', c='r')
    # ax[1].set(xlabel='time (s)', ylabel='voltage (kV)', title='PSCAD Breaker Voltage')
    # ax[1].legend(loc='best', fancybox=True, shadow=True)
    # ax[1].grid()
    # ax[1].set_xticks(np.arange(0.06, 0.09, step=0.005).tolist())
    # ax[1].set_ylim(-400,400)
    # ax[1].set_xlim(0.06, 0.09)
    # fig.set_size_inches(6.5,8)
    # fig.tight_layout()
    # plt.show()
    # fig.savefig('breaker_voltage_comparison_plots_zoom.pgf')
    # fig.savefig('breaker_voltage_comparison_plots_zoom.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 2 solution generator')

    args = parser.parse_args()

    main(args)