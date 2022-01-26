import argparse
import csv
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy

# Matplotlib export settings
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.size": 10 ,
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False     # don't setup fonts from rc parameters
})

# Constants

# indices for the results tuple
TIME_INDEX = 0
CURRENT_INDEX = 1
VOLTAGE_INDEX = 2

inductance = 0.02 # 20 mH
resistance = 10 # 10 Ohms
voltage = 10 # 10 Volts
simulation_time = 0.01 # 10 ms


GET_FLUX_KNEE = 0
GET_FLUX_INDUCTANCE = 1
GET_FLUX = 2
# Inductor flux anchor and inductance lookup table
def get_flux(current_input):

    # Normal region
    if (current_input >= -0.1) and (current_input <= 0.1):
        vals = (0.0, 0.02, (0.0 + 0.02*current_input))
    # Saturation region
    elif (current_input > 0.1) or (current_input < -0.1):
        vals = (1.956e-3, 4.444e-4, (1.956e-3 + 4.444e-4*current_input))
    else:
        raise ValueError('Invalid current')

    return vals

# Dataset Class
class Dataset:

    def __init__(self, description, linestyle):
        self.description = description
        self.data = [[],[],[]]
        self.linestyle = linestyle

# Plot generation
def generate_plots(name, datasets, xlims, delta_t):

    fig, ax = plt.subplots(2)
    for index, axis in enumerate(ax):
        # current
        if index == 0:
            for dataset in datasets:
                axis.plot(dataset.data[TIME_INDEX], dataset.data[CURRENT_INDEX], linestyle=dataset.linestyle, label=dataset.description)
            axis.set(xlabel='time (s)', ylabel=r'current (A)', title='Evaluations of $i_L(t)$ at $\Delta{}t=%gs$' % delta_t)
        # voltage
        elif index == 1:
            for dataset in datasets:
                axis.plot(dataset.data[TIME_INDEX], dataset.data[VOLTAGE_INDEX], linestyle=dataset.linestyle, label=dataset.description)
            axis.set(xlabel='time (s)', ylabel=r'voltage (V)', title='Evaluations of $v_L(t)$ at $\Delta{}t=%gs$' % delta_t)
        axis.legend(loc='best', fancybox=True, shadow=True)
        axis.grid()
        axis.set_xticks(numpy.arange(dataset.data[TIME_INDEX][0], dataset.data[TIME_INDEX][-1], delta_t/2), minor=True)
        axis.set_xlim(xlims)
        axis.grid(which='minor', alpha=0.2)
        axis.grid(which='major', alpha=0.5)

    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    fig.savefig(name + '.pgf')
    fig.savefig(name + '.png')   

# Flux/current plot generation
def generate_flux_plot():

    # piecewise flux linkage
    currents = [0.0,0.1,1.0]
    flux = [0.0,0.002,0.0024]
    # knee extrapolation
    knee_currents = [0.0, 0.1]
    knee_flux = [1.956e-3, 0.002]

    fig, axis = plt.subplots(1)
    axis.plot(currents, flux, color='blue')
    axis.plot(knee_currents, knee_flux, color='blue', linestyle='dotted', linewidth='0.5')
    axis.axhline(y=0.002, linestyle='dashed', c='green', linewidth='0.5')
    axis.text(0.5, 0.002-0.00008, '0.0020', fontsize=6, va='center', ha='center')
    axis.axhline(y=0.0024, linestyle='dashed', c='green', linewidth='0.5')
    axis.text(0.5, 0.0024-0.00008, '0.0024', fontsize=6, va='center', ha='center')
    axis.axvline(x=0.1, linestyle='dashed', c='green', linewidth='0.5')
    axis.axvline(x=1.0, linestyle='dashed', c='green', linewidth='0.5')
    axis.set(xlabel=r'current $(A)$', ylabel=r'$\lambda{}(i)$ $Wb-t$')
    axis.set_xlim([0.0,1.0])
    axis.set_xticks(list(axis.get_xticks()) + [0.1])
    yticks = list(axis.get_yticks())
    yticks.remove(0.002)
    axis.set_yticks(yticks + [1.956e-3])
    axis.set_ylim([0.0,0.0025])

    fig.set_size_inches(4,2.5)
    fig.tight_layout()
    fig.savefig('flux_plot.pgf')
    fig.savefig('flux_plot.png')

# Main function
def main(args):

    # Read in PSCAD .CSV data
    pscad_dat_0p0001 = Dataset('PSCAD Simulation', linestyle='dotted') # PSCAD data (non-linear inductance, delta_t = 0.1ms, sim_time = 10ms)
    print('*** Opening assignment 8 CSV data file...')
    with open('pscad-dat-0p0001.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                line_count = line_count + 1
                continue
            else:
                pscad_dat_0p0001.data[TIME_INDEX].append(float(row[0]))
                pscad_dat_0p0001.data[VOLTAGE_INDEX].append(float(row[1]))
                pscad_dat_0p0001.data[CURRENT_INDEX].append(float(row[2]))
                line_count = line_count + 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' lines.')
    pscad_dat_0p0003 = Dataset('PSCAD Simulation', linestyle='dotted') # PSCAD data (non-linear inductance, delta_t = 0.1ms, sim_time = 10ms)
    print('*** Opening assignment 8 CSV data file...')
    with open('pscad-dat-0p0001.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                line_count = line_count + 1
                continue
            else:
                pscad_dat_0p0003.data[TIME_INDEX].append(float(row[0]))
                pscad_dat_0p0003.data[VOLTAGE_INDEX].append(float(row[1]))
                pscad_dat_0p0003.data[CURRENT_INDEX].append(float(row[2]))
                line_count = line_count + 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' lines.')

    delta_ts = [0.0001, 0.00015, 0.0002, 0.0003]

    for delta_t in delta_ts:

        print("*** Running simulation ***")
        print("Simulation time = %g" % simulation_time)
        print("Time delta = %g" % delta_t)

        cont_results = Dataset('Continuous Ideal L', linestyle='dashed') # Continuous solution (ideal inductance)
        cont_nl_results = Dataset('Continuous Non-Linear L', linestyle='dashed') # Continuous solution (non-linear inductance)
        trap_results = Dataset('Trapezoidal Ideal L', linestyle='dashdot') # Trapezoidal discretized solution (ideal inductance)
        nonlin_results = Dataset('Trapezoidal Non-Linear L', linestyle='solid') # Trapezoidal discretized solution (non-linear inductance)
        cda_results = Dataset('Trapezoidal+CDA Non-Linear L', linestyle='solid') # Trapezoidal discretized solution + CDA (non-linear inductance)
        
        # PSCAD and continuous only really apply to delta_t == 0.0001
        if delta_t == 0.0001:
            datasets = [pscad_dat_0p0001, cont_results, cont_nl_results, trap_results, nonlin_results, cda_results]
        else:
            datasets = [pscad_dat_0p0003, cont_results, cont_nl_results, trap_results, nonlin_results, cda_results]

        # Set initial values
        i_0 = 0
        v_0 = 10

        # Calculate ideal inductance results using continuous solution
        # Calculate 3000 points for the continuous solution
        for step in range(3000):
            time = 0.0 + (step * simulation_time/3000)
            # calculate continuous solution
            # current
            i_cont = (10/resistance) - ((10/resistance)*math.exp(-(resistance/inductance)*time))
            # inductor voltage
            v_cont = 10 * math.exp((-(resistance)/(inductance))*time)
            # append the result
            cont_results.data[TIME_INDEX].append(0+(step*simulation_time/3000))
            cont_results.data[CURRENT_INDEX].append(i_cont)
            cont_results.data[VOLTAGE_INDEX].append(v_cont)

        # Calculate non-linear inductance results using continuous solution
        # Calculate 3000 points for the continuous solution
        i_cont = 0.0
        for step in range(3000):
            time = 0.0 + (step * simulation_time/3000)
            # calculate continuous solution
            # current
            i_cont = (10/resistance) - ((10/resistance)*math.exp(-(resistance/get_flux(i_cont)[GET_FLUX_INDUCTANCE])*time))
            # inductor voltage
            v_cont = 10 * math.exp((-(resistance)/(get_flux(i_cont)[GET_FLUX_INDUCTANCE]))*time)
            # append the result
            cont_nl_results.data[TIME_INDEX].append(0+(step*simulation_time/3000))
            cont_nl_results.data[CURRENT_INDEX].append(i_cont)
            cont_nl_results.data[VOLTAGE_INDEX].append(v_cont)

        # Calculate ideal inductance results using trapezoidal discretization
        sim_steps = int(simulation_time/delta_t)
        for step in range(sim_steps):
            if step == 0:
                # Use initial conditions
                i_prev = i_0
                v_prev = v_0
            else:
                # Use previous value
                i_prev = i_approx_next
                v_prev = v_approx_next

            # append the result
            trap_results.data[TIME_INDEX].append(0+(step*delta_t))
            trap_results.data[CURRENT_INDEX].append(i_prev)
            trap_results.data[VOLTAGE_INDEX].append(v_prev)

            # perform trapezoidal discretization step
            # current
            i_approx_next = ((i_prev*((2*inductance) - (resistance*delta_t))) + (20 * delta_t))/((resistance*delta_t)+(2*inductance))
            # voltage
            v_approx_next = (((2*inductance)/(delta_t))*(i_approx_next - i_prev)) - v_prev

        # Calculate non-linear inductance results using trapezoidal discretization
        # Nodal analysis
        switch_open_flag = False
        # Start at t = 0
        t = 0
        # Starting conditions
        v1 = 10.0
        v2 = 0.0
        v3 = 10.0
        i30 = 0.0
        r = resistance

        while t < simulation_time:

            # Get values for l, flux, and flux_knee
            l = get_flux(i30)[GET_FLUX_INDUCTANCE]
            flux = get_flux(i30)[GET_FLUX]
            flux_knee = get_flux(i30)[GET_FLUX_KNEE]

            # Check if breaker is opened or closed and set breaker resistance accordingly
            if switch_open_flag == False:
                # Breaker has very small resistance (closed)
                r_sw = 1e-9
            else:
                # Break has very large resistance (open)
                r_sw = 1e9

            # Voltage calculations
            # These voltage calculations are generated in MATLAB
            # by v_next[] = inv(GAA)*hA - inv(GAA)*GAB*vs

            v1_next = 10

            # trapezoidal
            v2_next_trap = (10*(2*l + delta_t*r))/(2*l + delta_t*r + delta_t*r_sw) - (delta_t*r_sw*(v3 + (2*flux)/delta_t - (2*flux_knee)/delta_t))/(2*l + delta_t*r + delta_t*r_sw)
            v3_next_trap = (20*l)/(2*l + delta_t*r + delta_t*r_sw) - (delta_t*(r + r_sw)*(v3 + (2*flux)/delta_t - (2*flux_knee)/delta_t))/(2*l + delta_t*r + delta_t*r_sw)

            # Next currents
            i30_trap_next = (delta_t*v3_next_trap)/(2*l) + (delta_t*(v3 + (2*flux)/delta_t - (2*flux_knee)/delta_t))/(2*l)

            # Append results
            nonlin_results.data[TIME_INDEX].append(t)
            nonlin_results.data[CURRENT_INDEX].append(i30)
            nonlin_results.data[VOLTAGE_INDEX].append(v3)

            # Next iteration
            i30 = i30_trap_next 
            v1 = v1_next
            v2 = v2_next_trap
            v3 = v3_next_trap

            # Go to next time step
            t = t + delta_t

        # Calculate non-linear inductance results using trapezoidal discretization + CDA
        # Nodal analysis
        switch_open_flag = False
        # Start at t = 0
        t = 0
        # Starting conditions
        v1 = 10.0
        v2 = 0.0
        v3 = 10.0
        i30 = 0.0
        r = resistance
        # Get an initial inductance
        l_prev = get_flux(i30)[GET_FLUX_INDUCTANCE]
        # CDA initially inactive
        cda_active = False
        cda_counter = 0
        while t < simulation_time:

            # Get values for l, flux, and flux_knee
            l = get_flux(i30)[GET_FLUX_INDUCTANCE]
            # Check if we need to perform CDA
            if (l != l_prev):
                # We've entered a new region, for the next two half time_deltas
                # Perform backward Euler instead
                cda_active = True
                cda_counter = 2

            flux = get_flux(i30)[GET_FLUX]
            flux_knee = get_flux(i30)[GET_FLUX_KNEE]

            # Check if breaker is opened or closed and set breaker resistance accordingly
            if switch_open_flag == False:
                # Breaker has very small resistance (closed)
                r_sw = 1e-9
            else:
                # Break has very large resistance (open)
                r_sw = 1e9

            # Voltage calculations
            # These voltage calculations are generated in MATLAB
            # by v_next[] = inv(GAA)*hA - inv(GAA)*GAB*vs

            v1_next = 10

            if (cda_active):
                # backward euler
                # voltages
                v2_next = (10*(l + delta_t*r))/(l + delta_t*r + delta_t*r_sw) - (delta_t*r_sw*(flux/delta_t - flux_knee/delta_t))/(l + delta_t*r + delta_t*r_sw)
                v3_next = (10*l)/(l + delta_t*r + delta_t*r_sw) - (delta_t*(r + r_sw)*(flux/delta_t - flux_knee/delta_t))/(l + delta_t*r + delta_t*r_sw)
                # current
                i30_next = (delta_t*(flux/delta_t - flux_knee/delta_t))/l + (delta_t*v3_next)/l
            else:
                # trapezoidal
                # voltages
                v2_next = (10*(2*l + delta_t*r))/(2*l + delta_t*r + delta_t*r_sw) - (delta_t*r_sw*(v3 + (2*flux)/delta_t - (2*flux_knee)/delta_t))/(2*l + delta_t*r + delta_t*r_sw)
                v3_next = (20*l)/(2*l + delta_t*r + delta_t*r_sw) - (delta_t*(r + r_sw)*(v3 + (2*flux)/delta_t - (2*flux_knee)/delta_t))/(2*l + delta_t*r + delta_t*r_sw)
                # current
                i30_next = (delta_t*v3_next)/(2*l) + (delta_t*(v3 + (2*flux)/delta_t - (2*flux_knee)/delta_t))/(2*l)    

            # Append results
            cda_results.data[TIME_INDEX].append(t)
            cda_results.data[CURRENT_INDEX].append(i30)
            cda_results.data[VOLTAGE_INDEX].append(v3)

            # Next iteration
            i30 = i30_next
            v1 = v1_next
            v2 = v2_next
            v3 = v3_next
            l_prev = l

            # Go to next time step
            if (cda_active):
                cda_counter = cda_counter - 1
                # Check if we're done with CDA
                if cda_counter == 0:
                    cda_active = False
                t = t + delta_t/2.0
            else:
                t = t + delta_t

        # Plots for publication
        # Comparisons between all techniques
        generate_flux_plot()
        generate_plots('compare_plots_' + str(delta_t).replace('.','p'), datasets, (0, 0.01), delta_t)
        # Change zoom level depending on time step
        if delta_t == 0.0001:
            x_zoom_min = 0.0002
            x_zoom_max = 0.0005
        elif delta_t == 0.00015:
            x_zoom_min = 0.0002
            x_zoom_max = 0.001
        elif delta_t == 0.0002:
            x_zoom_min = 0.0002
            x_zoom_max = 0.001
        elif delta_t == 0.0003:
            x_zoom_min = 0.0002
            x_zoom_max = 0.002
        else:
            # default zoom
            x_zoom_min = 0.0002
            x_zoom_max = 0.0005
        generate_plots('compare_plots_zoom_' + str(delta_t).replace('.','p'), datasets, (x_zoom_min, x_zoom_max), delta_t)    


if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 8 solution generator')

    args = parser.parse_args()

    main(args)