import argparse
import math
import matplotlib
import matplotlib.pyplot as plt

# Matplotlib export settings
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "font.size": 10 ,
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
    "pgf.rcfonts": False     # don't setup fonts from rc parameters
})

# Main function
def main(args):
    # Circuit parameters
    inductance_L = 0.02 # 20 mH
    resistance_R = 10 # 10 Ohms
    voltage_E = 10 # 10 Volts

    print("Simulation time = %g" % args.simulation_time)
    print("Will run for the following time deltas:")
    for delta in args.time_delta:
        print("Time delta = %g" % delta)

    cont_results = [[],[],[]]
    trap_results = [[[],[],[]]]
    back_results = [[[],[],[]]]
    # indices for the results tuple
    TIME_INDEX = 0
    CURRENT_INDEX = 1
    VOLTAGE_INDEX = 2

    # Calculate results using continuous solution
    print("Simulating using continuous solution")
    # Calculate 3000 points for the continuous solution
    for step in range(3000):
        time = 0.0 + (step * args.simulation_time/3000)
        # calculate continuous solution
        # current
        i_continuous_next = (10/resistance_R) - ((10/resistance_R)*math.exp(-(resistance_R/inductance_L)*time))
        # inductor voltage
        v_continuous_next = 10 * math.exp((-(resistance_R)/(inductance_L))*time)
        # append the result
        cont_results[TIME_INDEX].append(0+(step*args.simulation_time/3000))
        cont_results[CURRENT_INDEX].append(i_continuous_next)
        cont_results[VOLTAGE_INDEX].append(v_continuous_next)

    # Calculate results using trapezoidal discretization
    print("Simulating using trapezoidal discretization")
    for delta_index, delta in enumerate(args.time_delta):
        sim_steps = int(args.simulation_time/delta)
        for step in range(sim_steps):
            if step == 0:
                # Use initial conditions
                i_prev = args.i_0
                v_prev = args.v_0
            else:
                # Use previous value
                i_prev = trap_results[delta_index][CURRENT_INDEX][step - 1]
                v_prev = trap_results[delta_index][VOLTAGE_INDEX][step - 1]
            # perform trapezoidal discretization step
            # current
            i_approx_next = ((i_prev*((2*inductance_L) - (resistance_R*delta))) + (20 * delta))/((resistance_R*delta)+(2*inductance_L))
            # voltage
            v_approx_next = (((2*inductance_L)/(delta))*(i_approx_next - i_prev)) - v_prev
            # append the result
            trap_results[delta_index][TIME_INDEX].append(0+(step*delta))
            trap_results[delta_index][CURRENT_INDEX].append(i_approx_next)
            trap_results[delta_index][VOLTAGE_INDEX].append(v_approx_next)
        # append a new set of results for the next delta
        trap_results.append([[],[],[]])

    # Calculate results using backwards euler discretization
    print("Simulating using backward euler discretization")
    for delta_index, delta in enumerate(args.time_delta):
        sim_steps = int(args.simulation_time/delta)
        for step in range(sim_steps):
            if step == 0:
                # Use initial conditions
                i_prev = args.i_0
            else:
                # Use previous value
                i_prev = back_results[delta_index][CURRENT_INDEX][step - 1]
            # perform backward eurler discretization step
            # current
            i_approx_next = ((inductance_L*i_prev)+(10*delta))/((resistance_R*delta) + inductance_L)
            # voltage
            v_approx_next = ((inductance_L)/(delta))*(i_approx_next - i_prev)
            # append the result
            back_results[delta_index][TIME_INDEX].append(0+(step*delta))
            back_results[delta_index][CURRENT_INDEX].append(i_approx_next)
            back_results[delta_index][VOLTAGE_INDEX].append(v_approx_next)
        # append a new set of results for the next delta
        back_results.append([[],[],[]])

    # Plots for publication
    # Comparisons between all techniques
    for delta_index, delta in enumerate(args.time_delta):
        # Current plot
        fig, ax = plt.subplots(2)
        ax[0].plot(cont_results[TIME_INDEX], cont_results[CURRENT_INDEX], linestyle='dotted', label='Continuous')
        ax[0].plot(trap_results[delta_index][TIME_INDEX], trap_results[delta_index][CURRENT_INDEX], label='Trapezoidal')
        ax[0].plot(back_results[delta_index][TIME_INDEX], back_results[delta_index][CURRENT_INDEX], label='Backward Euler    ')
        ax[0].set(xlabel='time (s)', ylabel='current (A)', title='Step-By-Step Approximations of $i(t)$')
        ax[0].legend(loc='best', fancybox=True, shadow=True)
        ax[0].grid()
        # Voltage plot
        ax[1].plot(cont_results[TIME_INDEX], cont_results[VOLTAGE_INDEX], linestyle='dotted', label='Continuous')
        ax[1].plot(trap_results[delta_index][TIME_INDEX], trap_results[delta_index][VOLTAGE_INDEX], label='Trapezoidal')
        ax[1].plot(back_results[delta_index][TIME_INDEX], back_results[delta_index][VOLTAGE_INDEX], label='Backward Euler    ')
        ax[1].set(xlabel='time (s)', ylabel='voltage (V)', title='Step-By-Step Approximations of $v_L(t)$')
        ax[1].legend(loc='best', fancybox=True, shadow=True)
        ax[1].grid()
        fig.set_size_inches(6.5,8)
        fig.tight_layout()
        # plt.show()
        fig.savefig('compare_plot_' + str(delta).replace('.','p') + '.pgf')
        fig.savefig('compare_plot_' + str(delta).replace('.','p') + '.png')
    # Comparisons between individual techniques
    # Trapezoidal
    # Current plot
    fig, ax = plt.subplots(2)
    ax[0].plot(cont_results[TIME_INDEX], cont_results[CURRENT_INDEX], linestyle='dotted', label='Continuous')
    for delta_index, delta in enumerate(args.time_delta):
        ax[0].plot(trap_results[delta_index][TIME_INDEX], trap_results[delta_index][CURRENT_INDEX], label='Trapezoidal ($\Delta{}t = %g s$)'%delta)
        ax[0].set(xlabel='time (s)', ylabel='current (A)', title='Step-By-Step Approximations of $i(t)$')
        ax[0].legend(loc='best', fancybox=True, shadow=True)
    ax[0].grid()
    # Voltage plot
    ax[1].plot(cont_results[TIME_INDEX], cont_results[VOLTAGE_INDEX], linestyle='dotted', label='Continuous')
    for delta_index, delta in enumerate(args.time_delta):
        ax[1].plot(trap_results[delta_index][TIME_INDEX], trap_results[delta_index][VOLTAGE_INDEX], label='Trapezoidal ($\Delta{}t = %g s$)'%delta)
        ax[1].set(xlabel='time (s)', ylabel='voltage (V)', title='Step-By-Step Approximations of $v_L(t)$')
        ax[1].legend(loc='best', fancybox=True, shadow=True)
    ax[1].grid()
    # Backward Euler
    # Current plot
    fig, ax = plt.subplots(2)
    ax[0].plot(cont_results[TIME_INDEX], cont_results[CURRENT_INDEX], linestyle='dotted', label='Continuous')
    for delta_index, delta in enumerate(args.time_delta):
        ax[0].plot(back_results[delta_index][TIME_INDEX], back_results[delta_index][CURRENT_INDEX], label='Backward Euler ($\Delta{}t = %g s$)'%delta)
        ax[0].set(xlabel='time (s)', ylabel='current (A)', title='Step-By-Step Approximations of $i(t)$')
        ax[0].legend(loc='best', fancybox=True, shadow=True)
    ax[0].grid()
    # Voltage plot
    ax[1].plot(cont_results[TIME_INDEX], cont_results[VOLTAGE_INDEX], linestyle='dotted', label='Continuous')
    for delta_index, delta in enumerate(args.time_delta):
        ax[1].plot(back_results[delta_index][TIME_INDEX], back_results[delta_index][VOLTAGE_INDEX], label='Backward Euler ($\Delta{}t = %g s$)'%delta)
        ax[1].set(xlabel='time (s)', ylabel='voltage (V)', title='Step-By-Step Approximations of $v_L(t)$')
        ax[1].legend(loc='best', fancybox=True, shadow=True)
    ax[1].grid()
    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    # plt.show()
    fig.savefig('backeuler_plots.pgf')
    fig.savefig('backeuler_plots.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 1a solution generator')
    # allow the user to input multiple time steps if they desire (e.g. 0.1 ms and 0.8 ms simultaneously)
    parser.add_argument('-d', type=float, dest='time_delta', nargs='+', action='store',
        help='time delta to use in discretization')
    parser.add_argument('-l', type=float, dest='simulation_time', action='store',
        help='time to run the simulation for')
    parser.add_argument('-i', type=float, dest='i_0', action='store',
        help='i at i(0)')
    parser.add_argument('-v', type=float, dest='v_0', action='store',
        help='v_L at v_L(0)')
    parser.add_argument('-o', type=str, default='output.csv', dest='output_file', action='store',
        help='path and name of the file containing the results')

    args = parser.parse_args()

    main(args)