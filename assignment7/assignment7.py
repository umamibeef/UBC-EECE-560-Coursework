import argparse
import csv
import matplotlib
import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import numpy as np
 
 # Matplotlib export settings
matplotlib.use('pgf')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({
    'pgf.texsystem': 'pdflatex',
    'font.size': 10,
    'font.family': 'serif',  # use serif/main font for text elements
    'text.usetex': True,     # use inline math for ticks
    'pgf.rcfonts': False     # don't setup fonts from rc parameters
})

# Constants
C = 1
delta_t = 1
step_v = [0,0,0,1,1,1,1,1,1,1,1,1,1,1] # step voltages
ramp_v = [0,0,0,delta_t,2*delta_t,3*delta_t,4*delta_t,5*delta_t,6*delta_t,7*delta_t,8*delta_t,9*delta_t,10*delta_t,11*delta_t]
start_index = 2 # we start at t = 0

# Function that calculates the trapezoidal approximation of a capacitor
def trapezoidal(v_t, v_t_min_delta_t, i_t_min_delta_t):

    i_t = (2*C/delta_t)*v_t - (2*C/delta_t)*v_t_min_delta_t - i_t_min_delta_t

    return i_t

# Function that calculates the backward euler approximation of a capacitor
def backward_euler(v_t, v_t_min_delta_t):

    i_t = (C/delta_t)*v_t - (C/delta_t)*v_t_min_delta_t

    return i_t

# Function that calculates the forward euler discretized, frequency domain admittance of an inductor
def forward_euler(v_t, v_t_plus_delta_t):

    i_t = -(C/delta_t)*v_t + (C/delta_t)*v_t_plus_delta_t

    return i_t

# Function that calculates the gear's second order discretized, frequency domain impedance of an inductor
def gears_second_order(v_t, v_t_min_delta_t, v_t_min_2_delta_t):
    
    i_t = (3*C/(2*delta_t))*(v_t - (4/3)*v_t_min_delta_t + (1/3)*v_t_min_2_delta_t)

    return i_t

# Plot generator
def create_plot(name, y, modifier=None):
    print('Plotting ' + name)

    # x-index
    index_t = [-2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9] # time indices
    labels_t = [r'$%d\Delta{}t$'%index for index in index_t]

    current_ticks = []
    current_minor_ticks = []
    current_labels = []

    # plot specific trims
    if name == 'trapezoidal_step':
        current_ticks = [-2.0, 0, 2.0]
        current_labels = [r'$-\frac{2C}{\Delta{}t}$', r'$0$', r'$\frac{2C}{\Delta{}t}$']
    elif name == 'trapezoidal_ramp':
        current_ticks = [0, 2.0]
        current_labels = [r'$0$', r'$2C$']
        pass
    elif name == 'backward_euler_step':
        current_ticks = [0, 1.0]
        current_labels = [r'$0$', r'$\frac{C}{\Delta{}t}$']
        pass
    elif name == 'backward_euler_ramp':
        current_ticks = [0, 1.0]
        current_labels = [r'$0$', r'$C$']
        pass
    elif name == 'forward_euler_step':
        current_ticks = [0, 1.0]
        current_labels = [r'$0$', r'$\frac{C}{\Delta{}t}$']
        pass
    elif name == 'forward_euler_ramp':
        current_ticks = [0, 1.0]
        current_labels = [r'$0$', r'$C$']
        pass
    elif name == 'gears_second_order_step':
        current_ticks = [-0.5, 0, 1.5]
        current_labels = [r'$-\frac{C}{2\Delta{}t}$', r'$0$', r'$\frac{3C}{2\Delta{}t}$']
        pass
    elif name == 'gears_second_order_ramp':
        current_ticks = [0, 1, 1.5]
        current_labels = [r'$0$', r'$C$', r'$\frac{3C}{2}$']
        pass

    fig, ax = plt.subplots(2)

    if name[-4:] == 'step':
        ax[0].plot(index_t, step_v[:-2], marker='.')
        ax[0].set_yticks([0, 1])
        ax[0].set_yticks(np.arange(0, 1, 0.25), minor=True)
    else:
        ax[0].plot(index_t, ramp_v[:-2], marker='.')
        # slope is 1 so same as x, except for below zero
        ax[0].set_yticks(index_t[2:])
        ax[0].set_yticklabels(labels_t[2:])
    ax[0].set_xticks(index_t)
    ax[0].set_xticklabels(labels_t)
    ax[0].set_title('Voltage Input Stimulus')
    ax[0].set_xlabel('$time$')
    ax[0].set_ylabel('$v(t)$')
    ax[0].grid(which='minor', alpha=0.2)
    ax[0].grid(which='major', alpha=0.2)

    ax[1].plot(index_t, y, marker='.')
    ax[1].set_yticks(current_ticks)
    ax[1].set_yticks(np.arange(current_ticks[0], current_ticks[-1], (current_ticks[-1]-current_ticks[0])/4), minor=True)
    ax[1].set_yticklabels(current_labels)
    ax[1].set_xticks(index_t)
    ax[1].set_xticklabels(labels_t)
    ax[1].set_title('Current Output Response')
    ax[1].set_xlabel('$time$')
    ax[1].set_ylabel('$i(t)$')
    ax[1].grid(which='minor', alpha=0.2)
    ax[1].grid(which='major', alpha=0.2)

    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    fig.savefig(name+'.pgf')
    fig.savefig(name+'.png')

# Main function
def main(args):

    trapezoidal_step_current_vals = [0,0]
    trapezoidal_ramp_current_vals = [0,0]
    backward_euler_step_current_vals = [0,0]
    backward_euler_ramp_current_vals = [0,0]
    forward_euler_step_current_vals = [0,0]
    forward_euler_ramp_current_vals = [0,0]
    gears_second_order_step_current_vals = [0,0]
    gears_second_order_ramp_current_vals = [0,0]

    # plot 10 points
    for index in range(10):
        trapezoidal_step_current_vals.append(trapezoidal(step_v[start_index+index], step_v[start_index+index-1], trapezoidal_step_current_vals[start_index+index-1]))
        trapezoidal_ramp_current_vals.append(trapezoidal(ramp_v[start_index+index], ramp_v[start_index+index-1], trapezoidal_ramp_current_vals[start_index+index-1]))
        backward_euler_step_current_vals.append(backward_euler(step_v[start_index+index], step_v[start_index+index-1]))
        backward_euler_ramp_current_vals.append(backward_euler(ramp_v[start_index+index], ramp_v[start_index+index-1]))
        forward_euler_step_current_vals.append(forward_euler(step_v[start_index+index], step_v[start_index+index+1]))
        forward_euler_ramp_current_vals.append(forward_euler(ramp_v[start_index+index], ramp_v[start_index+index+1]))
        gears_second_order_step_current_vals.append(gears_second_order(step_v[start_index+index], step_v[start_index+index-1], step_v[start_index+index-2]))
        gears_second_order_ramp_current_vals.append(gears_second_order(ramp_v[start_index+index], ramp_v[start_index+index-1], ramp_v[start_index+index-2]))

    # generate plots
    create_plot('trapezoidal_step', trapezoidal_step_current_vals)
    create_plot('trapezoidal_ramp', trapezoidal_ramp_current_vals)
    create_plot('backward_euler_step', backward_euler_step_current_vals)
    create_plot('backward_euler_ramp', backward_euler_ramp_current_vals)
    create_plot('forward_euler_step', forward_euler_step_current_vals)
    create_plot('forward_euler_ramp', forward_euler_ramp_current_vals)
    create_plot('gears_second_order_step', gears_second_order_step_current_vals)
    create_plot('gears_second_order_ramp', gears_second_order_ramp_current_vals)

    print('Done.')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 7 solution generator')

    args = parser.parse_args()

    main(args)