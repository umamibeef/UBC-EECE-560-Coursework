import argparse
import csv
import matplotlib
import matplotlib.ticker as tck
import matplotlib.pyplot as plt
import numpy as np
import control
import sympy
 
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

# inductance is cancelled out algebraically
# delta can be anything, as long as the nyquist frequency chosen is based on it
delta_t = 1e-6
nyquist_frequency = 1.0/(2.0*delta_t)

# Function that calculates the continuous, frequency domain impedance of an inductor
def z_real(omega):

    return 1j*omega

# Function that calculates the continuous, frequency domain admittance of an inductor
def y_real(omega):

    return 1.0/z_real(omega)

# Function that calculates the trapezoidal discretized, frequency domain admittance of an inductor
def y_trapezoidal(omega):

    return (delta_t/2.0)*(np.exp(1j*omega*delta_t) + 1)/(np.exp(1j*omega*delta_t) - 1)

# Function that calculates the backward euler discretized, frequency domain admittance of an inductor
def y_backward_euler(omega):

    return (delta_t)*(np.exp(1j*omega*delta_t))/(np.exp(1j*omega*delta_t) - 1)

# Function that calculates the forward euler discretized, frequency domain admittance of an inductor
def y_forward_euler(omega):

    return (delta_t)*(1/(np.exp(1j*omega*delta_t) - 1))

# Function that calculates the gear's second order discretized, frequency domain impedance of an inductor
def z_gears_second_order(omega):
    
    return (1/(2*delta_t))*(np.exp(-1j*omega*delta_t) - 1)*(np.exp(-1j*omega*delta_t) - 3)

# Main function
def main(args):

    # Create a list of frequencies to plot against
    frequencies = np.linspace(0.1, nyquist_frequency, num = 1000)

    # Get magnitudes
    # Create lists
    y_real_vals = []
    z_real_vals = []
    y_trapezoidal_vals = []
    y_backward_euler_vals = []
    y_forward_euler_vals = []
    z_gears_second_order_vals = []
    # Fill lists
    for frequency in frequencies:
        omega = 2*np.pi*frequency
        y_real_vals.append(y_real(omega))
        z_real_vals.append(z_real(omega))
        y_trapezoidal_vals.append(y_trapezoidal(omega))
        y_backward_euler_vals.append(y_backward_euler(omega))
        y_forward_euler_vals.append(y_forward_euler(omega))
        z_gears_second_order_vals.append(z_gears_second_order(omega))
    # Convert lists to numpy arrays
    y_real_vals = np.array(y_real_vals)
    z_real_vals = np.array(z_real_vals)
    y_trapezoidal_vals = np.array(y_trapezoidal_vals)
    y_backward_euler_vals = np.array(y_backward_euler_vals)
    y_forward_euler_vals = np.array(y_forward_euler_vals)
    z_gears_second_order_vals = np.array(z_gears_second_order_vals)

    # Plots for publication
    legend_font_size = 6

    # Plot Z(w) magnitude and phase
    fig, ax = plt.subplots(2)

    ax[0].plot(frequencies*delta_t, np.abs(y_trapezoidal_vals/y_real_vals), label='Trapezoidal')
    ax[0].plot(frequencies*delta_t, np.abs(y_backward_euler_vals/y_real_vals), label='Backward Euler')
    ax[0].plot(frequencies*delta_t, np.abs(y_forward_euler_vals/y_real_vals), linestyle='--', label='Forward Euler')
    ax[0].plot(frequencies*delta_t, np.abs((1.0/z_gears_second_order_vals)/(1.0/z_real_vals)), label='Gear\'s Second Order')
    ax[0].xaxis.set_minor_locator(tck.MultipleLocator(0.01))
    ax[0].yaxis.set_major_locator(tck.MultipleLocator(0.5))
    ax[0].yaxis.set_minor_locator(tck.MultipleLocator(0.1))
    ax[0].set(xlabel='$\\frac{f}{f_{Nyquist}}$ Frequency (p.u.)', ylabel='$\\frac{H_e(\omega{})}{H(\omega{}}$ Magnitude Distortion (p.u.)', title='Magnitude and Phase Distortion vs. Frequency')
    ax[0].grid(b=True, which='major', color='gray', linestyle='-')
    ax[0].grid(b=True, which='minor', color='gainsboro', linestyle='dotted')
    ax[0].legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)

    ax[1].plot(frequencies*delta_t, np.imag(y_trapezoidal_vals/y_real_vals), label='Trapezoidal')
    ax[1].plot(frequencies*delta_t, np.imag(y_backward_euler_vals/y_real_vals), label='Backward Euler')
    ax[1].plot(frequencies*delta_t, np.imag(y_forward_euler_vals/y_real_vals), label='Forward Euler')
    ax[1].plot(frequencies*delta_t, np.imag((1.0/z_gears_second_order_vals)/(1.0/z_real_vals)), label='Gear\'s Second Order')
    ax[1].xaxis.set_minor_locator(tck.MultipleLocator(0.01))
    ax[1].yaxis.set_major_formatter(tck.FormatStrFormatter('%1.1f $\pi$'))
    ax[1].yaxis.set_major_locator(tck.MultipleLocator(base=1/4))
    ax[1].set(xlabel='$\\frac{f}{f_{Nyquist}}$ Frequency (p.u.)', ylabel='$\\frac{H_e(\omega{})}{H(\omega{}}$ Phase Distortion (radians)')
    ax[1].grid(b=True, which='major', color='gray', linestyle='-')
    ax[1].grid(b=True, which='minor', color='gainsboro', linestyle='dotted')
    ax[1].legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)

    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    fig.savefig('transfer_function_magnitude_and_phase_plot.pgf')
    fig.savefig('transfer_function_magnitude_and_phase_plot.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 6 solution generator')

    args = parser.parse_args()

    main(args)