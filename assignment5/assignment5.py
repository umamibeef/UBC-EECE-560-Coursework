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

# Main function
def main(args):

    # Circuit Constants
    C_zero = 7.5240e-03 * 1e-6 # Farads/km
    G_zero = 2.0000e-08 # Mhos/km

    # List indices
    # CSV data
    FREQ_INDEX = 0
    R_ZERO_INDEX = 1
    L_ZERO_INDEX = 2
    # Calculated impedance
    MAGNITUDE_INDEX = 0
    PHASE_INDEX = 1

    # ZEROES (in Hz)
    zeros = [1.8,107,3.5e3,2.27e5]
    # POLES (in Hz)
    poles = [1.5,90,3e3,2e5]

    # prepopulate data with a list of five empty lists
    data = [[] for i in range(3)]

    # Read in PSCAD .CSV data
    print('*** Opening assignment 4 CSV data file...')
    with open('data_assign04.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                print('Column names are: ' + ', '.join(row))
            else:
                data[FREQ_INDEX].append(float(row[0]))
                data[R_ZERO_INDEX].append(float(row[1])) # Ohms/km
                data[L_ZERO_INDEX].append(float(row[2]) * 1e-3) # Henries/km
            line_count += 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' data points.')

    # Interpolate additional data
    # - We have a total of eight decades, and we want 10 points per decade = 80 datapoints
    # - Create an array of X values to interpolate for, from 1e-1 to 1e-7, base 10 logarithmic increase
    interp_f_data = np.logspace(-1, 7, base=10, num=80)
    # Obtain omega values from frequencies
    interp_w_data = 2*np.pi*interp_f_data
    # Use numpy's interp function to interpolate from existing data
    interp_r_data = np.interp(interp_f_data, data[FREQ_INDEX], data[R_ZERO_INDEX])
    interp_l_data = np.interp(interp_f_data, data[FREQ_INDEX], data[L_ZERO_INDEX])

    num_interp_data_points = len(interp_f_data)

    # Prepare values for interpreted Z(w) magnitude and phase
    interp_impedance_zero = [[],[]]
    for index in range(num_interp_data_points):
        omega = 2*np.pi*interp_f_data[index]
        interp_impedance_zero_val = np.sqrt((interp_r_data[index] + (1j*omega*interp_l_data[index]))/(G_zero + (1j*omega*C_zero)))
        interp_impedance_zero[MAGNITUDE_INDEX].append(np.absolute(interp_impedance_zero_val))
        # print(interp_impedance_zero[MAGNITUDE_INDEX][-1])
        interp_impedance_zero[PHASE_INDEX].append(np.angle(interp_impedance_zero_val))
        # print(interp_impedance_zero[PHASE_INDEX][-1])

    #### Generate Bode plot ####

    # Covert poles and zeros into omega values
    poles_w = [2*np.pi*f for f in poles]
    zeros_w = [2*np.pi*f for f in zeros]

    # Python control library requires the numerator and denominator be expanded polynomials.
    # The following code will prepare the poles and zeros for expansion and extraction of the coeffecients.
    # Create the numerator string
    numerator = ''.join(['(s + %d)'%zero if i is 0 else '*(s + %d)'%zero for i, zero in enumerate(zeros_w)])
    # Create the denominator string
    denominator = ''.join(['(s + %d)'%pole if i is 0 else '*(s + %d)'%pole for i, pole in enumerate(poles_w)])

    # Tell sympy we're using s as a symbol in our equations
    s = sympy.symbols('s')

    # Change numerator and denominator strings into sympy symbolic expressions
    num_poly = sympy.poly(sympy.core.sympify(numerator))
    den_poly = sympy.poly(sympy.core.sympify(denominator))

    # Get polynomial coefficients
    num_coeffs = [float(num) for num in num_poly.coeffs()]
    den_coeffs = [float(den) for den in den_poly.coeffs()]

    # Create our control function from the coefficients
    H = control.tf(num_coeffs, den_coeffs)

    h_mag, h_phase, h_frequencies = control.bode_plot(H, interp_w_data)
    h_mag_scaled = np.array(h_mag) * interp_impedance_zero[MAGNITUDE_INDEX][0]/h_mag[0]
    h_scale_factor = interp_impedance_zero[MAGNITUDE_INDEX][0]/h_mag[0]
    print("Scaled Bode plot by a factor of: " + str(h_scale_factor))

    # Perform partial fraction expansion to obtain our equivalent circuit components
    symbolic_transfer_function = sympy.core.sympify(h_scale_factor)*num_poly/den_poly
    print("Transfer function factor form: \r\n" + str(sympy.latex(sympy.N(symbolic_transfer_function,5))))
    print("Simple fraction expansion of transfer function: \r\n" + str(sympy.latex(sympy.N(sympy.apart(symbolic_transfer_function),5))))

    # Plots for publication
    legend_font_size = 6

    # Plot Z(w) magnitude and phase
    fig, ax = plt.subplots(2)

    ax[0].plot(interp_f_data, interp_impedance_zero[MAGNITUDE_INDEX], label='original')
    ax[0].plot(interp_f_data, h_mag_scaled, color='g', label='rational function approximation')
    ax[0].set(xlabel='Frequency $Hz$', ylabel='Magnitude ($\Omega/km$)', title='$Z_c(\omega{})$ - Magnitude vs. Frequency')
    ax[0].grid(b=True, which='major', color='gray', linestyle='-')
    ax[0].grid(b=True, which='minor', color='gainsboro', linestyle='--')
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    # Plot zero locaitons
    for i in range(len(zeros)):
        ax[0].axvline(x=zeros[i], linestyle='dashed', c='red', linewidth='0.75')
    # Plot pole locations
    for i in range(len(poles)):
        ax[0].axvline(x=poles[i], linestyle='dashed', c='orange', linewidth='0.75')
    handles, labels = ax[0].get_legend_handles_labels()
    handles.append(plt.axvline(x=zeros[0], linestyle='dashed', c='red', linewidth='0.75'))
    labels.append('zero location')
    handles.append(plt.axvline(x=poles[0], linestyle='dashed', c='orange', linewidth='0.75'))
    labels.append('pole location')
    ax[0].legend(handles=handles, labels=labels, loc='upper right', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    
    ax[1].plot(interp_f_data, interp_impedance_zero[PHASE_INDEX], label='original')
    ax[1].plot(interp_f_data, h_phase, color='g', label='rational function approximation')
    ax[1].yaxis.set_major_formatter(tck.FormatStrFormatter('%1.2f $\pi$'))
    ax[1].yaxis.set_major_locator(tck.MultipleLocator(base=1/100))
    ax[1].set(xlabel='Frequency $Hz$', ylabel='Phase ($rad$)', title='$Z_c(\omega{})$ - Phase vs. Frequency')
    ax[1].grid(b=True, which='major', color='gray', linestyle='-')
    ax[1].grid(b=True, which='minor', color='gainsboro', linestyle='--')
    ax[1].set_xscale('log')
    # Plot zero locaitons
    for i in range(len(zeros)):
        ax[1].axvline(x=zeros[i], linestyle='dashed', c='red', linewidth='0.75')
    # Plot pole locations
    for i in range(len(poles)):
        ax[1].axvline(x=poles[i], linestyle='dashed', c='orange', linewidth='0.75')
    handles, labels = ax[1].get_legend_handles_labels()
    handles.append(plt.axvline(x=zeros[0], linestyle='dashed', c='red', linewidth='0.75'))
    labels.append('zero location')
    handles.append(plt.axvline(x=poles[0], linestyle='dashed', c='orange', linewidth='0.75'))
    labels.append('pole location')
    ax[1].legend(handles=handles, labels=labels, loc='lower right', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    fig.savefig('zc_magnitude_and_phase_plot.pgf')
    fig.savefig('zc_magnitude_and_phase_plot.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 5 solution generator')

    args = parser.parse_args()

    main(args)