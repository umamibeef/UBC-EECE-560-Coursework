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

# Main function
def main(args):

    C_zero = 7.5240e-03 * 1e-6 # Farads/km
    C_pos = 1.2027e-02 * 1e-6 # Farads/km
    G_zero = 2.0000e-08 # Mhos/km
    G_pos = 2.0000e-08 # Mhos/km
    length = 300 # km

    FREQ_INDEX = 0
    R_ZERO_INDEX = 1
    L_ZERO_INDEX = 2
    R_POS_INDEX = 3
    L_POS_INDEX = 4

    MAGNITUDE_INDEX = 0
    PHASE_INDEX = 1

    # prepopulate data with a list of five empty lists
    data = [[] for i in range(5)]

    # Read in PSCAD .CSV data
    print('*** Opening assignment 4 CSV data file...')
    with open('data_assign04.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        # Read in row data
        for row in csv_reader:
            if line_count == 0:
                print('Column names are: ' + ', '.join(row))
                line_count += 1
            else:
                data[FREQ_INDEX].append(float(row[0]))
                data[R_ZERO_INDEX].append(float(row[1])) # Ohms/km
                data[L_ZERO_INDEX].append(float(row[2]) * 1e-3) # Henries/km
                data[R_POS_INDEX].append(float(row[3])) # Ohms/km
                data[L_POS_INDEX].append(float(row[4]) * 1e-3) # Henries/km
                line_count += 1
        # Figure out when break switched
        print('Processed ' + str(line_count) + ' lines.')

    num_data_points = len(data[FREQ_INDEX])

    # Prepare values for Z(w) magnitude and phase
    impedance_zero = [[],[]]
    impedance_pos = [[],[]]
    for index in range(num_data_points):
        omega = 2*np.pi*data[FREQ_INDEX][index]
        impedance_zero_val = np.sqrt((data[R_ZERO_INDEX][index] + (1j*omega*data[L_ZERO_INDEX][index]))/(G_zero + (1j*omega*C_zero)))
        impedance_pos_val =  np.sqrt((data[R_POS_INDEX][index] +  (1j*omega*data[L_POS_INDEX][index])) /(G_pos +  (1j*omega*C_pos)))
        # print("F: " + str(data[FREQ_INDEX][index]))
        # print("Omega: " + str(omega))
        # print("R_0: " + str(data[R_ZERO_INDEX][index]))
        # print("L_0: " + str(data[L_ZERO_INDEX][index]))
        # print("C_0: " + str(C_zero))
        # print("G_0: " + str(G_zero))
        # print("R_+: " + str(data[R_POS_INDEX][index]))
        # print("L_+: " + str(data[L_POS_INDEX][index]))
        # print("C_+: " + str(C_pos))
        # print("G_+: " + str(G_pos))
        # print("Zc_0: " + str(impedance_zero_val))
        # print("Zc_0 mag: " + str(np.absolute(impedance_zero_val)))
        # print("Zc_+: " + str(impedance_pos_val))
        # print("Zc_+ mag: " + str(np.absolute(impedance_pos_val)))
        impedance_zero[MAGNITUDE_INDEX].append(np.absolute(impedance_zero_val))
        impedance_zero[PHASE_INDEX].append(np.angle(impedance_zero_val))
        impedance_pos[MAGNITUDE_INDEX].append(np.absolute(impedance_pos_val))
        impedance_pos[PHASE_INDEX].append(np.angle(impedance_pos_val))
        print("\r\n")

    # Prepare values for propagation function magnitude and phase as well
    # Prepare values for attenuation alpha(w) (nepers/km)
    # Prepare values for phase displacement beta(w) (radians/km)
    # Prepare values for propagation speed a(w) (km/s)
    propagation_zero = [[],[]]
    propagation_pos = [[],[]]
    attenuation_zero = []
    attenuation_pos = []
    phase_zero = []
    phase_pos = []
    propagation_speed_zero = []
    propagation_speed_pos = []
    for index in range(num_data_points):
        omega = 2*np.pi*data[FREQ_INDEX][index]
        gamma_zero = np.sqrt((data[R_ZERO_INDEX][index] + 1j*omega*data[L_ZERO_INDEX][index])*(G_zero + 1j*omega*C_zero))
        gamma_pos = np.sqrt((data[R_POS_INDEX][index] + 1j*omega*data[L_POS_INDEX][index])*(G_pos + 1j*omega*C_pos))
        # propagation function magnitude and phase
        propagation_zero[MAGNITUDE_INDEX].append(np.absolute(np.exp(-1*gamma_zero*length)))
        propagation_zero[PHASE_INDEX].append(-np.imag(gamma_zero)*length)
        propagation_pos[MAGNITUDE_INDEX].append(np.absolute(np.exp(-1*gamma_pos*length)))
        propagation_pos[PHASE_INDEX].append(-np.imag(gamma_pos)*length)
        # attenuation (real component of gamma) (nepers/km)
        attenuation_zero.append(np.real(gamma_zero))
        attenuation_pos.append(np.real(gamma_pos))
        # phase displacement (imaginary component of gamma) (radians/km)
        phase_zero.append(np.imag(gamma_zero))
        phase_pos.append(np.imag(gamma_pos))
        # propagation speed (omega/phase_displacement) (km/s)
        propagation_speed_zero.append(omega/phase_zero[-1])
        propagation_speed_pos.append(omega/phase_pos[-1])
        # propagation_speed_zero.append(1/np.sqrt(data[L_ZERO_INDEX][index]*C_zero))
        # propagation_speed_pos.append(1/np.sqrt(data[L_POS_INDEX][index]*C_pos))

    # Plots for publication
    legend_font_size = 6
    # Plot Z(w) magnitude and phase
    fig, ax = plt.subplots(2)
    ax[0].plot(data[FREQ_INDEX], impedance_zero[MAGNITUDE_INDEX], color='b', label='zero sequence')
    ax[0].plot(data[FREQ_INDEX], impedance_pos[MAGNITUDE_INDEX], color='g', label='positive sequence')
    ax[0].set(xlabel='Frequency $Hz$', ylabel='Magnitude ($\Omega/km$)', title='$Z_c$ - Magnitude vs. Frequency')
    ax[0].grid()
    ax[0].set_xscale('log') 
    ax[0].legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax[1].plot(data[FREQ_INDEX], impedance_zero[PHASE_INDEX], color='b', label='zero sequence')
    ax[1].plot(data[FREQ_INDEX], impedance_pos[PHASE_INDEX], color='g', label='positive sequence')
    ax[1].yaxis.set_major_formatter(tck.FormatStrFormatter('%1.1f $\pi$'))
    ax[1].yaxis.set_major_locator(tck.MultipleLocator(base=1/5))
    ax[1].set(xlabel='Frequency $Hz$', ylabel='Phase ($rad$)', title='$Z_c$ - Phase vs. Frequency')
    ax[1].grid()
    ax[1].set_xscale('log')
    ax[1].legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    fig.savefig('zc_magnitude_phase_plots.pgf')
    fig.savefig('zc_magnitude_phase_plots.png')

    # Plot propagation function magnitude and phase
    fig, ax = plt.subplots(2)
    ax[0].plot(data[FREQ_INDEX], propagation_zero[MAGNITUDE_INDEX], color='b', label='zero sequence')
    ax[0].plot(data[FREQ_INDEX], propagation_pos[MAGNITUDE_INDEX], color='g', label='positive sequence')
    ax[0].set(xlabel='Frequency $Hz$', ylabel=r'Magnitude $\left|e^{-\gamma{}l}\right|$', title='$e^{-\gamma{}l}$ - Magnitude vs. Frequency')
    ax[0].grid()
    ax[0].set_xscale('log')
    ax[0].legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    ax[1].plot(data[FREQ_INDEX], propagation_zero[PHASE_INDEX], color='b', label='zero sequence')
    ax[1].plot(data[FREQ_INDEX], propagation_pos[PHASE_INDEX], color='g', label='positive sequence')
    ax[1].set(xlabel='Frequency $Hz$', ylabel=r'Phase $\phi{}=\beta{}l=\omega{}\tau$ ($rad$)', title='$e^{-\gamma{}l}$ - Phase vs. Frequency')
    ax[1].grid()
    ax[1].set_xscale('log')
    ax[1].legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,8)
    fig.tight_layout()
    fig.savefig('prop_magnitude_phase_plots.pgf')
    fig.savefig('prop_magnitude_phase_plots.png')

    # Plot propagation function magnitude and phase (no long for frequency)
    fig, ax = plt.subplots(1)
    ax.plot(data[FREQ_INDEX], propagation_zero[PHASE_INDEX], color='b', label='zero sequence')
    ax.plot(data[FREQ_INDEX], propagation_pos[PHASE_INDEX], color='g', label='positive sequence')
    ax.set(xlabel='Frequency $Hz$', ylabel=r'Phase $\phi{}=\beta{}l=\omega{}\tau$ ($rad$)', title='$e^{-\gamma{}l}$ - Phase vs. Frequency')
    ax.grid()
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('prop_phase_plot_nolog.pgf')
    fig.savefig('prop_phase_plot_nolog.png')

    # Plot attenuation (real component of gamma) (nepers/km)
    fig, ax = plt.subplots()
    ax.plot(data[FREQ_INDEX], attenuation_zero, color='b', label='zero sequence')
    ax.plot(data[FREQ_INDEX], attenuation_pos, color='g', label='positive sequence')
    ax.set(xlabel='Frequency $Hz$', ylabel=r'Attenuation $\alpha{}(\omega)$ $(nepers/km)$', title=r'Attenuation $\alpha{}(\omega)$ vs. Frequency')
    ax.grid()
    ax.set_xscale('log')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('attenuation_plots.pgf')
    fig.savefig('attenuation_plots.png')

    # Plot phase displacement beta(w) (radians/km)
    fig, ax = plt.subplots()
    ax.plot(data[FREQ_INDEX], phase_zero, color='b', label='zero sequence')
    ax.plot(data[FREQ_INDEX], phase_pos, color='g', label='positive sequence')
    ax.set(xlabel='Frequency $Hz$', ylabel=r'Phase Displacement $\beta{}(\omega)$ $(rad/km)$', title=r'Phase Displacement $\beta{}(\omega)$ vs. Frequency')
    ax.grid()
    ax.set_xscale('log')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('phase_displacement_plots.pgf')
    fig.savefig('phase_displacement_plots.png')

    # Plot propagation speed a(w) (km/s)
    fig, ax = plt.subplots()
    ax.plot(data[FREQ_INDEX], propagation_speed_zero, color='b', label='zero sequence')
    ax.plot(data[FREQ_INDEX], propagation_speed_pos, color='g', label='positive sequence')
    ax.set(xlabel='Frequency $Hz$', ylabel=r'Propagation Speed $a(\omega)$ $(km/s)$', title=r'Propagation Speed $a(\omega)$ vs. Frequency')
    ax.grid()
    ax.set_xscale('log')
    ax.legend(loc='best', prop={'size':legend_font_size}, fancybox=True, shadow=True)
    fig.set_size_inches(6.5,3.5)
    fig.tight_layout()
    fig.savefig('propagation_speed_plots.pgf')
    fig.savefig('propagation_speed_plots.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 4 solution generator')

    args = parser.parse_args()

    main(args)