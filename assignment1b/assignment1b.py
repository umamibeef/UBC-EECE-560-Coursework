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

# 0.8 ms Current
# ----------------------------------------
# Domain, Current
pscad_dat_current_0p8 = [
    [0.0,0.236111111111],
    [0.0008,0.490740740741],
    [0.0016,0.66049382716],
    [0.0024,0.77366255144],
    [0.0032,0.849108367627],
    [0.004,0.899405578418],
    [0.0048,0.932937052279],
    [0.0056,0.955291368186],
    [0.0064,0.970194245457],
    [0.0072,0.980129496971],
    [0.008,0.986752997981],
    [0.0088,0.991168665321],
    [0.0096,0.994112443547],
    [0.0104,0.996074962365],
    [0.0112,0.996074962365],
]
# transpose
pscad_dat_current_0p8 = map(list, zip(*pscad_dat_current_0p8))

# 0.8 ms Voltage
# ----------------------------------------
# Domain, Voltage
pscad_dat_voltage_0p8 = [
    [0.0,7.63888888889],
    [0.0008,5.09259259259],
    [0.0016,3.3950617284],
    [0.0024,2.2633744856],
    [0.0032,1.50891632373],
    [0.004,1.00594421582],
    [0.0048,0.670629477214],
    [0.0056,0.447086318143],
    [0.0064,0.298057545428],
    [0.0072,0.198705030286],
    [0.008,0.13247002019],
    [0.0088,0.0883133467936],
    [0.0096,0.0588755645291],
    [0.0104,0.0392503763527],
    [0.0112,0.0392503763527],
]
# transpose
pscad_dat_voltage_0p8 = map(list, zip(*pscad_dat_voltage_0p8))

# 0.1 ms Current
# ----------------------------------------
# Domain, Current
pscad_dat_current_0p1 = [
    [0.0,0.0362879238548],
    [0.0001,0.0832982690327],
    [0.0002,0.128015426641],
    [0.0003,0.170551259488],
    [0.0004,0.211012173659],
    [0.0005,0.2494993847],
    [0.0006,0.286109170812],
    [0.0007,0.320933113699],
    [0.0008,0.354058327665],
    [0.0009,0.385567677535],
    [0.001,0.415539985948],
    [0.0011,0.444050230536],
    [0.0012,0.471169731486],
    [0.0013,0.49696632995],
    [0.0014,0.521504557757],
    [0.0015,0.544845798842],
    [0.0016,0.567048442801],
    [0.0017,0.588168030957],
    [0.0018,0.608257395301],
    [0.0019,0.627366790652],
    [0.002,0.645544020376],
    [0.0021,0.662834555967],
    [0.0022,0.679281650798],
    [0.0023,0.69492644832],
    [0.0024,0.709808084988],
    [0.0025,0.723963788159],
    [0.0026,0.737428969224],
    [0.0027,0.750237312189],
    [0.0028,0.762420857936],
    [0.0029,0.774010084378],
    [0.003,0.785033982701],
    [0.0031,0.795520129886],
    [0.0032,0.805494757697],
    [0.0033,0.814982818297],
    [0.0034,0.824008046673],
    [0.0035,0.832593020006],
    [0.0036,0.840759214152],
    [0.0037,0.848527057364],
    [0.0038,0.855915981395],
    [0.0039,0.862944470107],
    [0.004,0.869630105712],
    [0.0041,0.87598961275],
    [0.0042,0.882038899933],
    [0.0043,0.887793099937],
    [0.0044,0.893266607257],
    [0.0045,0.89847311422],
    [0.0046,0.903425645233],
    [0.0047,0.908136589368],
    [0.0048,0.91261773135],
    [0.0049,0.916880281041],
    [0.005,0.920934901478],
    [0.0051,0.924791735552],
    [0.0052,0.928460431379],
    [0.0053,0.931950166433],
    [0.0054,0.93526967051],
    [0.0055,0.938427247558],
    [0.0056,0.941430796458],
    [0.0057,0.944287830777],
    [0.0058,0.947005497568],
    [0.0059,0.949590595248],
    [0.006,0.952049590602],
    [0.0061,0.954388634962],
    [0.0062,0.956613579598],
    [0.0063,0.95872999035],
    [0.0064,0.960743161552],
    [0.0065,0.962658129281],
    [0.0066,0.964479683951],
    [0.0067,0.966212382294],
    [0.0068,0.967860558768],
    [0.0069,0.969428336389],
    [0.007,0.970919637053],
    [0.0071,0.972338191343],
    [0.0072,0.973687547863],
    [0.0073,0.974971082113],
    [0.0074,0.976192004937],
    [0.0075,0.97735337055],
    [0.0076,0.978458084182],
    [0.0077,0.979508909344],
    [0.0078,0.980508474741],
    [0.0079,0.981459280852],
    [0.008,0.982363706176],
    [0.0081,0.983224013192],
    [0.0082,0.984042354012],
    [0.0083,0.984820775767],
    [0.0084,0.98556122573],
    [0.0085,0.986265556182],
    [0.0086,0.986935529051],
    [0.0087,0.987572820317],
    [0.0088,0.988179024204],
    [0.0089,0.98875565717],
    [0.009,0.989304161698],
    [0.0091,0.989825909908],
    [0.0092,0.990322206985],
    [0.0093,0.99079429445],
    [0.0094,0.991243353257],
    [0.0095,0.991670506757],
    [0.0096,0.9920768235],
    [0.0097,0.992463319915],
    [0.0098,0.992830962846],
    [0.0099,0.993180671975],
    [0.01,0.993513322123],
    [0.0101,0.993513322123],
]
# transpose
pscad_dat_current_0p1 = map(list, zip(*pscad_dat_current_0p1))

# 0.1 ms Voltage
# ----------------------------------------
# Domain, Voltage
pscad_dat_voltage_0p1 = [
    [0.0,9.63712076145],
    [0.0001,9.16701730967],
    [0.0002,8.71984573359],
    [0.0003,8.29448740512],
    [0.0004,7.88987826341],
    [0.0005,7.505006153],
    [0.0006,7.13890829188],
    [0.0007,6.79066886301],
    [0.0008,6.45941672335],
    [0.0009,6.14432322465],
    [0.001,5.84460014052],
    [0.0011,5.55949769464],
    [0.0012,5.28830268514],
    [0.0013,5.0303367005],
    [0.0014,4.78495442243],
    [0.0015,4.55154201158],
    [0.0016,4.32951557199],
    [0.0017,4.11831969043],
    [0.0018,3.91742604699],
    [0.0019,3.72633209348],
    [0.002,3.54455979624],
    [0.0021,3.37165444033],
    [0.0022,3.20718349202],
    [0.0023,3.0507355168],
    [0.0024,2.90191915012],
    [0.0025,2.76036211841],
    [0.0026,2.62571030776],
    [0.0027,2.49762687811],
    [0.0028,2.37579142064],
    [0.0029,2.25989915622],
    [0.003,2.14966017299],
    [0.0031,2.04479870114],
    [0.0032,1.94505242303],
    [0.0033,1.85017181703],
    [0.0034,1.75991953327],
    [0.0035,1.67406979994],
    [0.0036,1.59240785848],
    [0.0037,1.51472942636],
    [0.0038,1.44084018605],
    [0.0039,1.37055529893],
    [0.004,1.30369894288],
    [0.0041,1.2401038725],
    [0.0042,1.17961100067],
    [0.0043,1.12206900063],
    [0.0044,1.06733392743],
    [0.0045,1.0152688578],
    [0.0046,0.965743547665],
    [0.0047,0.918634106316],
    [0.0048,0.873822686496],
    [0.0049,0.831197189593],
    [0.005,0.790650985223],
    [0.0051,0.75208264448],
    [0.0052,0.715395686213],
    [0.0053,0.680498335666],
    [0.0054,0.647303294902],
    [0.0055,0.615727524419],
    [0.0056,0.585692035423],
    [0.0057,0.557121692231],
    [0.0058,0.529945024318],
    [0.0059,0.504094047522],
    [0.006,0.479504093984],
    [0.0061,0.456113650375],
    [0.0062,0.433864204015],
    [0.0063,0.412700096502],
    [0.0064,0.392568384478],
    [0.0065,0.373418707186],
    [0.0066,0.355203160494],
    [0.0067,0.337876177056],
    [0.0068,0.321394412321],
    [0.0069,0.30571663611],
    [0.007,0.290803629471],
    [0.0071,0.27661808657],
    [0.0072,0.263124521371],
    [0.0073,0.250289178865],
    [0.0074,0.238079950628],
    [0.0075,0.2264662945],
    [0.0076,0.215419158183],
    [0.0077,0.204910906564],
    [0.0078,0.194915252585],
    [0.0079,0.185407191484],
    [0.008,0.176362938241],
    [0.0081,0.167759868082],
    [0.0082,0.159576459883],
    [0.0083,0.151792242328],
    [0.0084,0.144387742702],
    [0.0085,0.13734443818],
    [0.0086,0.130644709488],
    [0.0087,0.124271796831],
    [0.0088,0.118209757961],
    [0.0089,0.112443428304],
    [0.009,0.106958383021],
    [0.0091,0.101740900922],
    [0.0092,0.0967779301457],
    [0.0093,0.0920570555044],
    [0.0094,0.087566467431],
    [0.0095,0.0832949324344],
    [0.0096,0.0792317649986],
    [0.0097,0.0753668008523],
    [0.0098,0.0716903715424],
    [0.0099,0.0681932802477],
    [0.01,0.0648667787722],
    [0.0101,0.0648667787722],
]
# transpose
pscad_dat_voltage_0p1 = map(list, zip(*pscad_dat_voltage_0p1))

# concatenate lists together
pscad_dat_current = [pscad_dat_current_0p1, pscad_dat_current_0p8]
pscad_dat_voltage = [pscad_dat_voltage_0p1, pscad_dat_voltage_0p8]

# Main function
def main(args):
    # Circuit parameters
    i_0 = 0
    v_0 = 10
    inductance_L = 0.02 # 20 mH
    resistance_R = 10 # 10 Ohms
    voltage_E = 10 # 10 Volts
    simulation_time = 0.01 # 10 ms
    time_deltas = [0.0001, 0.0008] # 100 us, 800 us time deltas

    print("Simulation time = %g" % simulation_time)
    print("Will run for the following time deltas:")
    for delta in time_deltas:
        print("Time delta = %g" % delta)

    cont_results = [[],[],[]]
    trap_results = [[[],[],[]]]
    back_results = [[[],[],[]]]
    # indices for the results tuple
    TIME_INDEX = 0
    VALUE_INDEX = 1
    CURRENT_INDEX = 1
    VOLTAGE_INDEX = 2

    # Calculate results using continuous solution
    print("Simulating using continuous solution")
    # Calculate 3000 points for the continuous solution
    for step in range(3000):
        time = 0.0 + (step * simulation_time/3000)
        # calculate continuous solution
        # current
        i_continuous_next = (10/resistance_R) - ((10/resistance_R)*math.exp(-(resistance_R/inductance_L)*time))
        # inductor voltage
        v_continuous_next = 10 * math.exp((-(resistance_R)/(inductance_L))*time)
        # append the result
        cont_results[TIME_INDEX].append(0+(step*simulation_time/3000))
        cont_results[CURRENT_INDEX].append(i_continuous_next)
        cont_results[VOLTAGE_INDEX].append(v_continuous_next)

    # Calculate results using trapezoidal discretization
    print("Simulating using trapezoidal discretization")
    for delta_index, delta in enumerate(time_deltas):
        sim_steps = int(simulation_time/delta)
        for step in range(sim_steps):
            if step == 0:
                # Use initial conditions
                i_prev = i_0
                v_prev = v_0
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
    for delta_index, delta in enumerate(time_deltas):
        sim_steps = int(simulation_time/delta)
        for step in range(sim_steps):
            if step == 0:
                # Use initial conditions
                i_prev = i_0
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
    for delta_index, delta in enumerate(time_deltas):
        # Current plot
        fig, ax = plt.subplots(2)
        ax[0].plot(cont_results[TIME_INDEX], cont_results[CURRENT_INDEX], linestyle='dotted', label='Continuous')
        ax[0].plot(trap_results[delta_index][TIME_INDEX], trap_results[delta_index][CURRENT_INDEX], label='Trapezoidal')
        ax[0].plot(back_results[delta_index][TIME_INDEX], back_results[delta_index][CURRENT_INDEX], label='Backward Euler    ')
        ax[0].plot(pscad_dat_current[delta_index][TIME_INDEX], pscad_dat_current[delta_index][VALUE_INDEX], label='PSCAD')
        ax[0].set(xlabel='time (s)', ylabel='current (A)', title='Evaluations of $i(t)$ at $%gs$' % delta)
        ax[0].legend(loc='best', fancybox=True, shadow=True)
        ax[0].grid()
        # Voltage plot
        ax[1].plot(cont_results[TIME_INDEX], cont_results[VOLTAGE_INDEX], linestyle='dotted', label='Continuous')
        ax[1].plot(trap_results[delta_index][TIME_INDEX], trap_results[delta_index][VOLTAGE_INDEX], label='Trapezoidal')
        ax[1].plot(back_results[delta_index][TIME_INDEX], back_results[delta_index][VOLTAGE_INDEX], label='Backward Euler    ')
        ax[1].plot(pscad_dat_voltage[delta_index][TIME_INDEX], pscad_dat_voltage[delta_index][VALUE_INDEX], label='PSCAD')
        ax[1].set(xlabel='time (s)', ylabel='voltage (V)', title='Evaluations of $v_L(t)$ at $%gs$' % delta)
        ax[1].legend(loc='best', fancybox=True, shadow=True)
        ax[1].grid()
        fig.set_size_inches(6.5,8)
        fig.tight_layout()
        # plt.show()
        fig.savefig('compare_plot_' + str(delta).replace('.','p') + '.pgf')
        fig.savefig('compare_plot_' + str(delta).replace('.','p') + '.png')

if __name__ == '__main__':
    # the following sets up the argument parser for the program
    parser = argparse.ArgumentParser(description='Assignment 1a solution generator')

    args = parser.parse_args()

    main(args)