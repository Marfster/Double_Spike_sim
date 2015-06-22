__author__ = 'marf'

from parameter import *
from dspike_formulas import *
import numpy as np
from bokeh.plotting import figure, output_file, show


spike1 = [["116"],["120", "117", "122"]]
spike2 = [["120"],["116", "117", "122"]]
# Define Model Parameter#
#*** Natural Fractionation ***#
fnat_sim = -1
#*** Instrumental Fractionation ***#
fins_sim =  2.2
#*** Sample/Spike Ratio ***#
mix = 0.5
mix_2 = [0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98]
mix_3 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
# Dampening of simulated mixtures
damp = 0.9

q_range = [x * 0.1 for x in range(1, 10)]
q_opti_117 = [0.4897]
spike_q = collections.OrderedDict()

for q in q_opti_117:
    spike_mix = {}
    for isotope in spike_117:
        spike_mix[isotope] = spike_117[isotope] * q + spike_122[isotope] * (1-q)
    spike_obj = load_abundance_dict(spike_mix)
    spike_q.update({q : spike_obj})
    print spike_q

print spike_q[0.4897].get_all_abundances()
    # Spike Calculation#

sim1 = calc_dspike_sample(Sn_meas_obj,df_new,spike_obj,Sn_mass_obj,spike1,"120")
sim1_mix_2 = sim1.spike_sim_p_range(mix_2, fnat_sim,fins_sim,damp,3,6,-0.1,-2,'z')
sim1_q = spike_sim_q_range(q_range,spike_117,spike_122,Sn_meas_obj,df_new,Sn_mass_obj,spike1,mix,fnat_sim,fins_sim,damp,3,6,-0.1,-2,'z')
sim2 = calc_dspike_sample(Sn_meas_obj,df_new,spike_obj,Sn_mass_obj,spike2,"120")
sim2_mix_2 = sim2.spike_sim_p_range(mix_2, fnat_sim,fins_sim,damp,3,6,-0.1,-2,'z')



x1, y1 = sim1.error_vs_p(sim1_mix_2, 'z')
x2, y2 = sim2.error_vs_p(sim2_mix_2, 'z')
x3, y3 = error_vs_q(sim1_q, 'z')
# prepare some data

# output to static HTML file
output_file(path+'/Spike_Sim/'+"Spike_Sim117-122_116_Qvar.html", title="Error plot Spike 117-122, 116")

# create a new plot with a title and axis labels
plot = figure(title="Error plot Spike 117-122, 116", x_axis_label='Prop. of 117Sn in 117-122Sn Double Spike', y_axis_label='1s mean Fnat',
           x_axis_type="log", x_range=[0.1, 1], y_range=[0, 5000])

# add a line renderer with legend and line thickness
plot.line(x3, y3, legend="117-122,116", line_width=2, color='red')
#p.line(x2, y2, legend="117-122,120", line_width=2, color='blue')
# show the results
show(plot)

#log_df.to_csv(path + "Sn117-122_120_1367.csv")

p_rang = np.linspace(0.1, 0.9, 18)
q_rang = np.linspace(0.1, 0.9, 18)
xx, yy = np.meshgrid(p_rang, q_rang)

z = collections.OrderedDict()
for p in p_rang:
    counter = 0
    sim1_q = spike_sim_q_range(q_rang,spike_117,spike_122,Sn_meas_obj,df_new,Sn_mass_obj,spike1,p,fnat_sim,fins_sim,damp,3,6,-0.1,-2,'z')
    q_ls, z_ls = error_vs_q(sim1_q, 'z')
    z.update({p: collections.OrderedDict()})
    for q_value in q_ls:
        z[p].update({q_value : z_ls[counter]})
        counter+=1
    print z[p]

z_df = pd.DataFrame.from_dict(z,"index")
print z_df
z_matrix = z_df.as_matrix()


pl.contourf(p_rang, q_rang, z_matrix, 8, alpha=.75, cmap=pl.cm.hot)
C = pl.contour(p_rang, q_rang, z_matrix, 8, colors='black', linewidth=.5)
pl.clabel(C, inline=1, fontsize=10)

pl.show()