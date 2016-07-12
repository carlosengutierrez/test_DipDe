# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 13:22:48 2016

@author: carlosengutierrez
"""

import nest
import nest.raster_plot
import numpy as np
from dipde.internals.internalpopulation import InternalPopulation
from dipde.internals.externalpopulation import ExternalPopulation
from dipde.internals.simulation import Simulation
from dipde.internals.connection import Connection as Connection

# Spiking Network Model Parameters 
N_I = 200
N_E = 800
N_rec = 50  # recorded neurons
u_rest= 0.#-70. # resting potential
u_reset = 0.#-65.#10. #reset potential
u_th= 20.#-55.#20.  # firing threshold
t_m = 20. # time constante(ms)
C_m = 1. #Capacitance pF
R = t_m/C_m  # tau_m = R C , Resistance
delay = 1.5 #synapse delay
t_ref = 0. #2. #0.#2. #absolute refractory period
dt = 0.1
simtime = 1000. # simulation time
I_ext = 0.#External input
g = 5. # ratio inhibitory weight/excitatory weight
epsilon = 0.1  # connection probability
C_E = int(epsilon*N_E)  # number of excitatory synapses per neuron
C_I = int(epsilon*N_I)   # number of inhibitory synapses per neuron
neuron_params= {"C_m":        C_m, # dictionary with neuron parameters
                "tau_m":      t_m,
                "t_ref":      t_ref,
                "E_L":        u_rest,
                "V_reset":    u_reset,
                "V_m":        u_rest,
                "V_th":       u_th,
                "I_e":I_ext}
J     = 0.1     # postsynaptic amplitude in mV
J_ex  = J       # amplitude of excitatory postsynaptic potential
J_in  = -g*J_ex # amplitude of inhibitory postsynaptic potential
nu_th  = u_th/(J*C_E*t_m)
eta     = 2.0  # external rate relative to threshold rate
nu_ex  = eta*nu_th
p_rate = 1000.0*nu_ex*C_E # external rate [Hz] to add to poisson generator.

#### Nest Simulation
nest.ResetKernel() # reset NEST
nest.SetKernelStatus({'resolution': dt, "print_time": True, "overwrite_files": True }) # set resolution
nest.SetDefaults("iaf_psc_delta", neuron_params)
nest.SetDefaults("poisson_generator",{"rate": p_rate})

nodes_E = nest.Create("iaf_psc_delta",N_E)
nodes_I = nest.Create("iaf_psc_delta",N_I)
noise    = nest.Create("poisson_generator")
spikes_E  = nest.Create("spike_detector")
spikes_I  = nest.Create("spike_detector")

nest.SetStatus(spikes_E,[{"label": "Exc",
                         "withtime": True,
                         "withgid": True,
                         "to_file": True}])
nest.SetStatus(spikes_I,[{"label": "Inh",
                         "withtime": True,
                         "withgid": True,
                         "to_file": True}])
nest.CopyModel("static_synapse","excitatory",{"weight":J_ex, "delay":delay})
nest.CopyModel("static_synapse","inhibitory",{"weight":J_in, "delay":delay})

nest.Connect(noise,nodes_E, syn_spec="excitatory")
nest.Connect(noise,nodes_I, syn_spec="excitatory")
nest.Connect(nodes_E[:N_rec], spikes_E)#, syn_spec="excitatory") #
nest.Connect(nodes_I[:N_rec], spikes_I)#, syn_spec="excitatory") #

conn_params_ex = {'rule': 'fixed_indegree', 'indegree': C_E}
nest.Connect(nodes_E, nodes_E+nodes_I, conn_params_ex, "excitatory")
conn_params_in = {'rule': 'fixed_indegree', 'indegree': C_I}
nest.Connect(nodes_I, nodes_E+nodes_I, conn_params_in, "inhibitory")

nest.Simulate(simtime) #simulate
                                 
events_ex = nest.GetStatus(spikes_E,"n_events")[0]
events_in = nest.GetStatus(spikes_I,"n_events")[0]
rate_ex   = events_ex/simtime*1000./N_rec
rate_in   = events_in/simtime*1000./N_rec
print 'NEST   Mean Exc rate: ',rate_ex,'[Hz]'
print 'NEST   Mean Inh rate: ',rate_in,'[Hz]'


###DipDe 
# Settings:
t0 = 0.
dt = 0.0001
dv = 0.1
tf = 1.#.1 #simulation duration.
update_method = 'approx'
approx_order = 1
tol = 1e-14
verbose = False

# Create simulation:
b1 = ExternalPopulation(p_rate, record=True)
i1 = InternalPopulation(tau_m=t_m*1e-3,v_min=-0.02, v_max=u_th*1e-3,
                        dv=dv, update_method=update_method, approx_order=approx_order,
                        tol=tol,record=True, curr_firing_rate=0.0,norm='inf')
i2 = InternalPopulation(tau_m=t_m*1e-3,v_min=-0.02, v_max=u_th*1e-3,
                        dv=dv, update_method=update_method, approx_order=approx_order,
                        tol=tol,record=True, curr_firing_rate=0.0,norm='inf')

b1_i1 = Connection(b1, i1, C_E, weights=[J_ex*1e-3], probs=[1.], delay=delay*1e-3)
i1_i1 = Connection(i1, i1,C_E , weights=[J_ex*1e-3], probs=[1.], delay=delay*1e-3)

b1_i2 = Connection(b1, i2, C_E, weights=[J_ex*1e-3], probs=[1.], delay=delay*1e-3)
i2_i2 = Connection(i2, i2, C_I, weights=[J_in*1e-3], probs=[1.], delay=delay*1e-3)

i1_i2 = Connection(i1, i2, C_E, weights=[J_ex*1e-3], probs=[1.], delay=delay*1e-3)
i2_i1 = Connection(i2, i1, C_I, weights=[J_in*1e-3], probs=[1.], delay=delay*1e-3)


simulation = Simulation([b1, i1,i2], [b1_i1, i1_i1,b1_i2,i2_i2,i1_i2,i2_i1], verbose=verbose)
simulation.run(dt=dt, tf=tf, t0=t0)

i1 = simulation.population_list[1]
print 'DipDe  Mean Exc rate: ',np.array(i1.firing_rate_record).mean(),'[Hz]'
i2 = simulation.population_list[2]
print 'DipDe  Mean Inh rate: ',np.array(i2.firing_rate_record).mean(),'[Hz]'
