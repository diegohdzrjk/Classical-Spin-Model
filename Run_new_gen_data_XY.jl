# Version 1.0.5 de Julia
push!(LOAD_PATH, ".");
include("Lib_new_gen_data_Ising.jl")

# Hamiltonian to simulate is given by
# Hamiltonian = - H Sum_i Si cos(phi_i) - J Sum_ij Si Sj cos(phi_i-phi_j)

#    Variables of the function simulation
#        L           -> Grid size
#        Spin        -> Max spin in the system
#        Parameters  -> [ J, H ]
#        T_array     -> [ T_start , T_step , T_end ]
#        MC_array    -> [ Total MC steps, MC steps to save data,
#                         MC steps to save the correlation length]
#        doc_path    -> Folder where to save the data
#        flags       -> [Flag for correlation, Flag for save configuration,
#                        Flag apply threshold,flags for the helicity and the vorticity]

Spin = 1
MC_save = 500
MC_array = [5000, Int(MC_save), Int(MC_save/20.0)]
T_array = [0.1, 0.01, 1.0]

J, H = 1.0, 0.0
Parameters = [J, H]
flags = [false, true, true, false]
initial_path = "./"
repeat_simulation = 50

for L in [8, 16, 32, 64, 128]
    for j in 1:repeat_simulation
        try mkdir(initial_path*"DATA") catch; end
        try mkdir(initial_path*"DATA/L"*string(L)) catch; end

        E_threshold = [0.002*L*L/4.0, MC_array[1]/200]

        run_data = string("DATA/L",L,"/XY_S",Spin,"_J",J,"_H",H,"_MC",MC_array[1],"_Ti",T_array[1],"_Tf",T_array[3],"__",j)
        doc_path = initial_path*run_data
        # Creates the main folder
        try mkdir(doc_path) catch ; end

        message = string(Spin,"\n",L,"\n",T_array,"\n",MC_array,"\n",Parameters)
        data_save(doc_path*"/configuration_data.csv", "w", message)

        # Run the simulation
        simulacion(L, Spin, Parameters, T_array, MC_array, E_threshold, doc_path, flags)
    end
end
