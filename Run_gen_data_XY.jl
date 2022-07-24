push!(LOAD_PATH, ".");
using Lib_gen_data_XY

# Hamiltonian to simulate is given by
# Hamiltonian = - H Sum_i Si cos(phi_i) - J Sum_ij Si Sj cos(phi_i-phi_j) + K Sum_ij Si^2 Sj^2 cos(phi_i-phi_j)^2

#    Variables of the function simulation
#        L           -> Grid size
#        Spin        -> Max spin in the system
#        Parameters  -> [ J, H, K ]
#        T_array     -> [ T_start , T_step , T_end ]
#        MC_array    -> [ Total MC steps, MC steps to save data, MC steps to save the correlation length]
#        doc_path    -> Folder where to save the data
#        flags       -> [Flag for correlation, Flag for the helicty and vortcity, Mean phase as a func of MC]

L = 16
Spin = 0.5
T_array = [0.01, 0.02, 3.01]
MC_array = [20001, 20000, 200]

Parameters = [1.0,0.0,0.0]
flags = [true, true, false]

Param_var = -1.0:0.2:1.0
Param_idx = 3
Param_name = "K"

run_data = @sprintf("DATA/XY_BEG/XY_Square_S%.1f_L%d_J%.6f_H%.6f_K_%.f6_MC%d_T%.6f_%.6f", Spin, L, Parameters[1], Parameters[2], Parameters[3],
                                    MC_array[1], T_array[1], T_array[3])

doc_path = "/Users/diego/Documents/Modelo_Ising_XY/MultiParam_Ising_XY/"*run_data

# Creates the main folder
try
    mkdir(doc_path)
    println("The folder "*doc_path*" has been created.")
catch
    println("The folder "*doc_path*" existed already.")
end

if flags[1]
    c = "T"
elseif !flags[1]
    c = "F"
end

message = string(Spin,"\n",L,"\n",T_array,"\n",MC_array,"\n",Parameters,"\n",collect(Param_var),"\n",Param_name,"\n",Param_idx-1,"\n",c,"\nXY")
data_save(doc_path*"/data_run.csv", "w", message)

for j in Param_var
    Parameters[Param_idx] = j
    simulacion(L, Spin, Parameters, T_array, MC_array, doc_path, flags)
end
