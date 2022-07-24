__precompile__()

module Lib_gen_data_XY
    export simulacion, data_save

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# Hamiltonian to simulate is given by
# Hamiltonian = - H Sum_i cos(phi_i) - J Sum_ij cos(phi_i-phi_j)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
    function data_save(file::AbstractString, op_type::AbstractString, message::AbstractString)
        # Opening the file
        out = open(file, op_type)
        # Writing the labels of the data to store
        write(out, message)
        # Closing the file
        flush(out); close(out);
    end

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
    function create_data_folders(doc_path::AbstractString, Spin::Float64, L::Int64, Parameters::Array{Float64, 1},
                                T_array::Array{Float64,1}, MC_array::Array{Int64,1}, flags::Array{Bool, 1})
        """
                    L           -> Grid size
                    Spin        -> Max spin in the system
                    Parameters  -> [ J, H ]
                    T_array     -> [ T_start , T_step , T_end ]
                    MC_array    -> [ Total MC steps, MC steps to save data,
                                    MC steps to save the correlation length, MC to save the configuration]
                    doc_path    -> Folder where to save the data
                    flags       -> [Flag for correlation, flag for the helicity and vorticity,
                                    Flag for save configuration, Flag apply threshold]
        """
        # Name of the main folder where all the information is stored
        data_folder = doc_path*@sprintf("/XY_Square_S%.1f_L%d_J%.6f_H%.6f_K%.6f_MC%d_T%.6f_%.6f", Spin, L, Parameters[1], Parameters[2],
                                                                                        Parameters[3],MC_array[1], T_array[1], T_array[3])
        # If the folder is non existing it creates the folder
        try
            mkdir(data_folder)
            println("The folder "*data_folder*" has been created.")
        catch
            println("The folder "*data_folder*" existed already.")
        end

        # Name of the files where the information is stored
        File_logs = data_folder*"/logs.csv"
        File_data = data_folder*"/data.csv"
        File_correlation = data_folder*"/correlation.csv"
        Folder_matix_config = data_folder*"/Matrix_configurations"

        # Register the initial time of the program
        data_save(File_logs, "w", @sprintf("Start time %s \n", string(now())))
        # Register the labels of the data
        labels = "Temperature, Energy, Magnetization, Cv, Chi_m, Vortex_plus, Vortex_minus, Upsilon_x, Upsilon_y\n"
        data_save(File_data, "w", labels)
        data_save(File_phase,"w", "")

        # If the correlation length is to be calculated, it creates the file where the data is stored
        if flags[1]
            r_s = ""
            for r in 0:div(L,2)
                r_s *= @sprintf(",%d", r)
            end
            # Register the labels of the the correlation
            data_save(File_correlation, "w", "T"*r_s*"\n")
        end

        if flags[3]
            try
                mkdir(Folder_matix_config)
                println("The folder "*Folder_matix_config*" has been created.")

            catch
                println("The folder "*Folder_matix_config*" existed already.")
            end
        end

        # Resturns the name of the files
        return [File_logs, File_data, File_correlation, Folder_matix_config]
    end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
    function simulacion(L::Int64, Spin::Float64, Parameters::Array{Float64, 1}, T_array::Array{Float64,1},
         MC_array::Array{Int64,1}, E_threshold::Array{Float64, 1}, doc_path::AbstractString, flags::Array{Bool, 1})
        """
            L           -> Grid size
            Spin        -> Max spin in the system
            Parameters  -> [ J, H ]
            T_array     -> [ T_start, T_step, T_end]
            MC_array    -> [ Total MC steps, MC steps to save data,
                            MC steps to save the correlation length, MC to save the configuration]
            E_threshold -> Energy threshold for stopping the simulation, MC to check
            doc_path    -> Folder where to save the data
            flags       -> [Flag for correlation, flags for the helicity and the vorticity,
                            Flag for save configuration, Flag apply threshold]
        """

        # Temperatures array and posible spin array
        Temperatures = collect(T_array[1]:T_array[2]:T_array[3])
        Spins = [1.0]#collect(0:1:Spin)

        # Some information about the the to run
        println("Temperatures: -> ",T_array[1]," - ", T_array[3]," in steps of ", T_array[2])
        println("Total number of temperatures ", length(Temperatures))
        println("Parameters: J -> ",Parameters[1],"  H -> ",  Parameters[2],"  K -> ",  Parameters[3])
        println("Grid size ",L,"x",L)
        println("Spins ", Spins)
        println("Montecarlo Steps -> ", MC_array[1])
        println("Energy threshold -> ", flags[4], "\t ", E_threshold[1], "\t Check each ->", E_threshold[2], "MC Steps")
        println("Calculate the correlation length -> ", flags[1], " Each ", MC_array[3], "MC")
        println("Calculate the Helicity and Vorticity -> ", flags[2])

        # Create the folder and the data
        # Files = [logs, data, correlation, Phases]
        Files = create_data_folders(doc_path, Spin, L, Parameters, T_array, MC_array, flags)

        # Evaluate the system mean properties for diferent temperatures
        @sync @parallel for t in Temperatures
            # Display the progress on
            @printf "Starting T = %.5f\n" t

            # Create a file to write the matrix configuration as a function of the MC steps
            if flags[3]
                data_save(Files[4]*@sprintf("/T%.5f.csv", t),"w", "")
            end

            # Iniciate the Montecarlo simulation
            tic()
            # Create a random spin matrix of size LxL
            # And a phase matrix
            S_ini = rand(Spins, L, L)
            Phase_ini = 2*pi*rand(L,L)
            My_data, Correl_array = configuration_run(L, Spins, S_ini, Phase_ini, t, Parameters, MC_array, E_threshold, Files, flags)
            tot_time = toq()

            # Save the data
            message = string(t)
            for i in 1:length(My_data)
                message *= string(",", round(My_data[i], 12))
            end
            data_save(Files[2],"a", message*"\n")

            if flags[1]
                # Save the correlation length
                corr_str = @sprintf("%.5f", t)
                for r in 1:div(L,2)+1
                    corr_str *= @sprintf(",%.5f", Correl_array[r])
                end
                data_save(Files[3],"a", corr_str*"\n")
            end
            # Log the time of computation
            data_save(Files[1], "a", @sprintf("T=%.5f time=%.4f\n", t, tot_time))
        end

        # Log the final time of the simulation
        data_save(Files[1], "a", @sprintf("Final time %s \n", string(now())))

    end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
    function configuration_run(L::Int64, Spins::Array{Float64, 1}, S::Array{Float64, 2}, Phase::Array{Float64, 2}, T::Float64,
                            Parameters::Array{Float64, 1}, MC_array::Array{Int64,1}, E_threshold::Array{Float64, 1}, Files::Array{String,1},
                            flags::Array{Bool, 1})
        """
            L           -> Grid size
            Spins       -> List of spins in the system
            S           -> Spin matrix
            T           -> Temperature
            Parameters  -> [ J, H ]
            MC_array    -> [ Total MC steps, MC steps to save data,
                            MC steps to save the correlation length, MC to save the configuration]
            Files       -> Files
            E_threshold -> Energy threshold for stopping the simulation, MC to check
            flags       -> [Flag for correlation, flags for the helicity and the vorticity,
                            Flag for save configuration, Flag apply threshold]
            ------------------------------------------------------------------------------------------
            DATA = [ <E>, <E*E>, <M>, <M*M>]
        """

        # Define the variables that control the Montecarlo parameters
        N = L*L
        Evol_Steps = (MC_array[1]-MC_array[2])*N
        Save_steps = MC_array[2]*N

        E_sys_temp = - Parameters[2]*sum(cos(Phase)) - 0.5*Parameters[1]*sum(
                                    cos(Phase - circshift(Phase,(1,0)))+     #top
                                    cos(Phase - circshift(Phase,(-1,0)))+   #bottom
                                    cos(Phase - circshift(Phase,(0,1)))+     #right
                                    cos(Phase - circshift(Phase,(0,-1)))    #left
                                    )

        # Runs the montecarlo steps "until" equilibrium
        for step in 1:max(Evol_Steps,1)
            # Check if the simulation converged by checking the threshold
            if flags[3]
                if step%(N*E_threshold[2]) == 0
                    # Calculate the energy of the system
                    E_sys_new = - Parameters[2]*sum(cos(Phase)) - 0.5*Parameters[1]*sum(
                                                cos(Phase - circshift(Phase,(1,0)))+     #top
                                                cos(Phase - circshift(Phase,(-1,0)))+   #bottom
                                                cos(Phase - circshift(Phase,(0,1)))+     #right
                                                cos(Phase - circshift(Phase,(0,-1)))    #left
                                                )

                    if abs(E_sys_temp-E_sys_new) < E_threshold[1]
                        print("The method converged before, energy threshold apply. MC -> ", step/N)
                        println("\t Temperature -> ", T)
                        break
                    end
                    E_sys_temp = E_sys_new
                end
            end

            # Choose a random spin and phase of the grid and a new spin proyection
            i = rand(1:N)
            ph = 2*pi*rand()

            # For clarity we use these extra variables
            i_top = div(i-1,L)*L+mod(i-2,L)+1
            i_bottom = div(i-1,L)*L+mod(i,L)+1
            i_right = mod(i-1+L,N)+1
            i_left = mod(i-1-L,N)+1

            Phase_i = Phase[i]
            Phase_neig = Phase[[i_top, i_bottom, i_right, i_left]]

            # Calculate the Diference in energy
            DeltaE  =   Parameters[2]*( cos(Phase[i])-cos(ph)) +
                        Parameters[1]*( cos(Phase_i-Phase_neig[1])-cos(ph-Phase_neig[1])+
                                        cos(Phase_i-Phase_neig[2])-cos(ph-Phase_neig[2])+
                                        cos(Phase_i-Phase_neig[3])-cos(ph-Phase_neig[3])+
                                        cos(Phase_i-Phase_neig[4])-cos(ph-Phase_neig[4])
                                        )
            # If the state is accepted the spin and phase matrix gets updated
            if rand() < min(1.0, exp(-DeltaE/T))
                Phase[i] = ph
            end

        end

        # Up to this point the properties of the system are been calculated
        # The next part is in charge of calculating the properties of the system meanwhile
        # the system continues to evolve
        # Create the thermodynamic arrays of variables
        DATA = zeros(Complex64, 4)
        Corr_array = zeros(Float64, div(L,2)+1)
        # Create an array to count the number of vortex and antivortex
        VORTEX = zeros(Float64, 2)
        # Create the array for the helicity modulus in both x and y directions
        HELICITY_MOD = zeros(Float64, 2)

        # Calculate the Magnetization
        Mag = sum(cos(Phase)+im*sin(Phase))
        # Calculate the energy of the system
        E_sys = - Parameters[2]*sum(cos(Phase)) - 0.5*Parameters[1]*sum(
                                    cos(Phase - circshift(Phase,(1,0)))+     #top
                                    cos(Phase - circshift(Phase,(-1,0)))+   #bottom
                                    cos(Phase - circshift(Phase,(0,1)))+     #right
                                    cos(Phase - circshift(Phase,(0,-1)))    #left
                                    )

        if flags[2]
            # Calculate the vorticity of the phase matrix
            Vortex_matrix = round((rem(Phase - circshift(Phase,(0,1)), pi) +
                                    rem(circshift(Phase,(0,1)) - circshift(Phase,(1,1)), pi) +
                                    rem(circshift(Phase,(1,1)) - circshift(Phase,(1,0)), pi) +
                                    rem(circshift(Phase,(1,0)) - Phase, pi)
                                    )/pi)
        end

        # Runs the montecarlo steps taking the averages
        for step in 1:Save_steps
            # Save the data in the array
            DATA += [E_sys, E_sys*E_sys, Mag, Mag*Mag']

            if step%(500*N)==1 && flags[3]
                message = @sprintf("%.5f,%.5f,%.5f,%.5f\n", mean(Phase), var(Phase), norm(Mag)/N, E_sys/N)
                data_save(Files[5]*@sprintf("/T%.5f.csv", T), "a", message)

                data_save(Files[6]*@sprintf("/T%.5f.csv", T), "a", string(reshape(Phase,1, N))[2:end-1]*"\n")
            end

            if flags[2]
                # Count the number of vortex with carge +1 and -1
                VORTEX += [N-countnz(Vortex_matrix-1), N-countnz(Vortex_matrix+1)]

                # Calculate the helicity modulus in both x and y directions
                Upsilon_x = sin(Phase - circshift(Phase,(0,1))) - sin(Phase - circshift(Phase,(0,-1)))
                Upsilon_y = sin(Phase - circshift(Phase,(1,0))) - sin(Phase - circshift(Phase,(-1,0)))
                HELICITY_MOD += [sum(0.5*Upsilon_x)^2, sum(0.5*Upsilon_y)^2]
            end
            # Choose a random spin of the grid and a new spin proyection
            i = rand(1:N)
            ph = 2*pi*rand()

            # For clarity we use these extra variables
            i_top = div(i-1,L)*L+mod(i-2,L)+1
            i_bottom = div(i-1,L)*L+mod(i,L)+1
            i_right = mod(i-1+L,N)+1
            i_left = mod(i-1-L,N)+1
            i_top_right = div(i_right-1,L)*L+mod(i_right-2,L)+1

            Phase_i = Phase[i]
            Phase_neig = Phase[[i_top, i_bottom, i_right, i_left, i_top_right]]

            # Calculate the Diference in energy
            DeltaE = Parameters[2]*(cos(Phase_i)-cos(ph))+Parameters[1]*(
                                    cos(Phase_i-Phase_neig[1])-cos(ph-Phase_neig[1])+
                                    cos(Phase_i-Phase_neig[2])-cos(ph-Phase_neig[2])+
                                    cos(Phase_i-Phase_neig[3])-cos(ph-Phase_neig[3])+
                                    cos(Phase_i-Phase_neig[4])-cos(ph-Phase_neig[4])
                                    )


            # If the state is accepted the spin and phase matrix, the energy, the magnetization,
            # and the concentration of spin 0 gets updated
            if rand() < min(1.0, exp(-DeltaE/T))
                E_sys += DeltaE
                Mag += cos(ph)+im*sin(ph)-(cos(Phase[i])+im*sin(Phase[i]))
                Phase[i] = ph
                if flags[2]
                    # Calculate the vorticity
                    Vortex_matrix[i] = round((rem(Phase_i-Phase_neig[3], pi)+
                                              rem(Phase_neig[3]-Phase_neig[5], pi)+
                                              rem(Phase_neig[5]-Phase_neig[1], pi)+
                                              rem(Phase_neig[1]-Phase_i, pi)
                                            )/pi)
                end
            end

            # The correlation length is calculated
            if flags[1] && step%(N*MC_array[3])==1
                for r in 0:div(L,2)
                    Corr =  cos(Phase - circshift(Phase,(0,r))) +
                            cos(Phase - circshift(Phase,(0,-r))) +
                            cos(Phase - circshift(Phase,(r,0))) +
                            cos(Phase - circshift(Phase,(-r,0)))
                    Corr_array[r+1] += 0.25*mean(Corr)
                end
            end

            if flags[3] && step%(N*MC_array[4])==0
                data_save(Files[4]*@sprintf("/T%.5f.csv", T), "a", string(reshape(S,1, N))[2:end-1]*"\n")
            end

        end

        # Take averages
        DATA = DATA/Save_steps

        if flags[2]
            VORTEX = VORTEX/Save_steps
            HELICITY_MOD = HELICITY_MOD/Save_steps
            HELICITY_MOD = -0.5*DATA[1] - HELICITY_MOD/T
        end

        Corr_array = Corr_array*MC_array[3]/MC_array[2]

        # Calculate the suceptibilities and the helicty
        Cv = (DATA[2] - DATA[1]*DATA[1])/(T*T)
        Chim = (DATA[4] - norm(DATA[3])^2)/T

        return real([DATA[1], norm(DATA[3]), Cv, Chim, VORTEX[1], VORTEX[2], HELICITY_MOD[1], HELICITY_MOD[2]])/N, Corr_array

    end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
end
