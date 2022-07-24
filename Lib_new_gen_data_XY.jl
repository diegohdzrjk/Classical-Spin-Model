# Version 1.0.5 de Julia
__precompile__()
module Lib_gen_data_Ising
    using Dates
    using Distributed
    export simulacion, data_save
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# Hamiltonian to simulate is given by
# Hamiltonian = - H Sum_i Si cos(phi_i) - J Sum_ij Si Sj cos(phi_i-phi_j)
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
    @everywhere function data_save(file::AbstractString, op_type::AbstractString, message::AbstractString)
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
    @everywhere function create_data_folders(doc_path::AbstractString, Spin::Float64, L::Int64, Parameters::Array{Float64, 1},
                                T_array::Array{Float64,1}, MC_array::Array{Int64,1}, flags::Array{Bool, 1})
        """
                    L           -> Grid size
                    Spin        -> Max spin in the system
                    Parameters  -> [ J, H ]
                    T_array     -> [ T_start , T_step , T_end ]
                    MC_array    -> [ Total MC steps, MC steps to save data,
                                    MC steps to save the correlation length]
                    doc_path -> Folder where to save the data
                    flags       -> [Flag for correlation, Flag for save configuration, Flag apply threshold]
        """
        # Name of the main folder where all the information is stored
        data_folder = doc_path*"/DATA"
        # If the folder is non existing it creates the folder
        try mkdir(data_folder) catch ; end

        # Name of the files where the information is stored
        File_logs = data_folder*"/logs.csv"
        File_data = data_folder*"/data.csv"
        File_correlation = data_folder*"/correlation.csv"
        Folder_matix_config = data_folder*"/Matrix_configurations"

        # Register the initial printf time of the program
        data_save(File_logs, "w", string("Start time %s \n", time()))
        # Register the labels of the data
        # T, E, Mx, My, Cv, Chi_m, Vx, Vy, Hx, Hy
        data_save(File_data, "w", "Temperature,Energy,Mag_x,Mag_y,Cv,Chi_m,Vort_x,Vort_y,Hel_x,Hel_y\n")

        # If the correlation length is to be calculated, it creates the file where the data is stored
        if flags[1]
            r_s = ""
            for r in 0:div(L,2)
                r_s *= string(",", r)
            end
            # Register the labels of the the correlation
            data_save(File_correlation, "w", "T"*r_s*"\n")
        end

        if flags[2]
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
    @everywhere function simulacion(L::Int64, Spin::Float64, Parameters::Array{Float64, 1}, T_array::Array{Float64,1},
         MC_array::Array{Int64,1}, E_threshold::Array{Float64, 1}, doc_path::AbstractString, flags::Array{Bool, 1}=[false, false, false])
        """
            L           -> Grid size
            Spin        -> Max spin in the system
            Parameters  -> [ J, H ]
            T_array     -> [ T_start, T_step, T_end]
            MC_array    -> [ Total MC steps, MC steps to save data,
                            MC steps to save the correlation length]
            E_threshold -> Energy threshold for stopping the simulation, MC to check
            doc_path    -> Folder where to save the data
            flags       -> [Flag for correlation, Flag for save configuration
                            Flag apply threshold,flags for the helicity and the vorticity]
        """

        # Temperatures array and posible spin array
        Temperatures = collect(T_array[1]:T_array[2]:T_array[3])
        Spins = [Spin]#collect(0:1:Spin)

        # Some information about the the to run
        println("Temperatures: -> ",T_array[1]," - ", T_array[3]," in steps of ", T_array[2])
        println("Total number of temperatures ", length(Temperatures))
        println("Parameters: J -> ",Parameters[1],"  H -> ",  Parameters[2])
        println("Grid size ",L,"x",L)
        println("Spins ", Spins)
        println("Montecarlo Steps -> ", MC_array[1])
        println("Energy threshold -> ", flags[3], "\t ", E_threshold[1], "\t Check each ->", E_threshold[2], "MC Steps")
        println("Calculate the correlation length -> ", flags[1])
        println("Calculate the Helicity and Vorticity -> ", flags[4])

        # Create the folder and the data
        # Files = [logs, data, correlation]
        Files = create_data_folders(doc_path, Spin, L, Parameters, T_array, MC_array, flags)

        # Create a random spin matrix of size LxL
        S_ini = rand(Spins, L, L)

        # Evaluate the system mean properties for diferent temperatures
        @sync @distributed for t in Temperatures
            # Display the progress on
            println(string("Starting T = ",t,"\n"))

            # Create a file to write the matrix configuration as a function of the MC steps
            if flags[2]
                data_save(Files[4]*string("/T", t,".csv"),"w", "")
            end

            # Iniciate the Montecarlo simulation
            start = time()
            S = copy(S_ini)
            Phase_ini = 2*pi*rand(L,L)

            # Energy per site
            # Magnetization x per site
            # Magnetization y per site
            # Cv per site
            # Chim per site
            # Vortex x per site
            # Vortex y per site
            # Helicity x per site
            # Helicity y per site
            # Correlation length
            E, Mx, My, Cv, Chi_m, Vx, Vy, Hx, Hy, Correl_array = configuration_run(L, Spins, S, Phase_ini, t, Parameters, MC_array, E_threshold, Files, flags)
            tot_time = time() - start

            # Save the data
            m = string(t)
            for var in [E, Mx, My, Cv, Chi_m, Vx, Vy, Hx, Hy]
                m *= string(",",var)
            end
            m *= "\n"
            data_save(Files[2],"a", m)

            if flags[1]
                # Save the correlation length
                corr_str = string(t)
                for r in 1:div(L,2)+1
                    corr_str *= string(",", Correl_array[r])
                end
                data_save(Files[3],"a", corr_str*"\n")
            end
            # Log the time of computation
            data_save(Files[1], "a", string("T=",t," time=",tot_time,"\n"))
        end

        # Log the final time of the simulation
        data_save(Files[1], "a", string("Final time ",time()," \n"))

    end
    #--------------------------------------------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------------------------------------------
    @everywhere function configuration_run(L::Int64, Spins::Array{Float64, 1}, S::Array{Float64, 2}, Phase::Array{Float64, 2},T::Float64,
                            Parameters::Array{Float64, 1}, MC_array::Array{Int64,1}, E_threshold::Array{Float64, 1}, Files::Array{String,1},
                            flags::Array{Bool, 1})
        """
            L           -> Grid size
            Spins       -> List of spins in the system
            S           -> Spin matrix
            Phase       -> Phase matrix
            T           -> Temperature
            Parameters  -> [ J, H ]
            MC_array    -> [ Total MC steps, MC steps to save data,
                            MC steps to save the correlation length, MC to save the configuration]
            E_threshold -> Energy threshold for stopping the simulation, MC to check
            Files       -> Files
            flags       -> [Flag for correlation, Flag for save configuration]
            ------------------------------------------------------------------------------------------
            DATA = [ <E>, <E*E>, <M>, <M*M> ]
        """

        # Define the variables that control the Montecarlo parameters
        N = L*L
        Evol_Steps = (MC_array[1]-MC_array[2])*N
        Save_steps = MC_array[2]*N

        # Define the variables that control the exchange of rows and columns
        Step_change = Int(round(0.01*N*E_threshold[2]))
        n_change = Int(round(0.5*L))

        # Calculate the Magnetization
        # M = Mx + i My
        Mag = sum(S.*(cos(Phase)+im*sin(Phase)))
        # Calculate the energy of the system
        E_sys_temp = - Parameters[2]*sum(S.*cos(Phase)) - 0.5*Parameters[1]*sum(
                                    S.*circshift(S,(1,0)).*cos(Phase - circshift(Phase,(1,0)))+     #top
                                    S.*circshift(S,(-1,0)).*cos(Phase - circshift(Phase,(-1,0)))+   #bottom
                                    S.*circshift(S,(0,1)).*cos(Phase - circshift(Phase,(0,1)))+     #right
                                    S.*circshift(S,(0,-1)).*cos(Phase - circshift(Phase,(0,-1))))    #left
                                    #)

        # Runs the montecarlo steps "until" equilibrium
        for step in 1:Evol_Steps
            # Check if the simulation converged by checking the threshold
            if flags[3]
                if step%(N*E_threshold[2]) == 0
                    Mag = sum(S.*(cos(Phase)+im*sin(Phase)))
                    # Calculate the energy of the system
                    E_sys_new = - Parameters[2]*sum(S.*cos(Phase)) - 0.5*Parameters[1]*sum(
                                                S.*circshift(S,(1,0)).*cos(Phase - circshift(Phase,(1,0)))+     #top
                                                S.*circshift(S,(-1,0)).*cos(Phase - circshift(Phase,(-1,0)))+   #bottom
                                                S.*circshift(S,(0,1)).*cos(Phase - circshift(Phase,(0,1)))+     #right
                                                S.*circshift(S,(0,-1)).*cos(Phase - circshift(Phase,(0,-1)))    #left
                                                )

                    if abs(E_sys_temp-E_sys_new) < E_threshold[1]
                        print("The method converged before, energy threshold apply. MC -> ", step/N)
                        println("\t Temperature -> ", T)
                        break
                    end
                    E_sys_temp = E_sys_new
                end
            end

            # Choose a random spin of the grid and a new spin proyection
            i = rand(1:N)
            sp = rand(Spins)
            ph = 2*pi*rand()

            # For clarity we use these extra variables
            i_top = div(i-1,L)*L+mod(i-2,L)+1
            i_bottom = div(i-1,L)*L+mod(i,L)+1
            i_right = mod(i-1+L,N)+1
            i_left = mod(i-1-L,N)+1

            Phase_i = Phase[i]
            s_i = S[i]
            Phase_neig = Phase[[i_top, i_bottom, i_right, i_left]]
            s_neig = S[[i_top, i_bottom, i_right, i_left]]

            # Calculate the Diference in energy
            DeltaE  =   Parameters[2]*( s_i*cos(Phase_i)-s*cos(ph)) +
                        Parameters[1]*( s_i*s_neig[1]*cos(Phase_i-Phase_neig[1])-sp*s_neig[1]*cos(ph-Phase_neig[1])+
                                        s_i*s_neig[2]*cos(Phase_i-Phase_neig[2])-sp*s_neig[2]*cos(ph-Phase_neig[2])+
                                        s_i*s_neig[3]*cos(Phase_i-Phase_neig[3])-sp*s_neig[3]*cos(ph-Phase_neig[3])+
                                        s_i*s_neig[4]*cos(Phase_i-Phase_neig[4])-sp*s_neig[4]*cos(ph-Phase_neig[4])
                                        )

            # If the state is accepted the spin matrix gets updated
            if rand() < min(1.0, exp(-DeltaE/T))
                Phase[i] = ph
                S[i] = sp
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
        Mag = sum(S.*(cos(Phase)+im*sin(Phase)))

        # Calculate the energy of the system
        E_sys = - Parameters[2]*sum(S.*cos(Phase)) - 0.5*Parameters[1]*sum(
                                    S.*circshift(S,(1,0)).*cos(Phase - circshift(Phase,(1,0)))+     #top
                                    S.*circshift(S,(-1,0)).*cos(Phase - circshift(Phase,(-1,0)))+   #bottom
                                    S.*circshift(S,(0,1)).*cos(Phase - circshift(Phase,(0,1)))+     #right
                                    S.*circshift(S,(0,-1)).*cos(Phase - circshift(Phase,(0,-1)))    #left
                                    )

        if flags[4]
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

            if flags[4]
                # Count the number of vortex with carge +1 and -1
                VORTEX += [N-countnz(Vortex_matrix-1), N-countnz(Vortex_matrix+1)]

                # Calculate the helicity modulus in both x and y directions
                Upsilon_x = sin(Phase - circshift(Phase,(0,1))) - sin(Phase - circshift(Phase,(0,-1)))
                Upsilon_y = sin(Phase - circshift(Phase,(1,0))) - sin(Phase - circshift(Phase,(-1,0)))
                HELICITY_MOD += [sum(0.5*Upsilon_x)^2, sum(0.5*Upsilon_y)^2]
            end

            # Choose a random spin of the grid and a new spin proyection
            i = rand(1:N)
            sp = rand(Spins)
            ph = 2*pi*rand()

            # For clarity we use these extra variables
            i_top = div(i-1,L)*L+mod(i-2,L)+1
            i_bottom = div(i-1,L)*L+mod(i,L)+1
            i_right = mod(i-1+L,N)+1
            i_left = mod(i-1-L,N)+1

            Phase_i = Phase[i]
            s_i = S[i]
            Phase_neig = Phase[[i_top, i_bottom, i_right, i_left]]
            s_neig = S[[i_top, i_bottom, i_right, i_left]]

            # Calculate the Diference in energy
            DeltaE  =   Parameters[2]*( s_i*cos(Phase_i)-s*cos(ph)) +
                        Parameters[1]*( s_i*s_neig[1]*cos(Phase_i-Phase_neig[1])-sp*s_neig[1]*cos(ph-Phase_neig[1])+
                                        s_i*s_neig[2]*cos(Phase_i-Phase_neig[2])-sp*s_neig[2]*cos(ph-Phase_neig[2])+
                                        s_i*s_neig[3]*cos(Phase_i-Phase_neig[3])-sp*s_neig[3]*cos(ph-Phase_neig[3])+
                                        s_i*s_neig[4]*cos(Phase_i-Phase_neig[4])-sp*s_neig[4]*cos(ph-Phase_neig[4])
                                        )

            # If the state is accepted the spin matrix, the energy, the magnetization,
            # and the concentration of spin 0 gets updated
            if rand() < min(1.0, exp(-DeltaE/T))
                E_sys += DeltaE
                Mag += sp*(cos(ph)+im*sin(ph))-S[i]*(cos(Phase[i])+im*sin(Phase[i]))
                Phase[i] = ph
                S[i] = sp
                if flags[4]
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
                    Corr =  S.*circshift(S,(0,r)).*cos(Phase - circshift(Phase,(0,r))) +
                            S.*circshift(S,(0,-r)).*cos(Phase - circshift(Phase,(0,-r))) +
                            S.*circshift(S,(r,0)).*cos(Phase - circshift(Phase,(r,0))) +
                            S.*circshift(S,(-r,0)).*cos(Phase - circshift(Phase,(-r,0)))
                    Corr_array[r+1] += 0.25*mean(Corr)
                end
            end

        end

        # Save the final configuration
        if flags[2]
            data_save(Files[4]*string("/T",T,".csv"), "a", string(reshape(S,1, N))[2:end-1]*"\n")
        end

        # Take averages
        DATA = DATA/Save_steps
        Corr_array = Corr_array/MC_array[3]

        if flags[4]
            VORTEX = VORTEX/Save_steps
            HELICITY_MOD = HELICITY_MOD/Save_steps
            HELICITY_MOD = -0.5*DATA[1] - HELICITY_MOD/T
        end

        # Calculate the suceptibilities
        Cv = (DATA[2] - DATA[1]*DATA[1])/(T*T)
        Chim = (DATA[4] - norm(DATA[3])^2)/T

        # The function returns
        # Energy per site
        # Magnetization x per site
        # Magnetization y per site
        # Cv per site
        # Chim per site
        # Vortex x per site
        # Vortex y per site
        # Helicity x per site
        # Helicity y per site
        # Correlation length
        return real(DATA[1]/N), real(norm(DATA[3])/N), real(Cv/N), real(Chim/N), real(VORTEX[1]/N), real(VORTEX[2]/N), real(HELICITY_MOD[1]/N), real(HELICITY_MOD[2]]/N), Corr_array
    end
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
end
