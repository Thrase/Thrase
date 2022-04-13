using DelimitedFiles


function make_ss(fault_coord, n, δNp, input_file)

    num = 1
    dir_name = string("../output_files/BP1_data", num)
    while isdir(dir_name)
        num += 1
        dir_name = string("../output_files/BP1_data", num)
    end
    
    mkdir(dir_name)
    cp(input_file, string(dir_name, "/input_file.dat"))

    slip_file = string(dir_name,  "/slip.dat")
    to_write = copy(fault_coord)
    pushfirst!(to_write, 0.0, 0.0)
    
    io = open(slip_file, "w")
    writedlm(io, [δNp])
    writedlm(io, to_write', " ")
    close(io)
    
    slip_rate_file = string(dir_name, "/slip_rate.dat")
    io = open(slip_rate_file, "w")
    writedlm(io, to_write', " ")
    close(io)
    
    stress_file = string(dir_name,  "/stress.dat")
    io = open(stress_file, "w")
    writedlm(io, to_write', " ")
    close(io)
    
    state_file = string(dir_name, "/state.dat")
    io = open(state_file, "w")
    writedlm(io, to_write', " ")
    close(io)
    
    return dir_name, slip_file, slip_rate_file, stress_file, state_file

end

function make_stations(dir_name)

    stat_depth = append!(collect(0.0:2.5:20.0), collect(25.0:5.0:35.0))
    file_names = [string(dir_name, "/station_", stat_depth[i]) for i in 1:length(stat_depth)]
        
    header = "t slip slip_rate shear_stress state\n"

    for file in file_names
        io = open(file, "w")
        write(io, header)
        close(io)
    end
    
    return file_names

end


function write_out(δ, V, τ, θ, fault_coord, Lw, t, file_names, η=nothing)


    stat_depth = append!(collect(0.0:2.5:20.0), collect(25.0:5.0:35.0))
    for i in 1:length(stat_depth)

        file_name = file_names[i]
        depth = stat_depth[i]
        d_ind = 0
        d_val = Lw
        for j in 1:length(fault_coord)
            if abs(depth-fault_coord[j]) < d_val
                d_ind = j
                d_val = abs(depth-fault_coord[j])
            end
        end

        @assert d_ind != 0

        x1 = fault_coord[d_ind]
        δ1 = δ[d_ind]
        V1 = V[d_ind]
        τ1 = τ[d_ind]
        θ1 = θ[d_ind]
        
        if fault_coord[d_ind] <= depth
            x2 = fault_coord[d_ind + 1]
            δ2 = δ[d_ind + 1]
            V2 = V[d_ind + 1]
            τ2 = τ[d_ind + 1]
            θ2 = θ[d_ind + 1]
        end
        
        if fault_coord[d_ind] > depth
            x2 = fault_coord[d_ind - 1]
            δ2 = δ[d_ind - 1]
            V2 = V[d_ind - 1]
            τ2 = τ[d_ind - 1]
            θ2 = θ[d_ind - 1]
        end

        δw = l_interp(depth, x1, x2, δ1, δ2)
        Vw = l_interp(depth, x1, x2, V1, V2)
        if η != nothing
            τw = l_interp(depth, x1, x2, τ1, τ2) - (Vw * η)
        else
            τw = l_interp(depth, x1, x2, τ1, τ2)
        end
        
        θw = l_interp(depth, x1, x2, θ1, θ2)
        

        dat = [t, δw, log(10,abs(Vw)), τw, log(10, θw)]
        io = open(file_name, "a")
        writedlm(io, dat', " ")
        close(io)
    end 

end



function write_out_ss(δ, V, τ, ψ, t, slip_file, stress_file, slip_rate_file, state_file, η=nothing)
    
    max_V = log(10, maximum(V))
    to_write = copy(δ)
    pushfirst!(to_write, t, max_V)
    io = open(slip_file, "a")
    writedlm(io, to_write')
    close(io)
    
    to_write = copy(τ)
    if η != nothing
        to_write = to_write .- (η .* V)
    end
    pushfirst!(to_write, t, max_V)
    io = open(stress_file, "a")
    writedlm(io, to_write')
    close(io)

    io = open(slip_rate_file, "a")
    writedlm(io, V')
    close(io)

    io = open(state_file, "a")
    writedlm(io, ψ')
    close(io)

    
end


function l_interp(x, x1, x2, y1, y2)

    return (y2 - y1)/(x2 - x1) * x - (y2 - y1)/(x2 - x1) * x1 + y1
    
end



