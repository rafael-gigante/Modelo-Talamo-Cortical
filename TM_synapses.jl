module TMsynapses

export TMsynE_inst_D, TMsynE_inst_F, TMsynE_inst, TMsynI_inst

function TMsynE_inst_D(r, x, Is, sp_event)
    dt = 0.1
    
    tauf = [670, 17, 326]
    taud = [138, 671, 329]
    U = [0.09, 0.5, 0.29]
    A = [0, 1, 0] # Somente a sinapse depressora (D) contribui
    
    taus = 3
    
    r = r .+ dt .* ((-r ./ tauf ) .+ (U .* (1 .- r) .* sp_event))
    x = x .+ dt .* (((1 .- x) ./ taud) .- ((r .+ U .* (1 .- r)) .* x .* sp_event))
    Is = Is .+ dt .* ((-Is ./ taus) .+ A .* ((r .+ U .* (1 .- r)) .* x .* sp_event))
    
    Ipost = sum(Is)
    return r, x, Is, Ipost
end

function TMsynE_inst_F(r, x, Is, sp_event)
    dt = 0.1
    
    tauf = [670, 17, 326]
    taud = [138, 671, 329]
    U = [0.09, 0.5, 0.29]
    A = [1, 0, 0]  # Somente a sinapse facilitadora (F) contribui
    taus = 3
    
    r = r .+ dt .* ((-r ./ tauf ) .+ (U .* (1 .- r) .* sp_event))
    x = x .+ dt .* (((1 .- x) ./ taud) .- ((r .+ U .* (1 .- r)) .* x .* sp_event))
    Is = Is .+ dt .* ((-Is ./ taus) .+ A .* ((r .+ U .* (1 .- r)) .* x .* sp_event))
    
    Ipost = sum(Is)
    return r, x, Is, Ipost
end

function TMsynE_inst(r, x, Is, sp_event; dt=0.1)
    tauf = [670, 17, 326]
    taud = [138, 671, 329]
    U = [0.09, 0.5, 0.29]
    A = [0.20, 0.63, 0.17]
    
    taus = 3
    
    r = r .+ dt .* ((-r ./ tauf ) .+ (U .* (1 .- r) .* sp_event))
    x = x .+ dt .* (((1 .- x) ./ taud) .- ((r .+ U .* (1 .- r)) .* x .* sp_event))
    Is = Is .+ dt .* ((-Is ./ taus) .+ A .* ((r .+ U .* (1 .- r)) .* x .* sp_event))
    
    Ipost = sum(Is)
    return r, x, Is, Ipost
end

function TMsynI_inst(r, x, Is, sp_event; dt=0.1)
    tauf = [376, 21, 62]
    taud = [45, 706, 144]
    U = [0.016, 0.25, 0.32]
    A = [0.08, 0.75, 0.17]
    
    taus = 11
    
    r = r .+ dt .* ((-r ./ tauf ) .+ (U .* (1 .- r) .* sp_event))
    x = x .+ dt .* (((1 .- x) ./ taud) .- ((r .+ U .* (1 .- r)) .* x .* sp_event))
    Is = Is .+ dt .* ((-Is ./ taus) .+ A .* ((r .+ U .* (1 .- r)) .* x .* sp_event))

    Ipost = sum(Is)
    return r, x, Is, Ipost
end

end