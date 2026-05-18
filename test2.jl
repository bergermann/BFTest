
using BoostFractor, Plots
using BoostFractor: e_field_dimensions




function axion_induced_modes(coords::CoordinateSystem, modes::Modes; B=nothing, diskR=0.15)
    if B === nothing
        if e_field_dimensions(modes) == 1
            B = ones(length(coords.X), length(coords.Y), e_field_dimensions(modes))
        elseif e_field_dimensions(modes) == 3
            B = zeros(length(coords.X), length(coords.Y), e_field_dimensions(modes))
            B[:,:,2] = ones(length(coords.X), length(coords.Y))
        end
    end

    B .*= [sqrt(x^2 + y^2) <= diskR for x in coords.X, y in coords.Y]
    B ./= sqrt.( sum(abs2.(B)) )

    modes_initial = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1)))
    for m in 1:modes.M, l in -modes.L:modes.L
        modes_initial[(m-1)*(2modes.L+1)+l+modes.L+1] =
                sum( conj.(modes.mode_patterns[m,l+modes.L+1,:,:,:]) .* B )
    end

    return modes_initial
end


function propagation_matrix(dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords::CoordinateSystem, modes::Modes;
        is_air=(real(eps)==1), onlydiagonal=false, prop=propagator)

    matching_matrix = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1),modes.M*(2modes.L+1)))

    k0 = 2pi/lambda*sqrt(eps)

    propfunc = nothing # initialize
    if is_air
        function propagate(x)
            return prop(copy(x), dz, diskR, eps, tilt_x, tilt_y, surface, lambda, coords)
        end
        propfunc = propagate
    else
        propfunc(efields) = efields.*[exp(-1im*k0*tilt_x*x) * exp(-1im*k0*tilt_y*y) for x in coords.X, y in coords.Y].*exp.(-1im*k0*surface)
    end
    
    for m_prime in 1:modes.M, l_prime in -modes.L:modes.L
        for i in 1:e_field_dimensions(modes)
            propagated = propfunc(modes.mode_patterns[m_prime,l_prime+modes.L+1,:,:,i])

            for m in (onlydiagonal ? [m_prime] : 1:modes.M), l in (onlydiagonal ? [l_prime] : -modes.L:modes.L)
                matching_matrix[(m-1)*(2modes.L+1)+l+modes.L+1, (m_prime-1)*(2modes.L+1)+l_prime+modes.L+1] +=
                        sum( conj.(modes.mode_patterns[m,l+modes.L+1,:,:,i]) .* propagated )
            end
        end
    end

    if !is_air
        propagation_matrix = Array{Complex{Float64}}(zeros(modes.M*(2modes.L+1),modes.M*(2modes.L+1)))
        for m in 1:modes.M, l in -modes.L:modes.L
            kz = sqrt(k0^2 - modes.mode_kt[m,l+modes.L+1]^2)
            propagation_matrix[(m-1)*(2modes.L+1)+l+modes.L+1, (m-1)*(2modes.L+1)+l+modes.L+1] = exp(-1im*kz*dz)
        end
        
        matching_matrix = propagation_matrix*matching_matrix
    end

    return matching_matrix
end

function get_boundary_matrix(n_left, n_right, diffprop::Array{Complex{T}}, modes::Modes) where T<:Real
    G = (( (1. /(2*n_right)).*[(n_right+n_left)*modes.id (n_right-n_left)*modes.id ; (n_right-n_left)*modes.id (n_right+n_left)*modes.id] ))


    return G * [diffprop modes.zeromatrix; modes.zeromatrix inv(Array{Complex{T}}(diffprop))]
end

function axion_contrib(T,n1,n0, initial, modes::Modes)
    axion_beam = axion_S_factor(n1,n0) .* (T[index(modes,2),index(modes,1)]*(copy(initial)) + T[index(modes,2),index(modes,2)]*(copy(initial)))
    return axion_beam
end

function axion_S_factor(n1,n0)
    return (1. /n1^2 - 1. /n0^2)/2
end

function calc_propagation_matrices(bdry::SetupBoundaries, coords::CoordinateSystem, modes::Modes; f=10.0e9, prop=propagator, diskR=0.15)
    Nregions = length(bdry.eps)
    lambda = wavelength(f)
    return [ propagation_matrix(bdry.distance[i], diskR, bdry.eps[i],
            bdry.relative_tilt_x[i], bdry.relative_tilt_y[i], bdry.relative_surfaces[i,:,:], lambda, coords, modes;
            prop=prop) for i in 1:(Nregions) ]
end

function transformer(bdry::SetupBoundaries, coords::CoordinateSystem, modes::Modes;
        f=10.0e9, prop=propagator,
        propagation_matrices::Array{Array{Complex{T},2},1}=Array{Complex{Float64},2}[],
        diskR=0.15, emit=axion_induced_modes(coords,modes;B=nothing,diskR=diskR),
        reflect=nothing) where T<:Real
    
    bdry.eps[isnan.(bdry.eps)] .= 1e30

    transmissionfunction_complete = [modes.id modes.zeromatrix ; modes.zeromatrix modes.id ]
    lambda = wavelength(f)

    initial = emit

    axion_beam = Array{Complex{T}}(zeros((modes.M)*(2modes.L+1)))

    Nregions = length(bdry.eps)
    idx_reg(s) = Nregions-s+1


    for s in (Nregions-1):-1:1
        axion_beam .+= axion_contrib(transmissionfunction_complete, sqrt(bdry.eps[idx_reg(s+1)]), sqrt(bdry.eps[idx_reg(s)]), initial, modes)

        diffprop = (isempty(propagation_matrices) ?
                        propagation_matrix(bdry.distance[idx_reg(s)], diskR, bdry.eps[idx_reg(s)], bdry.relative_tilt_x[idx_reg(s)], bdry.relative_tilt_y[idx_reg(s)], bdry.relative_surfaces[idx_reg(s),:,:], lambda, coords, modes; prop=prop) :
                        propagation_matrices[idx_reg(s)])

        transmissionfunction_complete *= get_boundary_matrix(sqrt(bdry.eps[idx_reg(s)]), sqrt(bdry.eps[idx_reg(s+1)]), diffprop, modes)
    end

    boost =  - (transmissionfunction_complete[index(modes,2),index(modes,2)]) \ (axion_beam)

    if reflect === nothing
        return boost
    end

    refl = - transmissionfunction_complete[index(modes,2),index(modes,2)] \
           ((transmissionfunction_complete[index(modes,2),index(modes,1)]) * (reflect))
    return boost, refl
end






freqs = range(18e9,24e9,10_000);


coords = SeedCoordinateSystem(X=-0.5:0.02:0.5,Y=-0.5:0.02:0.5)
diskR = 0.15

eps = Complex{Float64}[NaN,1,24.0,1]
distance = [0,7,1,0]*1e-3
tiltx = deg2rad.([0,1,0,-1])
tilty = deg2rad.([0,0,0,0])

sbdry = SeedSetupBoundaries(coords; diskno=1,distance=distance,epsilon=eps,
    relative_tilt_x=tiltx,relative_tilt_y=tilty)

M = 1; L = 1
modes = SeedModes(coords, ThreeDim=true, Mmax=M, Lmax=L, diskR=diskR)


ax = axion_induced_modes(coords,modes)