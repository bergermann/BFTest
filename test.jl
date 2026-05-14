
using BoostFractor, Plots



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

B = zeros(ComplexF64,M*(2L+1),length(freqs))

@time for i in eachindex(freqs)
    B[:,i] .= transformer(sbdry,coords,modes; prop=propagator,f=freqs[i],diskR=diskR)
end

plot(freqs/1e9,abs2.(B)'; label=["L=-1" "L= 0" "L= 1"])