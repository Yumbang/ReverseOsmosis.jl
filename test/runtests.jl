using ReverseOsmosis
using Unitful
using Test

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
const R_GAS = 8.3145e3                        # [J/(kmol·K)]
const H_CH  = 8.64e-4                         # channel height [m]
const L_MOD = 1.016                           # module length [m]
const W_CH  = 37.0 / L_MOD                    # channel width = area / length [m]
const RM_1  = 1 / (5.295039332042272 / 3600 / 1e3 / 1e5)   # A = 5.295 LMH/bar → R_m [Pa·s/m]
const KFP_1 = 0.674150595374451e9             # [Pa·s/m²]

vessel1(; n=100, Rm=RM_1, kfp=KFP_1, K=16.0, r=0.995) =
    pristine_vessel(7, n, L_MOD, H_CH, W_CH, Rm, kfp, K, r)
vessel2(; n=50, Rm=6e14, B=3.5e-8, kfp=0.8, K=16.0) =
    pristine_vessel2(7, n, L_MOD, 7e-4, W_CH, Rm, B, kfp, K)

"Run one single-pass step; `P_gauge` is the applied pressure above atmosphere [Pa]."
function run1(feed, P_gauge; vessel=vessel1(), dt=3600.0, fouling=false)
    proc = SinglePassRO(0.8, vessel)
    b, p, _, _, pw = process_singlepass_RO!(feed, (P_ATM + P_gauge) * u"Pa";
                                            process=proc, dt=dt * u"s", fouling=fouling)
    return b, p, pw, proc
end

salt_error(feed, b, p) = (p.Q * p.C + b.Q * b.C - feed.Q * feed.C) / (feed.Q * feed.C)
water_error(feed, b, p) = (p.Q + b.Q - feed.Q) / feed.Q

all_R_m(v) = [e.R_m for m in v.modules_array for e in m.elements_array]
all_R_c(v) = [e.R_c for m in v.modules_array for e in m.elements_array]

# -------------------------------------------------------------------
@testset "ReverseOsmosis.jl physics" begin

@testset "Osmotic pressure (van't Hoff)" begin
    # 35 kg/m³ NaCl at 25 °C: 2·C/M·R·T ≈ 29.7 bar
    π_sw = osmo_p(35.0, 25.0)
    @test π_sw ≈ 2 * 35.0 / 58.44 * R_GAS * 298.15
    @test 29.0e5 < π_sw < 30.5e5
    @test osmo_p(70.0, 25.0) ≈ 2 * osmo_p(35.0, 25.0)                # linear in C
    @test osmo_p(35.0, 50.0) / osmo_p(35.0, 25.0) ≈ 323.15 / 298.15  # linear in absolute T
    @test osmo_p(0.0, 25.0) == 0.0
end

@testset "Stream operations" begin
    a = Water(10.0, 20.0, 1.0, 5e5)
    b = Water(30.0, 30.0, 3.0, 5e5)
    m = mix(a, b)
    @test m.Q ≈ 40.0
    @test m.Q * m.C ≈ a.Q * a.C + b.Q * b.C        # salt mass conserved
    @test m.T ≈ (10 * 20 + 30 * 30) / 40           # flow-weighted temperature

    w = Water(10.0, 20.0, 1.0, P_ATM, 74.55)       # KCl
    @test pressurize(w, 5e5).m_avg == 74.55         # pressurize keeps ion mass
    @test pressurize(w, 5e5).P ≈ P_ATM + 5e5

    # Pure-water streams (e.g. two permeates with C = 0) should mix without NaN
    @test !isnan(mix(Water(1.0, 20.0, 0.0, P_ATM), Water(1.0, 20.0, 0.0, P_ATM)).m_avg)

    pw, power = pump(Water(36.0, 20.0, 1.0, P_ATM), 10e5; efficiency=0.8)
    @test pw.P ≈ P_ATM + 10e5
    @test ustrip(u"W", power) ≈ 36.0 / 3600 * 10e5 / 0.8
end

@testset "Model 1: single element" begin
    dx   = L_MOD / 100
    elem = MembraneElement(H_CH, W_CH, dx, 16.0, RM_1, KFP_1, 0.995)
    feed = Water(12.0, 20.0, 0.5, P_ATM + 9e5)
    b, p, ΔR = element_filtration(feed, elem; dt=3600.0)

    @test water_error(feed, b, p) ≈ 0 atol = 1e-12
    @test salt_error(feed, b, p) ≈ 0 atol = 1e-12
    @test p.C ≈ (1 - 0.995) * feed.C
    @test p.P == P_ATM
    @test b.P < feed.P                                  # friction loss

    # Flux follows v_w = (P - P_ATM - π) / R_m at 20 °C (TCF_A = 1).
    # For a thin element the outlet concentration ≈ feed concentration.
    v_w = p.Q / 3600 / (W_CH * dx)
    @test v_w ≈ (feed.P - P_ATM - osmo_p(feed.C, 20.0)) / RM_1 rtol = 1e-3
    @test ΔR ≈ KFP_1 * v_w * 3600.0 rtol = 1e-12

    # No driving force → no permeate, feed passes through
    b0, p0, _ = element_filtration(Water(12.0, 20.0, 35.0, P_ATM + 5e5), elem)
    @test p0.Q == 0.0
    @test b0.Q ≈ 12.0
    @test b0.C ≈ 35.0

    # Perfect rejection → salt-free permeate
    elem_r1 = MembraneElement(H_CH, W_CH, dx, 16.0, RM_1, KFP_1, 1.0)
    _, p1, _ = element_filtration(feed, elem_r1)
    @test p1.C == 0.0

    # No spacer friction → no pressure drop
    elem_k0 = MembraneElement(H_CH, W_CH, dx, 0.0, RM_1, KFP_1, 0.995)
    bk, _, _ = element_filtration(feed, elem_k0)
    @test bk.P ≈ feed.P
end

@testset "Model 1: pressure vessel" begin
    feed = Water(1000 / 84, 20.0, 0.5, P_ATM)
    b, p, pw, _ = run1(feed, 9e5)

    @test water_error(feed, b, p) ≈ 0 atol = 1e-10
    @test salt_error(feed, b, p) ≈ 0 atol = 1e-10
    @test 0 < p.Q < feed.Q
    @test p.C < feed.C < b.C
    @test ustrip(u"W", pw) ≈ feed.Q / 3600 * 9e5 / 0.8

    # Regression against an independent implementation (ROMARL simulator with
    # cumulative salt balance and flow-weighted permeate concentration).
    @test p.Q ≈ 10.0816 rtol = 1e-4
    @test b.C ≈ 3.23436 rtol = 1e-4
    @test (b.P - P_ATM) / 1e5 ≈ 7.943 rtol = 1e-3

    # Monotonic responses
    Qp(; P=9e5, C=0.5, T=20.0) = run1(Water(1000 / 84, T, C, P_ATM), P)[2].Q
    @test Qp(P=7e5) < Qp(P=9e5) < Qp(P=11e5)            # more pressure → more permeate
    @test Qp(C=1.0) < Qp(C=0.5)                          # more salt → less permeate
    @test Qp(T=15.0) < Qp(T=20.0) < Qp(T=25.0)           # warmer → more permeable

    # Spatial discretisation: result should converge with segment count
    Q50  = run1(feed, 9e5; vessel=vessel1(n=50))[2].Q
    Q100 = run1(feed, 9e5; vessel=vessel1(n=100))[2].Q
    @test Q50 ≈ Q100 rtol = 1e-3
end

@testset "Model 1: fouling" begin
    feed = Water(1000 / 84, 20.0, 0.5, P_ATM)
    v    = vessel1(n=20)
    proc = SinglePassRO(0.8, v)
    Rm0  = all_R_m(v)
    Qs   = Float64[]
    for _ in 1:5
        _, p, _, _, _ = process_singlepass_RO!(feed, (P_ATM + 8e5) * u"Pa";
                                               process=proc, dt=21600.0u"s", fouling=true)
        push!(Qs, p.Q)
    end
    @test all(all_R_m(v) .> Rm0)                # resistance only increases
    @test issorted(Qs; rev=true)                # flux declines over time

    v_off = vessel1(n=20)
    run1(feed, 8e5; vessel=v_off, fouling=false)
    @test all_R_m(v_off) == fill(RM_1, length(all_R_m(v_off)))   # no fouling when disabled
end

@testset "Model 2: solution-diffusion + cake layer" begin
    feed = Water2(10.0, 25.0, 35.0, 0.1, P_ATM)
    b, p, _, proc = run1(feed, 54e5; vessel=vessel2())

    @test water_error(feed, b, p) ≈ 0 atol = 1e-10
    @test salt_error(feed, b, p) ≈ 0 atol = 1e-10
    @test p.Cf == 0.0                                         # foulant fully rejected
    @test b.Cf * b.Q ≤ feed.Cf * feed.Q                       # foulant cannot be created

    # Higher salt permeability → saltier permeate
    _, pB, _, _ = run1(feed, 54e5; vessel=vessel2(B=7e-8))
    @test pB.C > p.C

    # Below critical flux there is no cake; above it a cake forms
    _, _, _, low = run1(feed, 54e5; vessel=vessel2(Rm=8e14), fouling=true)
    @test all(all_R_c(low.pressure_vessel) .== 0.0)
    @test all(all_R_c(proc.pressure_vessel) .== 0.0)          # fouling=false leaves the cake untouched
    _, _, _, high = run1(feed, 54e5; vessel=vessel2(Rm=4e14), fouling=true)
    @test maximum(all_R_c(high.pressure_vessel)) > 0.0

    # Cake resistance must not depend on the segment count
    Rc1(n) = all_R_c(run1(feed, 54e5; vessel=vessel2(n=n, Rm=4e14), fouling=true)[4].pressure_vessel)[1]
    @test Rc1(25) ≈ Rc1(100) rtol = 0.02

    # High-flux case that used to trip the iteration limit right as it converged
    b_hi, p_hi, _, _ = run1(feed, 54e5; vessel=vessel2(Rm=2e14))
    @test salt_error(feed, b_hi, p_hi) ≈ 0 atol = 1e-10

    # Driving pressure below osmotic pressure: permeate concentration must stay finite
    b_lo, p_lo, _, _ = run1(Water2(10.0, 25.0, 35.0, 0.1, P_ATM), 20e5; vessel=vessel2())
    @test p_lo.Q ≥ 0.0
    @test p_lo.C ≤ feed.C
end

end
