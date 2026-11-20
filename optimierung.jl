using JuMP, Ipopt
using NLsolve

include("lagrange.jl")

# Konstanten bestimmen
const_map = Dict(
    L => 0.6, l => 0.5, l_C => 0.1,
    d => 0.2, d_x => 0.15, d_z => 0.13,
    a_z => 0.8, a_x => 0.1,
    b => 0.3, b_x => 0.2, b_z => 0.22,
    l_F => 0.52, gam_F => 0.2, gam_D => 0.3, gam_B => 0.25,
    l_S4x => 0.0, l_S4z => 0.1, l_S2 => 0.25, gam_S2 => 0.1,
    m_1 => 2.0, m_2 => 3.0, m_3 => 2.5, m_4 => 1.5,
    J_1 => 0.05, J_2 => 0.08, J_3 => 0.06,
    l_0 => 0.25, k => 300.0, g => 9.81
)

# Aus const_map die reinen Float64-Werte ziehen
ax_val    = const_map[a_x]
az_val    = const_map[a_z]
L_val     = const_map[L]
l_val     = const_map[l]
b_val     = const_map[b]
gam_B_val = const_map[gam_B]
l_F_val   = const_map[l_F]
gam_F_val = const_map[gam_F]
m_ges_val = const_map[m_1] + const_map[m_2] + const_map[m_3] + const_map[m_4]
g_val     = const_map[g]

#konstanten einsetzen in Gleichungen
M_q_num             = substitute.(M_q, (const_map,))
f_vektor_num        = substitute.(f_vektor, (const_map,))
Q_ZB_num            = substitute.(Q_ZB, (const_map,))
Q_ZF_num            = substitute.(Q_ZF_mit_Fnorm, (const_map,))
g_flug_x_num        = substitute(g_flug[1], const_map)
g_flug_z_num        = substitute(g_flug[2], const_map)
g_kontakt_boden_num = substitute(g_kontakt[3], const_map)
F_N_num             = substitute(F_N, const_map)

# Funktionen kompilieren
M_func, _             = Symbolics.build_function(M_q_num, [q; qP], expression = Val{false})
f_func, _             = Symbolics.build_function(f_vektor_num, [q; qP; M], expression = Val{false})
g_flug_x_func         = Symbolics.build_function(g_flug_x_num, q, expression = Val{false})
g_flug_z_func         = Symbolics.build_function(g_flug_z_num, q, expression = Val{false})
g_kontakt_boden_func  = Symbolics.build_function(g_kontakt_boden_num, q, expression = Val{false})
F_N_func              = Symbolics.build_function(F_N_num, [q; qP; qPP], expression = Val{false})


#Bewegungsgleichung umgestellt auf linke Seite - rechte Seite = Residuum =!= 0 (Nebenbedingung für optimizer)
res_flug_expr    = [sum(M_q_num[i, j] * qPP[j] for j in 1:4) - (f_vektor_num[i] + Q_ZB_num[i])                for i in 1:4]
res_kontakt_expr = [sum(M_q_num[i, j] * qPP[j] for j in 1:4) - (f_vektor_num[i] + Q_ZB_num[i] + Q_ZF_num[i]) for i in 1:4]

#Bewegungsgleichungen("Equations of Motion - EOM") kompilieren
eom_flug_1    = Symbolics.build_function(res_flug_expr[1],    [q; qP; qPP; M; F_Bx; F_Bz], expression = Val{false})
eom_flug_2    = Symbolics.build_function(res_flug_expr[2],    [q; qP; qPP; M; F_Bx; F_Bz], expression = Val{false})
eom_flug_3    = Symbolics.build_function(res_flug_expr[3],    [q; qP; qPP; M; F_Bx; F_Bz], expression = Val{false})
eom_flug_4    = Symbolics.build_function(res_flug_expr[4],    [q; qP; qPP; M; F_Bx; F_Bz], expression = Val{false})

eom_kontakt_1 = Symbolics.build_function(res_kontakt_expr[1], [q; qP; qPP; M; F_Bx; F_Bz; F_norm_var], expression = Val{false})
eom_kontakt_2 = Symbolics.build_function(res_kontakt_expr[2], [q; qP; qPP; M; F_Bx; F_Bz; F_norm_var], expression = Val{false})
eom_kontakt_3 = Symbolics.build_function(res_kontakt_expr[3], [q; qP; qPP; M; F_Bx; F_Bz; F_norm_var], expression = Val{false})
eom_kontakt_4 = Symbolics.build_function(res_kontakt_expr[4], [q; qP; qPP; M; F_Bx; F_Bz; F_norm_var], expression = Val{false})

#Ipopt-Modell bauen
model = Model(Ipopt.Optimizer)
set_attribute(model, "max_iter", 500)

N = 100  # Anzahl Zeitschritte
dt = 0.1 # Zeitschrittgröße

#Variablen einführen für Optimizer
JuMP.@variable(model, q[1:4, 1:N])
JuMP.@variable(model, qP[1:4, 1:N])
JuMP.@variable(model, M[1:N])
JuMP.@variable(model, F_Bx[1:N])
JuMP.@variable(model, F_Bz[1:N])
JuMP.@variable(model, F_norm[1:N])


#Bewegungsgleichungen dem Optimizer registrieren
register(model, :eom_flug_reg_1,    15, eom_flug_1;    autodiff = true)
register(model, :eom_flug_reg_2,    15, eom_flug_2;    autodiff = true)
register(model, :eom_flug_reg_3,    15, eom_flug_3;    autodiff = true)
register(model, :eom_flug_reg_4,    15, eom_flug_4;    autodiff = true)

register(model, :eom_kontakt_reg_1, 16, eom_kontakt_1; autodiff = true)
register(model, :eom_kontakt_reg_2, 16, eom_kontakt_2; autodiff = true)
register(model, :eom_kontakt_reg_3, 16, eom_kontakt_3; autodiff = true)
register(model, :eom_kontakt_reg_4, 16, eom_kontakt_4; autodiff = true)

#q an qPunkt koppeln für alle Zeitschritte
for k in 1:(N-1), i in 1:4
    @constraint(model, q[i, k+1] - q[i, k] == 0.5 * dt * (qP[i, k] + qP[i, k+1]))
end

#Festen Abflugzeitpunkt festlegen
abflug_zeit = 3.0 #sekunden
k_abflug = Int(abflug_zeit / dt) #Iterationsschritt, bei dem der Abflug geschieht


#Bewegunsgleichungen als Nebenbedingungen einführen
for k in 2:(N-1)
    qPP1 = (q[1,k+1] - 2*q[1,k] + q[1,k-1]) / (dt^2)
    qPP2 = (q[2,k+1] - 2*q[2,k] + q[2,k-1]) / (dt^2)
    qPP3 = (q[3,k+1] - 2*q[3,k] + q[3,k-1]) / (dt^2)
    qPP4 = (q[4,k+1] - 2*q[4,k] + q[4,k-1]) / (dt^2)

    if k <= k_abflug
        @NLconstraint(model, eom_kontakt_reg_1(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k], F_norm[k]) == 0)
        @NLconstraint(model, eom_kontakt_reg_2(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k], F_norm[k]) == 0)
        @NLconstraint(model, eom_kontakt_reg_3(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k], F_norm[k]) == 0)
        @NLconstraint(model, eom_kontakt_reg_4(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k], F_norm[k]) == 0)
    else
        @NLconstraint(model, eom_flug_reg_1(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k]) == 0)
        @NLconstraint(model, eom_flug_reg_2(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k]) == 0)
        @NLconstraint(model, eom_flug_reg_3(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k]) == 0)
        @NLconstraint(model, eom_flug_reg_4(q[1,k], q[2,k], q[3,k], q[4,k], qP[1,k], qP[2,k], qP[3,k], qP[4,k], qPP1, qPP2, qPP3, qPP4, M[k], F_Bx[k], F_Bz[k]) == 0)
    end
end

#andere Funktionen in JuMP registrieren
register(model, :g_flug_x_reg, 4, g_flug_x_func; autodiff = true)
register(model, :g_flug_z_reg, 4, g_flug_z_func; autodiff = true)
register(model, :g_kontakt_reg, 4, g_kontakt_boden_func; autodiff = true)
register(model, :F_N_reg, 12, F_N_func; autodiff = true)


#geometrische Nebenbedingungen
for k in 1:N
    # Gelenk B immer geschlossen
    @NLconstraint(model, g_flug_x_reg(q[1,k], q[2,k], q[3,k], q[4,k]) == 0)
    @NLconstraint(model, g_flug_z_reg(q[1,k], q[2,k], q[3,k], q[4,k]) == 0)
end
for k in 1:k_abflug
    #Fuß am Boden bis zum Abflug-Zeitpunkt
    @NLconstraint(model, g_kontakt_reg(q[1,k], q[2,k], q[3,k], q[4,k]) == 0)
end

#F_norm[1] = m_ges * erdbeschl_g
for k in 2:(N-1)
    if k <= k_abflug #Kontaktphase
        @NLconstraint(model, F_norm[k] >= 0)
        #=
        #muss nicht definiert werden, F_norm wird schon über Zwangsbedingungen hergeleitet vom Optimizer
        @NLconstraint(model, F_norm[k] == F_N_reg(
            q[1,k], q[2,k], q[3,k], q[4,k],
            qP[1,k], qP[2,k], qP[3,k], qP[4,k],
            # qPP[1...4] näherungsweise eingesetzt:
            (q[1,k+1] - 2*q[1,k] + q[1,k-1]) / (dt^2),
            (q[2,k+1] - 2*q[2,k] + q[2,k-1]) / (dt^2),
            (q[3,k+1] - 2*q[3,k] + q[3,k-1]) / (dt^2),
            (q[4,k+1] - 2*q[4,k] + q[4,k-1]) / (dt^2)
        ))
        =#
    else #Flugphase
        @constraint(model, F_norm[k] == 0)
    end
end
#zum Abflugzeitpunkt ist die Kraft dann genau 0 → Roboter steuert auf diesen Absprungzeitpunkt hin
@NLconstraint(model, F_norm[k_abflug] == 0)


#Drehmoment einschränken
M_max = 150.0
for k in 1:N
    set_lower_bound(M[k], -M_max)
    set_upper_bound(M[k], M_max)
end

#Startwerte u.ä. übergeben
#erzwingt aus startwert phi (in p) andere startwerte z, b30, b20
function start_kinematik!(F, x, phi_val)
    q_try = [phi_val, x[1], x[2], x[3]]
    F[1] = g_flug_x_func(q_try)
    F[2] = g_flug_z_func(q_try)
    F[3] = g_kontakt_boden_func(q_try)
end

phi_start = deg2rad(30.0)
sol = nlsolve((F, x) -> start_kinematik!(F, x, phi_start), [1.0, 0.0, 0.0])
@show sol.f_converged sol.x_converged sol.residual_norm
z_start, b30_start, b20_start = sol.zero

#übergeben von schätzungswerten für den solver: Startwerte für jeden Zeitschritt + konst. Drehmoment
for k in 1:N
    set_start_value(q[1, k], phi_start) 
    set_start_value(q[2, k], z_start)
    set_start_value(q[3, k], b30_start)
    set_start_value(q[4, k], b20_start)
    set_start_value(M[k], 10.0)
end

#Startwerte für Koord. und Geschw. übergeben
@constraint(model, q[:, 1] == [phi_start, z_start, b30_start, b20_start])
@constraint(model, qP[:, 1] == zeros(4))

#Endwert: Sprunghöhe (muss noch wirklich ausgewählt werden)
ziel_hoehe = 2.0
JuMP.@variable(model, z_max)

# z_max steigt dann an, während z steigt und bleibt beim höchstwert wenn der roboter fällt
for k in 1:N
    @constraint(model, z_max >= q[2, k])
end
@constraint(model, z_max >= ziel_hoehe)

#Zielfunktion J: Energieminimierung
@objective(model, Min, sum(0.5 * (M[k]^2 + M[k+1]^2) * dt for k in 1:(N-1)))

optimize!(model)
println("Termination status: ", termination_status(model))
report = primal_feasibility_report(model)
for (con, viol) in sort(collect(report), by = x -> -x[2])[1:min(8, length(report))]
    println("Verletzung = ", viol, "  |  ", con)
end

println("q[:, 2]  = ", value.(q[:, 2]))
println("qP[:, 2] = ", value.(qP[:, 2]))
println("M[2]     = ", value(M[2]))
println("F_Bx[2]  = ", value(F_Bx[2]))
println("F_Bz[2]  = ", value(F_Bz[2]))
println("F_norm[2]= ", value(F_norm[2]))

optimal_M = value.(M)
optimal_trajectory = value.(q)

println("\nProgrammende\n")