using JuMP, Ipopt
using NLsolve

include("lagrange.jl")

#Symbolic-Formeln in funktionen kompilieren, sodass sie vom Solver verwendet werden können
#TODO: Symbolische Variablen die eh konstant bleiben (z.B. l, m_1) vorher zu konstanten mit substitution über dein Dict machen
M_func, _ = Symbolics.build_function(M_q, [q; qP; params...])
f_func, _ = Symbolics.build_function(f_vektor, [q; qP; M; params...])
g_flug_x_func, _ = Symbolics.build_function(g_flug[1], [q; params...])
g_flug_z_func, _ = Symbolics.build_function(g_flug[2], [q; params...])
g_kontakt_boden_func, _ = Symbolics.build_function(g_kontakt[3], [q; params...])
F_N_func, _ = Symbolics.build_function(F_N, [q; qP; qPP]) #hier auf jeden Fall Symbolics mit Konstanten ersetzen! 


#Ipopt-Modell bauen
model = Model(Ipopt.Optimizer)
set_attribute(model, "max_iter", 1000)

N = 1000  # Anzahl Zeitschritte
dt = 0.01 # Zeitschrittgröße

@variable(model, q[1:4, 1:N]) #
@variable(model, qP[1:4, 1:N])
@variable(model, M[1:N])

#Festen Abflugzeitpunkt festlegen
abflug_zeit = 3.0 #sekunden
k_abflug = Int(abflug_zeit / dt) #Iterationsschritt, bei dem der Abflug geschieht

#Funktionen in JuMP registrieren
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

#Normalkraft, unterschieden zwischen Phasen
@variable(model, F_norm[1:N])

#F_norm[1] = m_ges * erdbeschl_g
for k in 2:(N-1)
    if k <= k_abflug #Kontaktphase
        @NLconstraint(model, F_norm[k] >= 0)
        @NLconstraint(model, F_norm[k] == F_N_reg(
            q[1,k], q[2,k], q[3,k], q[4,k],
            qP[1,k], qP[2,k], qP[3,k], qP[4,k],
            # qPP eingesetzt via Finite Differenzen:
            (q[1,k+1] - 2*q[1,k] + q[1,k-1]) / (dt^2),
            (q[2,k+1] - 2*q[2,k] + q[2,k-1]) / (dt^2),
            (q[3,k+1] - 2*q[3,k] + q[3,k-1]) / (dt^2),
            (q[4,k+1] - 2*q[4,k] + q[4,k-1]) / (dt^2)
        ))
    else #Flugphase
        @constraint(model, F_norm[k] == 0)
    end
end

#Drehmoment einschränken
M_max = 150.0
@variable(model, M[1:N])
for k in 1:N
    set_lower_bound(M[k], -M_max)
    set_upper_bound(M[k], M_max)
end

#Startwerte u.ä. übergeben
#erzwingt aus startwert phi (in p) andere startwerte z, b30, b20
function mechanismus!(F, x, p)
    z, b30, b20 = x[1], x[2], x[3]
    phi, ax, az, L, l, b, gam_B, l_F, gam_F = p 
    
    F[1] = ax + L * sin(b30) - l * sin(phi) - b * sin(b20 + gam_B) #[cite: 1]
    F[2] = az + L * cos(b30) - l * cos(phi) - b * cos(b20 + gam_B) #[cite: 1]
    F[3] = z - l * cos(phi) - l_F * cos(b20 + gam_F)
end

phi_start = deg2rad(30.0)
p_werte = [phi_start, a_x, a_z, L, l, b, gam_B, l_F, gam_F]
start_guess = [1.0, 0.0, 0.0]
ergebnis = nlsolve((F, x) -> mechanismus!(F, x, p_werte), start_guess)
z_start   = ergebnis.zero[1]
b30_start = ergebnis.zero[2]
b20_start = ergebnis.zero[3]

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
@constraint(model, q[2, N] = ziel_hoehe)

#Nebenbedingungen (Bewegungsgl., geom. Bedingungen)

#TODO:
#Nebenbedingungen einbringen: Bewegungsgleichungen und Zwangsbedingungen g
#WICHTIG: dem Solver muss k_abflug von vorne herein übergeben werden,
#sowie die Bedingung, dass davor F>0 und bei k_abflug F==0 sein muss.
#Der Solver kann nicht spontan entscheiden, wann Abflug ist.

#Zielfunktion J: Energieminimierung
@objective(model, Min, sum(0.5 * (M[k]^2 + M[k+1]^2) * dt for k in 1:(N-1)))

optimize!(model)
optimal_M = value.(M)
optimal_trajectory = value.(q)

println("\nProgrammende\n")