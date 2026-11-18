using JuMP, Ipopt
using NLsolve

include("lagrange.jl")

#Symbolic-Formeln kompilieren, sodass sie vom Solver verwendet werden können
M_func, _ = Symbolics.build_function(M_q, [q; qP; params...])
f_func, _ = Symbolics.build_function(f_vektor, [q; qP; M; params...])
g_flug_func, _ = Symbolics.build_function(g_flug, [q; params...])
g_kontakt_func, _ = Symbolics.build_function(g_kontakt, [q; params...])
F_N_func, _ = Symbolics.build_function(F_N, [q; qP; qPP; params...])


#Ipopt-Modell bauen
model = Model(Ipopt.Optimizer)
set_attribute(model, "max_iter", 1000)

N = 1000  # Anzahl Zeitschritte
dt = 0.01 # Zeitschrittgröße

@variable(model, q[1:4, 1:N]) #
@variable(model, qP[1:4, 1:N])
@variable(model, M[1:N])


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
@constraint(model, q[:, 1] == [start_phi, start_z, start_b30, start_b20])
@constraint(model, qP[:, 1] == zeros(4))

#Endwert: Sprunghöhe (muss noch wirklich ausgewählt werden)
zie_hoehe = 2.0
@constraint(model, q[2, N] = ziel_hoehe)

abflug_zeit = 3.0 #sekunden
k_abflug = (int)abflug_zeit / dt #Iterationsschritt, bei dem der Abflug geschieht

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