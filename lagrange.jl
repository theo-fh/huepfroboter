using Symbolics
using Latexify

#generalisierte Koordinaten
@variables phi, z, b30, b20
#ableitungen 
@variables phiP, zP, b30P, b20P
@variables phiPP, zPP, b30PP, b20PP
#eingeprägte Kräfte
@variables M, F2_Bx, F2_Bz, F3_Bx, F3_Bz


#Konstanten
@variables L, l, l_C, d, d_x, d_z, a_z, a_x, b, b_x, b_z, l_F, gam_F, gam_D, gam_B, l_S4x, l_S4z, l_S2, gam_S2, m_1, m_2, m_3, m_4, J_1, J_2, J_3, l_0, k, g

#=
#Konstanten für Schwerpunkt-/Gelenkkoordinaten
d = sqrt(d_x^2 + d_z^2)
gam_D = atan(d_x/d_z)
b = sqrt(b_x^2 + b_z^2)
gam_B = atan(b_x/b_z)
l_F = sqrt(l^2 + d_x^2)
gam_F = atan(d_x/l)
=#

#potentielle Energie (Schwerkraft)
V_1 = m_1*g*(z - l/2 * cos(phi))
V_2 = m_2*g*(z - l*cos(phi) - l_S2 * cos(b20 + gam_S2))
V_3 = m_3*g*(z - a_z - L/2 * cos(b30))
V_4 = m_4*g*(z - l_S4z)

#potentielle Energie (Feder)
V_F = 1/2 * k * (sqrt((l-l_C)^2 + d^2 - 2*(l-l_C)*d*cos(phi - b20 - gam_D)) - l_0)^2

#kinetische Energie
T_1 = 1/2 * J_1*phiP^2 + 1/2 * m_1 * ((-l/2 * phiP * cos(phi))^2 + (zP + l/2 * phiP * sin(phi))^2)
T_2 = 1/2 * J_2*b20P^2 + 1/2 * m_2 * ((-l * phiP * cos(phi) - l_S2 * b20P * cos(b20 + gam_S2))^2 + (zP + l * phiP * sin(phi) + l_S2 * b20P * sin(b20 + gam_S2))^2)
T_3 = 1/2 * J_3*b30P^2 + 1/2 * m_3 * ((-L/2 * b30P * cos(b30))^2 + (zP + L/2 * b30P * sin(b30))^2)
T_4 = 1/2 * m_4 * zP^2

#Lagrange-Funktion L = T - V
@variables L_func
L_func = T_1 + T_2 + T_3 + T_4 - V_1 - V_2 - V_3 - V_4 - V_F

#generalisierte Kräfte
@variables Q_z, Q_phi, Q_b30, Q_b20
Q_phi = M + F2_Bz * l * sin(phi) - F2_Bx * l * cos(phi) 
Q_z = F2_Bz + F3_Bz
Q_b30 = -F3_Bx * L * cos(b30) + F3_Bz * L * sin(b30)
Q_b20 = F2_Bz * b * sin(b20 + gam_B) - F2_Bx * b * cos(b20 + gam_B)


#generalisierte Koordinaten ordnen
q = [phi, z, b30, b20]
q_formel = ["\\varphi", "z", "\\beta_{30}", "\\beta_{20}"]
qP = [phiP, zP, b30P, b20P]
qPP = [phiPP, zPP, b30PP, b20PP]
Q = [Q_phi, Q_z, Q_b30, Q_b20]

dL_dq = Vector{Any}(undef, length(q))
dL_dqP = Vector{Any}(undef, length(q))
d_dt_dL_dqP = Vector{Any}(undef, length(q))
dL_diff = Vector{Any}(undef, length(q))


#einmal alles ableiten
for i in 1:4
    dL_dq[i] = Symbolics.derivative(L_func, q[i])
    dL_dqP[i] = Symbolics.derivative(L_func, qP[i])
    d_dt_dL_dqP[i] = 0
    for j in 1:4
        d_dt_dL_dqP[i] += Symbolics.derivative(dL_dqP[i], q[j]) * qP[j] + Symbolics.derivative(dL_dqP[i], qP[j]) * qPP[j]
    end

    differenz = d_dt_dL_dqP[i] - dL_dq[i]
    dL_diff[i] = simplify(differenz)
end


#=
#probe: formeln drucken
ersetzungen = [
    "0.5" => "\\frac12",
    "0.25" => "\\frac14",

    "\\mathtt{phiPP}" => "\\ddot{\\varphi}",
    "\\mathtt{phiP}" => "\\dot{\\varphi}",
    "\\mathtt{phi}"   => "\\varphi",

    "\\mathtt{zPP}" => "\\ddot{z}",
    "\\mathtt{zP}" => "\\dot{z}",
    "\\mathtt{z}"   => "z",
    
    "\\mathtt{b30PP}" => "\\ddot{\\beta_{30}}",
    "\\mathtt{b30P}" => "\\dot{\\beta_{30}}",
    "\\mathtt{b30}"   => "\\beta_{30}",

    "\\mathtt{b20PP}" => "\\ddot{\\beta_{20}}",
    "\\mathtt{b20P}" => "\\dot{\\beta_{20}}",
    "\\mathtt{b20}"   => "\\beta_{20}",

    "\\mathtt{m\\_1}" => "m_1",
    "\\mathtt{m\\_2}" => "m_2",
    "\\mathtt{m\\_3}" => "m_3",
    "\\mathtt{m\\_4}" => "m_4",

    "\\mathtt{J\\_1}" => "J_1",
    "\\mathtt{J\\_2}" => "J_2",
    "\\mathtt{J\\_3}" => "J_3",

    "\\mathtt{gam\\_D}" => "\\gamma_D",
    "\\mathtt{gam\\_S2}" => "\\gamma_{S2}",
    "\\mathtt{gam\\_B}" => "\\gamma_B",

    "\\mathtt{l\\_S2}" => "l_{S2}",
    "\\mathtt{l\\_C}" => "l_C",
    "\\mathtt{l\\_0}" => "l_0",

    "\\mathtt{F2\\_Bx}" => "F^{(2)}_{Bx}",
    "\\mathtt{F2\\_Bz}" => "F^{(2)}_{Bz}",
    "\\mathtt{F3\\_Bx}" => "F^{(3)}_{Bx}",
    "\\mathtt{F3\\_Bz}" => "F^{(3)}_{Bz}"

]

println("############ FORMELN: ###############")
for i in 1:4
    print("\\large{\\textbf{Nach \$q_", i, " = ", q_formel[i], "\$:}}\n")
    string = "\\begin{dmath}\n" * latexify(Q[i], env=:raw) * " = " * latexify(dL_diff[i], env=:raw) * "\n\\end{dmath}\n\n\n"
    string_clean = reduce(replace, ersetzungen; init=string)
    print(string_clean)
end
=#

f_vektor = Vector{Any}(undef, length(q))
f_vektor[1] = simplify(expand(dL_diff[1] - M - (J_1 + 0.25(l^2)*m_1 + (l^2)*m_2)*phiPP - (0.5l*m_1*sin(phi) + l*m_2*sin(phi))*zPP - b20PP*(cos(phi)*cos(b20 + gam_S2) + sin(b20 + gam_S2)*sin(phi))*l*l_S2*m_2))
f_vektor[2] = simplify(expand(dL_diff[2] - zPP * (m_1 + m_2 + m_3 + m_4) - (1//2)*L*b30PP*m_3*sin(b30) - b20PP*l_S2*m_2*sin(b20 + gam_S2) - ((1//2)*l*m_1*sin(phi) + l*m_2*sin(phi))*phiPP))
f_vektor[3] = simplify(expand(dL_diff[3] - (J_3 + 0.25(L^2)*m_3)*b30PP - 0.5L*zPP*m_3*sin(b30)))
f_vektor[4] = simplify(expand(dL_diff[4] - zPP*l_S2*m_2*sin(b20 + gam_S2) - (cos(phi)*cos(b20 + gam_S2) + sin(b20 + gam_S2)*sin(phi))*l*l_S2*m_2*phiPP - b20PP*(J_2 + (l_S2^2)*m_2)))

M_q = Matrix{Any}(undef, length(q), length(q))
M_q[1,1] = (J_1 + 0.25(l^2)*m_1 + (l^2)*m_2)
M_q[1,2] = (0.5l*m_1*sin(phi) + l*m_2*sin(phi))
M_q[1,3] = 0
M_q[1,4] = (cos(phi)*cos(b20 + gam_S2) + sin(b20 + gam_S2)*sin(phi))*l*l_S2*m_2
M_q[2,1] = (1//2)*l*m_1*sin(phi) + l*m_2*sin(phi)
M_q[2,2] = m_1 + m_2 + m_3 + m_4
M_q[2,3] = (1//2)*L*m_3*sin(b30)
M_q[2,4] = l_S2*m_2*sin(b20 + gam_S2)
M_q[3,1] = 0
M_q[3,2] = 0.5L*m_3*sin(b30)
M_q[3,3] = J_3 + 0.25(L^2)*m_3
M_q[3,4] = 0
M_q[4,1] = (cos(phi)*cos(b20 + gam_S2) + sin(b20 + gam_S2)*sin(phi))*l*l_S2*m_2
M_q[4,2] = l_S2*m_2*sin(b20 + gam_S2)
M_q[4,3] = 0
M_q[4,4] = J_2 + (l_S2^2)*m_2


#=
print("##############\n")
for i in 1:4
    println(f_vektor[i], "\n\n")
end
print("##############\n")
=#

@variables F_Bx, F_Bz

Q_ZB = Vector{Any}(undef, length(q))
Q_ZB[1] = -F_Bz * l * sin(phi) + F_Bx * l * cos(phi)
Q_ZB[2] = 0
Q_ZB[3] = F_Bz * L * sin(b30) - F_Bx * L * cos(b30)
Q_ZB[4] = -F_Bz * b * sin(b20 + gam_B) + F_Bz * b * cos(b20 + gam_B)

@variables F, z_S

z_S = (1/(m_1 + m_2 + m_3 + m_4)) * (m_1 * (z - cos(phi) * l/2) + m_2 * (z - cos(phi) * l - cos(b20 + gam_S2) * l_S2) + m_3 * (z - a_z - cos(b30 * L/2)) + m_4 * (z - l_S4z))
z_SP = 0
z_SPP = 0

#=
for j in 1:4
    z_SP += Symbolics.derivative(z_S, q[j]) * qP[j]
end
z_SPP = 0
for j in 1:4
    z_SPP += Symbolics.derivative(z_SP, q[j]) * qP[j] + Symbolics.derivative(z_SP, qP[j]) * qPP[j]
end
=#
z_SP += Symbolics.derivative(z_S, q[1]) * qP[1]
z_SP += Symbolics.derivative(z_S, q[2]) * qP[2]
z_SP += Symbolics.derivative(z_S, q[3]) * qP[3]
z_SP += Symbolics.derivative(z_S, q[4]) * qP[4]

z_SPP += Symbolics.derivative(z_SP, q[1]) * qP[1] + Symbolics.derivative(z_SP, qP[1]) * qPP[1]
z_SPP += Symbolics.derivative(z_SP, q[2]) * qP[2] + Symbolics.derivative(z_SP, qP[2]) * qPP[2]
z_SPP += Symbolics.derivative(z_SP, q[3]) * qP[3] + Symbolics.derivative(z_SP, qP[3]) * qPP[3]
z_SPP += Symbolics.derivative(z_SP, q[4]) * qP[4] + Symbolics.derivative(z_SP, qP[4]) * qPP[4]

F = (m_1+m_2+m_3+m_4) * (z_SPP + g)


Q_ZF = Vector{Any}(undef, length(q))
Q_ZF[1] = F * l * sin(phi)
Q_ZF[2] = F
Q_ZF[3] = 0
Q_ZF[4] = F * l_F * sin(b20 + gam_F)

K_matrix = Matrix{Any}(undef, length(q)*2, length(q)*2)
for i in 1:length(q)*2
    for j in 1:length(q)*2
        K_matrix[i,j] = 0
    end
end
for i in 1:4
    K_matrix[i,i] = 1
    for j in 1:4
        K_matrix[i+4,j+4] = M_q[i,j]
    end
end

flug_vektor = Vector{Any}(undef, length(q)*2)
for i in 1:4
    flug_vektor[i] = qP[i]
    flug_vektor[i+4] = f_vektor[i] + Q_ZB[i]
end

kontakt_vektor = Vector{Any}(undef, length(q)*2)
for i in 1:4
    kontakt_vektor[i] = qP[i]
    kontakt_vektor[i+4] = f_vektor[i] + Q_ZB[i] + Q_ZF[i]
end

#Nebenbedingungen
gelenkB_x = -sin(phi) * l - sin(b20 + gam_B) * b ~ -a_x - sin(b20 + gam_B) * b
gelenkB_z = -cos(phi) * l - cos(b20 + gam_B) * b ~ -a_z - cos(b30) * L
bodenkontaktF = z - cos(phi) * l - cos(b20 + gam_F) * l_F ~ 0

println("\nProgrammende\n")