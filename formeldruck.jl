#dieser code kann in das julia-skript eingefügt werden,
#um gleichungen in latex-form auszugeben


#Probe: Formeln drucken
#um die Formelzeichen in LaTex richtig darzustellen
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
    string = "\\begin{dmath}\n" * latexify(simplify(expand(Q[i])), env=:raw) * " = " * latexify(dL_diff[i], env=:raw) * "\n\\end{dmath}\n\n\n"
    string_clean = reduce(replace, ersetzungen; init=string)
    print(string_clean)
end