var documenterSearchIndex = {"docs":
[{"location":"Functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"Functions/#Clausen-functions","page":"Functions","title":"Clausen functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"The Clausen function is given by","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Cl_2(phi)=-int_0^phi log2sin(x2)  textdx","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"using Plots, FewSpecialFunctions, LaTeXStrings\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,15,1000)\nplot(x,Clausen.(x), label=L\"Cl_2(ϕ)\")\nxlabel!(L\"ϕ\")\ntitle!(\"Clausen function\")","category":"page"},{"location":"Functions/#Debye-functions","page":"Functions","title":"Debye functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"The Debye functions are given by","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"    D_n(x)= fracnx^n int_0^x fract^ntexte^t-1  textdx","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"And","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"using Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(0,25,1000)\n\nplot(x,Debye_function.(1,x),label=L\"D_1(x)\")\nplot!(x,Debye_function.(2,x),label=L\"D_2(x)\")\nplot!(x,Debye_function.(3,x), label=L\"D_3(x)\")\ntitle!(\"Debye Functions\")\nxlabel!(L\"x\")","category":"page"},{"location":"Functions/#Regular-Coulomb-wave-functions","page":"Functions","title":"Regular Coulomb wave functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"The Coulomb wave equation for a charged particle with arbitrary angular momentum and charge is given by ","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"    nabla^2psi +left( k^2-frac2muhbar^2V(r)right)psi = 0","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"where mu is the reduced mass of the system. The radial wave function u(r) satisfies the following differential equation","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\tfractextd^2 u_elltextdr^2+left( k^2-fracell(ell+1)r^2-frac2muhbar^2fracZe^2rright)u_ell=0","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"where Z is the product of the charges. Two independent solutions can be found to this equation – these are called the regular and irregular Coulomb wave functions denoted F_ell(r) and G_ell(r) respectively. The regular Coulomb wave function F_ell(r) is a real function that vanishes at r=0 and the behaviour of the function is described using a parameter eta which describes how strongly the Coulomb interaction is","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\teta = fracZmcalpha hbar k","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"where m is the mass of the particle, k is the wave number and alpha is the fine structure constant. The solution to is given by","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\tF_ell(etakr) = C_ell (eta) (kr)^ell+1texte^-ikr  _1 F_1(ell+1-ieta2ell+22ikr)","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"where _1F_1(kr) is a confluent hypergeometric function and C_ell(eta) is a normalization constant given by ","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\tC_ell(eta) = frac2^ell texte^-pieta2Gamma(ell+1+ieta)(2ell+1)","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"where Gamma is the gamma function. For numerical purposes, it is useful to use the integral representation of the regular Coulomb wave function","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\tF_ell(etarho) = fracrho^ell+12^ell e^irho-(pieta2)Gamma(ell+1+ieta) int_0^1 e^-2irho tt^ell+ieta(1-t)^ell-ieta  textdt","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"This implementation need the gamma function from SpecialFunctions.jl","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,25,1000)\n\nplot(x,regular_coulomb.(0,0.3,x), label=L\"F_0(0.3,ρ)\")\nplot!(x,regular_coulomb.(0,-0.3,x), label=L\"F_0(-0.3,ρ)\")\nxlabel!(L\"ρ\")\ntitle!(\"Regular Coulomb Wave Functions\")","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Use a similar approach to plot the regular Coulomb functions for different a ell","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,25,1000)\n\nplot(x,regular_coulomb.(1e-5,5,x), label=L\"F_0(5,ρ)\")\nplot!(x,regular_coulomb.(1,5,x), label=L\"F_1(5,ρ)\")\nplot!(x,regular_coulomb.(2,5,x), label=L\"F_2(5,ρ)\")\nplot!(x,regular_coulomb.(3,5,x), label=L\"F_3(5,ρ)\")\ntitle!(\"Regular Coulomb Wave Functions\")\nxlabel!(L\"ρ\")","category":"page"},{"location":"Functions/#Struve-functions","page":"Functions","title":"Struve functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"The Struve functions are solutions of the non-homogeneous Bessel's differential equation","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"    x^2 fractextd^2 ytextdx^2 + x fractextdytextdx+(x^2-alpha^2)y = frac4(x2)^alpha+1sqrtpiGamma(alpha+12)","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"The Struve functions are implemented using the following integral representation","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"    mathbfH_nu(z) = frac2(z2)^nusqrtpiGamma(nu+12) int_0^1 (1-t)^nu-12sin(zt)  textdt","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"And","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"    mathbfH_nu(z) = frac2(z2)^nusqrtpiGamma(nu+12) int_0^pi2 sin(zcos(theta)) sin^2nu(theta)  textdtheta","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Here is an example","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(-5,5,1000)\n\nplot(x,Struve.(0,x),label=L\"H_0(x)\")\nplot!(x,Struve.(1,x),label=L\"H_1(x)\")\nplot!(x,Struve.(2,x),label=L\"H_2(x)\")\nplot!(x,Struve.(3,x),label=L\"H_3(x)\")\nplot!(x,Struve.(4,x),label=L\"H_4(x)\")\nplot!(x,Struve.(5,x),label=L\"H_5(x)\")\nxlabel!(L\"x\")\ntitle!(\"Struve Functions\")","category":"page"},{"location":"Functions/#Fresnel-functions","page":"Functions","title":"Fresnel functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(-25,25,5000)\n\nplot(x,Fresnel_C_integral.(x),label=L\"C(x)\")\nplot!(x,Fresnel_C_err.(x), ls=:dash, lw=1.5, label=L\"\\tilde{C}(x)\")\ntitle!(\"Fresnel Integral\")\nxlabel!(L\"x\")\n","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"and","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(-25,25,5000)\n\nplot(x,Fresnel_S_integral.(x),label=L\"S(x)\")\nplot!(x,Fresnel_S_err.(x), ls=:dash, lw=1.5, label=L\"\\tilde{S}(x)\")\ntitle!(\"Fresnel Integral\")\nxlabel!(L\"x\")","category":"page"},{"location":"Functions/#Benchmarks","page":"Functions","title":"Benchmarks","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Using the integral approach","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing FewSpecialFunctions,BenchmarkTools\n\n@benchmark Fresnel_C_err.($x)\n","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Using the error function","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing FewSpecialFunctions,BenchmarkTools\n\n@benchmark Fresnel_C_integral.($x)","category":"page"},{"location":"#FewSpecialFunctions.jl","page":"Home","title":"FewSpecialFunctions.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Some special functions.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"FewSpecialFunctions.jl provides implementations of a few special functions. So far this includes","category":"page"},{"location":"","page":"Home","title":"Home","text":"Clausen functions\nDebye functions\nRegular Coulomb wave functions\nStruve functions\nFresnel functions","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Get the latest stable release with Julia's package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia ] add FewSpecialFunctions","category":"page"},{"location":"#Quick-example","page":"Home","title":"Quick example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Here is how to generate an Euler spiral using FewSpecialFunctions.jl. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Plots, FewSpecialFunctions, LaTeXStrings\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(-25,25,5000)\n\nplot(Fresnel_C_err.(x),Fresnel_S_err.(x))\nxlabel!(L\"C(x)\")\nylabel!(L\"S(x)\")\ntitle!(\"Euler Spiral\")","category":"page"}]
}