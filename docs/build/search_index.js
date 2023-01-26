var documenterSearchIndex = {"docs":
[{"location":"Functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"Functions/#Clausen-functions","page":"Functions","title":"Clausen functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"The Clausen function is given by","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"Cl_2(phi)=-int_0^phi log2sin(x2) dx","category":"page"},{"location":"Functions/#Examples","page":"Functions","title":"Examples","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"using Plots, FewSpecialFunctions, LaTeXStrings\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,15,1000)\nplot(x,Clausen.(x), label=L\"Cl_2(ϕ)\")\nxlabel!(L\"ϕ\")\ntitle!(\"Clausen function\")","category":"page"},{"location":"Functions/#Debye-functions","page":"Functions","title":"Debye functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"using Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,25,1000)\n\nplot(x,Debye_function.(1,x),label=L\"D_1(x)\")\nplot!(x,Debye_function.(2,x),label=L\"D_2(x)\")\nplot!(x,Debye_function.(3,x), label=L\"D_3(x)\")\ntitle!(\"Debye Functions\")\nxlabel!(L\"x\")","category":"page"},{"location":"Functions/#Regular-Coulomb-wave-functions","page":"Functions","title":"Regular Coulomb wave functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"This implementation need the gamma function from SpecialFunctions.jl","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,25,1000)\n\nplot(x,regular_coulomb.(0,0.3,x), label=L\"F_0(0.3,ρ)\")\nplot!(x,regular_coulomb.(0,-0.3,x), label=L\"F_0(0.3,ρ)\")\nxlabel!(L\"ρ\")\ntitle!(\"Regular Coulomb Wave Functions\")","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"And ","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,25,1000)\n\nplot(x,regular_coulomb.(1e-5,5,x), label=L\"F_0(5,ρ)\")\nplot!(x,regular_coulomb.(1,5,x), label=L\"F_1(5,ρ)\")\nplot!(x,regular_coulomb.(2,5,x), label=L\"F_2(5,ρ)\")\nplot!(x,regular_coulomb.(3,5,x), label=L\"F_3(5,ρ)\")\ntitle!(\"Regular Coulomb Wave Functions\")\nxlabel!(L\"ρ\")","category":"page"},{"location":"Functions/#Struve-functions","page":"Functions","title":"Struve functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\n\nx = range(0,15,1000)\n\nplot(x,Struve.(0,x),label=L\"H_0(x)\")\nplot!(x,Struve.(1,x),label=L\"H_1(x)\")\nplot!(x,Struve.(2,x),label=L\"H_3(x)\")\nxlabel!(L\"x\")\ntitle!(\"Struve Functions\")","category":"page"},{"location":"Functions/#Fresnel-functions","page":"Functions","title":"Fresnel functions","text":"","category":"section"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(-25,25,5000)\n\nplot(x,Fresnel_C_integral.(x),label=L\"C(x)\")\nplot!(x,Fresnel_C_err.(x), ls=:dash, lw=1.5, label=L\"\\tilde{C}(x)\")\ntitle!(\"Fresnel Integral\")\nxlabel!(L\"x\")\n","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"and","category":"page"},{"location":"Functions/","page":"Functions","title":"Functions","text":"\nusing Plots, FewSpecialFunctions\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(-25,25,5000)\n\nplot(x,Fresnel_S_integral.(x),label=L\"S(x)\")\nplot!(x,Fresnel_S_err.(x), ls=:dash, lw=1.5, label=L\"\\tilde{S}(x)\")\ntitle!(\"Fresnel Integral\")\nxlabel!(L\"x\")\n","category":"page"},{"location":"#FewSpecialFunctions.jl","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"","category":"section"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"Some special functions.","category":"page"},{"location":"#Overview","page":"FewSpecialFunctions.jl","title":"Overview","text":"","category":"section"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"FewSpecialFunctions.jl provides implementations of a few special functions. So far this includes","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"Clausen functions\nDebye functions\nRegular Coulomb wave functions\nStruve functions\nFresnel functions","category":"page"},{"location":"#Installation","page":"FewSpecialFunctions.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"Get the latest stable release with Julia's package manager:","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"julia ] add FewSpecialFunctions","category":"page"},{"location":"#Quick-example","page":"FewSpecialFunctions.jl","title":"Quick example","text":"","category":"section"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"Here is how to generate an Euler spiral using FewSpecialFunctions.jl. ","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"using Plots, FewSpecialFunctions, LaTeXStrings\n\nplot_font = \"Computer Modern\"\ndefault(\n    fontfamily=plot_font,\n    linewidth=2.5, \n    framestyle=:box, \n    label=nothing, \n    grid=true,\n    palette=:tab10,\n)\n\nx = range(-25,25,5000)\n\nplot(Fresnel_C_err.(x),Fresnel_S_err.(x))\nxlabel!(L\"C(x)\")\nylabel!(L\"S(x)\")\ntitle!(\"Euler Spiral\")","category":"page"},{"location":"#About","page":"FewSpecialFunctions.jl","title":"About","text":"","category":"section"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"MIT License","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"Copyright (c) 2023 Martin Mikkelsen","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.","category":"page"},{"location":"","page":"FewSpecialFunctions.jl","title":"FewSpecialFunctions.jl","text":"THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"}]
}
