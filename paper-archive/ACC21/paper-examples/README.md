# Replication Instructions
This folder contains the examples presented in:
T. Anevlavis, Z. Liu, N. Ozay and P. Tabuada, "An enhanced hierarchy for (robust) controlled invariance", In 2021 American Control Conference (ACC). (Accepted)

### Requirements:
	* MATLAB (written and tested in R2018a or later).
	* Parallel Computation Toolbox (optional for improved performance).
	* Multi-Parametric Toolbox (MPT) for MATLAB, available at:
		https://www.mpt3.org
	* YALMIP:
		included in MPT3 above.
	* MOSEK v8.1.0.75 (or newer) available at:
		https://www.mosek.com/downloads/8.1.0.75/
	* Julia 1.0.5, available at:
		https://julialang.org/downloads/
	* Jupyter Notebook.
	* Tested in: macOS v10.14 or higher, and Windows 10.

### Installation:
	* For MPT the instructions at:
      			https://www.mpt3.org/Main/Installation
      	suffice.
	* For MOSEK the simple instructions at:
			https://docs.mosek.com/8.1/toolbox/install-interface.html
		suffice. Make sure to get a free license and add it before using MOSEK. 
		A trial license can be found here:
		https://www.mosek.com/products/trial/
		and an academic license here:
		https://www.mosek.com/products/academic-licenses/. 
		Simple instructions to place the license into the correct directory can be found here: 
		https://docs.mosek.com/9.1/licensing/quickstart.html.
	* For Julia:
		1) Download and install the executable file depending on your OS.
		2) Open Julia.
		3) Type the following commands:
			* import(Pkg)                     \hit enter
			* Pkg.add("IJulia")               \hit enter
			* Pkg.add("SwitchOnSafety")       \hit enter
			* Pkg.add("JuMP")                 \hit enter
			* Pkg.add("Polyhedra")            \hit enter
			* Pkg.add("MathematicalSystems")  \hit enter
			* Pkg.add("MosekTools")           \hit enter
			* Pkg.add("LinearAlgebra")        \hit enter
			* Pkg.add("Plots")                \hit enter
		   This adds all the required packages in Julia.
		4) Type:
			using IJulia	\hit enter
			notebook()	\hit enter
		   Here you might be prompted to install Conda and Jupyter.
		   We recommend that you do so, as to smoothly replicate the results.
		   Successfully completing this step will open Jupyter Notebook in your browser.
