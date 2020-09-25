# Replication Instructions

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

### Example Replication Ιnstructions. (require all the above installed)

1) Download the Repeatability Evaluation Package (REP) from:
	https://github.com/janis10/cis2m.
2) Copy the downloaded folder to your MATLAB work directory.
3) Open MATLAB, and while at the directory where `cis2m-master` is copied, type the following in MATLAB command-line to add all the necessary folders into the MATLAB path:
    1. `addpath(genpath('./cis2m-master/algorithms'))`		\hit enter.
    2. `addpath(genpath('./cis2m-master/support_functions'))`		\hit enter.
    3. `addpath(genpath('./cis2m-master/paper-archive/ACC21'))`		\hit enter.
	
4) To use MOSEK with MPT3 for solving Linear Programs (optional) type in the following in MATLAB command-line:
	
	`mptopt('lpsolver', 'MOSEK')`		\hit enter.
	
	To make sure MPT uses MOSEK type in `mptopt` and hit enter. MOSEK should appear as the selected LP solver. Alternatively, by typing in `mpt_init` and hitting enter, MPT decides which solver to use from the available solvers in the MATLAB path.

#### Instructions for each specific example:
* Section IV- Tables I & II.
    1) In MATLAB:
    
	* Type in MATLAB command-line: `example1`	\hit enter.

	* Successfully running this script file generates the values for Tables I & II. The appropriate number of trailers (columns of Tables I & II) are selected by choosing the value of `N` in line 48. 
	
	* Type in MATLAB command-line: `example1_mpt`	\hit enter.
	
	* Successfully running this script file generates the values for Table II, row of MPT3.
	
    2) In Julia/Jupyter Notebook:
	* If Jupyter Notebook is not open:
		Open Julia and type in the following commands:
    
			* using IJulia	\hit enter
			* notebook()		\hit enter
      
		This opens Jupyter Notebook in your browser.
    
	* Go to Jupyter Notebook and navigate to 
		`/cis2m-master/ACC21/paper-examples/example1/`.

	* Open notebook `example1_ellips_times.ipynb` by clicking on it.
	From the menu bar click on `Cell`, and then on `Run All`.
	
	* Upon successful termination the output of each cell displays the times (in seconds) corresponding
	to the times reported in: Table II, row of "Alg. [20]".

	* Note 1: sometimes due to the initial 'using' commands in the first cell, 
		the time corresponding to the first cell is much larger. By going back to the menu bar
		click on `Cell`, and then on `Run All`, one obtains the reported times.
	
	
* Section IV- Tables III & IV.

	* Type in MATLAB command-line: `example2_disturb`	\hit enter.

	* Successfully running this script file generates the values for Table III. The disturbance bounds `|r_d|` (columns of Table III) are selected by choosing the values of `param.rd_min` and `param.rd_max` in lines 36 and 37 respectively of the MATLAB file `/cis2m-master/ACC21/paper-examples/example2/constant_lk4.m`. 

	* Type in MATLAB command-line: `example2`	\hit enter.

	* Successfully running this script file generates the values for Table IV.
