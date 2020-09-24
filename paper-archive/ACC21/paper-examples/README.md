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
    1. `addpath(genpath('./cis2m-master/support_functions'))`		\hit enter.
    2. `addpath(genpath('./cis2m-master/HSCC20'))`		\hit enter.
	
4) To use MOSEK with MPT3 for solving Linear Programs (optional) type in the following in MATLAB command-line:
	
	`mptopt('lpsolver', 'MOSEK')`		\hit enter.
	
	To make sure MPT uses MOSEK type in `mptopt` and hit enter. MOSEK should appear as the selected LP solver. Alternatively, by typing in `mpt_init` and hitting enter, MPT decides which solver to use from the available solvers in the MATLAB path.

#### Instructions for each specific example:
* Section 2 - Figure 1. AND Section 5 - Table 1.

	* Type in MATLAB command-line: `example_1_2`	\hit enter.

	* Successfully running this script file generates Fig. 1a-1d.
	
	* In addition: 
		- resulting variable `Times` contains the values reported in
	Table 1, row 1 (page 6, bottom-right, Section 5).
		- resulting variable `VolumesPercentage` contains the values reported in
	Table 1, row 2 (page 6, bottom-right, Section 5).

	* Saves in `mat2julia.mat` the matrices corresponding to the computed controlled invariant sets
	to a `.mat` file in order to export them, and use them in Julia. 
	(commented out currently due to importing issues on Julia's end; we have manually copied the
	required matrices in the corresponding notebooks)

* Section 5 - Figure 2.
	
  	* Open Julia and type in the following commands:
  
		* using IJulia	\hit enter
		* notebook()		\hit enter
    
		This opens Jupyter Notebook in your browser.
  
	* Navigate to `/cis2m-master/HSCC20/paper-examples/example_1_2/`.

	* Open notebook `example_2_figure_5.ipynb` by clicking on it.
	From the menu bar click on `Cell`, and then on `Run All`.
	
	Upon successful termination Figure 5 is generated.

* Section 5 - Table 2.
    1) In MATLAB:
  
	 * Type in MATLAB command-line: `example_3`	\hit enter.

	 * Resulting variable `Times` contains the values reported in Table 2, row 1 (Time) 
	  (page 7, top-right, Section 5).

	 * Resulting variable `Volumes` contains the values reported in Table 2, row 1 (Volume) 
  	  (page 7, top-right, Section 5).

	 * Resulting variable `VolElls` contains the values reported in Table 2, row 2 (Volume) 
	  (page 7, top-right, Section 5).
  
	Note 1: we have opted out of computing the volume for N=4 here since its computation exceeded 5 hours. 
	
	Note 2: the values of matrices used to compute the contents of 'VolElls' can be found in the
	following notebook:
	
    2) In Julia/Jupyter Notebook:
	* If Jupyter Notebook is not open:
		Open Julia and type in the following commands:
    
			* using IJulia	\hit enter
			* notebook()		\hit enter
      
		This opens Jupyter Notebook in your browser.
    
	* Go to Jupyter Notebook and navigate to 
		`/cis2m-master/HSCC20/paper-examples/example_3/`.

	* Open notebook `example_3_ellips_times.ipynb` by clicking on it.
	From the menu bar click on `Cell`, and then on `Run All`.
	
	* Upon successful termination the output of each cell displays the times (in seconds) corresponding
	to the times reported in: Table 2, row 2 (Times) (page 7, top-right, Section 5).

	* Note 1: sometimes due to the initial 'using' commands in the first cell, 
		the time corresponding to the first cell is much larger. By going back to the menu bar
		click on `Cell`, and then on `Run All`, one obtains the reported times.
	* Note 2: the values of `P1`, `P2`, `P3` matrices in `example_3.m` can be found in the outputs of each cell.
	
	* Note 3: we have opted out of computing the result for `N=4` in this package, since using the interior point method as for `N=1,2,3` results in error, and using the symmetric method that computes an ellipsoidal set around the origin results in a negligibly small set since the origin lies on a facet of the safe set.
	
    3) In MATLAB:
  
	 * Type in MATLAB command-line: `example_3_mpt`	\hit enter.

	 * Resulting variable `TimeMCIS` contains the values reported in Table 2, row 3 (Time) 
	  (page 7, top-right, Section 5). See note below.
  
	 * Note: For `n=7,9` we did not compute the maximal controlled invariant set: For `n=7`, after 50 iterations and 363.2 seconds the algorithm did not converge. For `n=9` computation was aborted after 2 hours without completing even the 50 iterations. Thus, we do not compute the corresponding volumes.
	
* Section 5 - Table 3.

	* Type in MATLAB command-line: `example3_OuterApprox_N1n3`	\hit enter.

		- The values reported in Table 3 for `n=3` are displayed. 
	
	* Type in MATLAB command-line: `example3_OuterApprox_N2n5`	\hit enter.

	    - The values reported in Table 3 for `n=5` are displayed. 

	* Note: depending on the machine that the code runs on, the last case `(n=5, d=8)` might take extremely long time to run. This is partially due to the command in line 104 `v1 = replace(v,x,f)`. On these grounds, we measure and report in Table 3 only the time to solve the optimization problem (lines 162-164), and not formulating it, since there might be a faster way to implement the formulation. As reported in the paper and in the corresponding MATLAB files, the main code formulating the optimization problem was obtained from https://homepages.laas.fr/henrion/ (under Software -> ROA), and adapted to our example.
