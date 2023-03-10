# adhesions_under_a_fiber
This code simulated the assembly and elongation of integrin adhesions under an actin fiber. 
Input variables are myosin force and probability of actin bundling. 
Outputs are percentage of ligated integrins, Integrin-ligand bond lifetimes and adhesion orientation


The inputs can be changed in the config file, whihc needs to be in the same directory of the prohram
The outputs are generated from the print_output.java file

Mechanosensing_2DApp.java contains:
initialize()
main()

It first initializes positions of integrins and ligands in the domain; then performs code iteration throug the step function
The step() function, applies force on itegrins and displaces them. Forces result from thermal effects, substrate stiffness (when integrins are bound to a ligand), and actomyosin contractility (when integrins are ligated under the fiber region.
At the end of the step() function, displacements are computed and integrin states are updated (ligand-bound or unligated)

