# adhesions_under_a_fiber
This code simulated the assembly and elongation of integrin adhesions under an actin fiber. 
Input variables are myosin force and probability of actin bundling. 
Outputs are percentage of ligated integrins, Integrin-ligand bond lifetimes and adhesion orientation


The inputs can be changed in the Mechanosensing_3DApp.java file
The outputs are generated from the same file, in the second part of the step() function

Mechanosensing_3DApp.java contains the main(). It first initializes positions of integrins and ligands in the domain; then performs code iteration throug the step( function)
The step() function, fist applies parallel force on itegrins and displaces them, second it applies substrate force and displaces them one more time.
At the end of the step() function integrin states are updated (ligand-bound or unligated)

