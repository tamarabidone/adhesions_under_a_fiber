package plos_comp_bio;


import java.util.ArrayList;

import java.util.Random;

/**
 * This program simulated brownian dynamics of integrin particles interactin with ligands 
 * 
 * @author tamarabidone & samuel campbell
 *
 * 
 */
public class Mechanosensing_2DApp {

	/* An ARRAY of integrin and ligands instances */
	
	ArrayList<Integrin_2D> sIntegrin = new ArrayList<Integrin_2D>();
	ArrayList<Ligand_2D> sLigand = new ArrayList<Ligand_2D>();

	double eta; // Viscosity
	double timestep;
	int nOfIntegrins;
	Random rand = new Random();
	double initialtime = 0;
	

	//INITIALIZE FUNCTION
	public void initialize() {

		timestep = 1E-4; //(s)
		nOfIntegrins = 300;
		
// INTEGRIN INITIALIZATION
		for (int i = 0; i < nOfIntegrins; i++) {
			sIntegrin.add(new Integrin_2D());
			sIntegrin.get(i).kon = 1; //(s-1) integrin activation rate
			sIntegrin.get(i).MyoForce=0; //(pN) force from actomyosin contractility
			sIntegrin.get(i).Pbundling =0.0; // probability of bundling. varies between 0 and 1
			sIntegrin.get(i).dt = timestep;
			sIntegrin.get(i).temperature = 300;
			sIntegrin.get(i).zetaI = 0.0142*5;// pN*s/um //friction corresponding to free diffusion of beta-3 integrin
			sIntegrin.get(i).DomainSize = 1; //(micrometer) lateral dimension of the square domain
			sIntegrin.get(i).BoundLigandYN = 0; // 0 for unbound; 1 for bound
			sIntegrin.get(i).TimesUnbound = 0; // sum how manu times i-th integrin unbinds a ligand
			sIntegrin.get(i).contaAT = 0; // counts for how many time steps the i-th integrin is attached
			sIntegrin.get(i).kub = 0; // (s-1) unbinding rate. initialized to zero
			sIntegrin.get(i).FinalLife = 0; // store the final lifetime before unbinding
			sIntegrin.get(i).V = 0.0; //actin flow velocity
			
			// x and y coordinates are chosen randomly within the square domain centered at 0
			double[] initialpos = new double[2];
			initialpos[0] = -sIntegrin.get(i).DomainSize / 2 + Math.random() * (sIntegrin.get(i).DomainSize);
			initialpos[1] = -sIntegrin.get(i).DomainSize / 2 + Math.random() * (sIntegrin.get(i).DomainSize);
			sIntegrin.get(i).start(initialpos);
		}


// LIGANDS INITIALIZATION	
			double x=-0.5-0.05;
			double y;
			int i=-1;
			for (int r = 0; r < 21; r++) { //create a lattice of 20 cells with side 0.05 micrometers
				x=x+0.05;
		        y=-0.5-0.05;
				for (int c = 0; c < 21; c++) {
					y=y+0.05;
					double[] initialpos = new double[2];
					i=i+1;
					sLigand.add(new Ligand_2D());
					sLigand.get(i).k =1.4*1000;  //(pN/micrometer) substrate spring constant
					sLigand.get(i).DomainSize = 1.0; //(micrometers)
					sLigand.get(i).dt = timestep; //(s) time step
					sLigand.get(i).integrinBoundNum = -1; //store the ID of the bound integrin. -1 when no integrin is bound
					sLigand.get(i).zetaL = 800;
					  
			initialpos[0] = x;
			initialpos[1] = y;
			sLigand.get(i).start(initialpos);
		}
				
				;
	}
	}



	/* STEP FUNCTION
	 * The step is done at each time step */
	public void doStep() {

		/* Initiate the class with in interactions between integrins and ligands */
		Force_Interaction forceinter = new Force_Interaction(sIntegrin, sLigand, nOfIntegrins, 441);
		Print_output print = new Print_output (sIntegrin, sLigand, nOfIntegrins, 441);
		
		for (int j = 0; j < 441; j++)
			sLigand.get(j).step(); //maintain ligands in fixed positions

		for (int j = 0; j < nOfIntegrins; j++) 
			sIntegrin.get(j).step(); // impose thermal force and update positions
	
		forceinter.CatchBond(); //calculate force on integrin-ligand bonds and calculate koff. update bound or unbound state
		forceinter.IL_interact(); //form new bonds between free integrins and ligands
		
		print.print_all(); // print output into files

	}


	public static void main(String[] args) {
		Mechanosensing_2DApp simulation = new Mechanosensing_2DApp();
		simulation.initialize();
		for (int i = 0; i < (int) (500 / simulation.timestep); i++) {
			simulation.doStep();
		}

		System.exit(0);
	}
}

