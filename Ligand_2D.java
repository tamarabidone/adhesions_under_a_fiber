package plos_comp_bio;

import java.util.Random;

public class Ligand_2D {

	double time = 0, dt, k;
	double Bead_PositionL[];  
	double initialp[], substrate[]=new double[3];	//Initial positions
	double radius, FF_interact, zetaL, DomainSize;
	int  integrinBoundNum,  BoundYN;
	
    Random rand = new Random();
    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_PositionL=new double[2];
	    Bead_PositionL[0]=initialp[0];
	    Bead_PositionL[1]=initialp[1];
	   substrate[0]=initialp[0];
	   substrate [1]=initialp[1];
	 
	
	  }
    /*Keep tracks of time */
	public void step(){		
		time+=dt;	
		

    }
}
