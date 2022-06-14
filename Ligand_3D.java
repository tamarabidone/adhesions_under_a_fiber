package nov_15;

import java.util.Random;

public class Ligand_3D {

	double time = 0, dt;
	double Bead_PositionL[];  
	double initialp[], substrate[]=new double[3];	//Initial positions
	double radius, FF_interact, zetaL, DomainSize, zLigand;
	int  integrinBoundNum, Y, BoundYN;
	
    Random rand = new Random();
    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_PositionL=new double[3];
	    Bead_PositionL[0]=initialp[0];
	    Bead_PositionL[1]=initialp[1];
	    Bead_PositionL[2]=initialp[2];   
	   substrate[0]=initialp[0];
	   substrate [1]=initialp[1];
	    substrate[2]=initialp[2]-Math.abs(initialp[2]);
	
	  }
    /*Keep tracks for time and re-assign initial positions to the ligand, do make sure it does not move*/
	public void step(){		
		time+=dt;	substrate[0]=initialp[0];
	    substrate [1]=initialp[1];
	    substrate[2]=initialp[2]-Math.abs(initialp[2]);
	 Bead_PositionL[0]=initialp[0];
	Bead_PositionL[1]=initialp[1];
   Bead_PositionL[2]=initialp[2];  
    }
}
