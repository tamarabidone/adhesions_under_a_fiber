package nov_15;

import java.util.*;

/**
 * The class Integrin_3D includes the properties of a single
 * integrin. Thermal and velocity motions 
 * are included. Main class is in Mechanosensing_3DApp.
 * @author tamarabidone
 *
 */

public class Integrin_3D{
	double V;
	double time = 0, dt;
	double Bead_Position[];  
	double Bead_Force[];	
	double initialp[];	//Initial positions
	double kub=0;
	int LigandNum, TimesUnbound, MyoForce;	
	int contaAT, BoundLigandYN;
	double zetaI,  contaBoundLifetimeTot, FinalLife, finalTension=0;
	double DomainSize;
	double temperature;
    Random rand = new Random();
    double zIntegrin;
    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_Position=new double[3];
	    Bead_Position[0]=initialpos[0];
	    Bead_Position[1]=initialpos[1];
	    Bead_Position[2]=initialpos[2];   
    }
    

    /*Update forces and move along the 2D surface*/
	public void step(){		
	    Bead_Force=new double[3];
		time+=dt;
		Boundary_Force_Update();
		if (BoundLigandYN==1 && LigandNum>=189 && LigandNum<252)
			Bead_Force[1]+=MyoForce; // only the ligand-bound integrins in the fiber region are subjected to force
		
		if (BoundLigandYN==0){
			Thermal_Force_Update();
			kub=0; //when integrin is unbound from a ligand kub is set to 0, since it depends on substrate force and there is no force from substrate
			}
		else{
		Thermal_Force_Update();
		Flow_Force_Update();
		}

		move();
    }

	
/* Thermal force is applied in 2D */
	protected void Thermal_Force_Update(){
			double thermal_Force=Math.sqrt(2*0.0000138*temperature*zetaI/dt);
			Bead_Force[0]+=thermal_Force*rand.nextGaussian();
			Bead_Force[1]+=thermal_Force*rand.nextGaussian();
		
	}	

	protected void Flow_Force_Update(){
		Bead_Force[1]+=-zetaI*V;
		
	}

	protected void Boundary_Force_Update(){
		if (Bead_Position[0]>=DomainSize/2 && BoundLigandYN==0){
			Bead_Position[0]=Bead_Position[0]-DomainSize;
		}
		if (Bead_Position[0]<-DomainSize/2 && BoundLigandYN==0){
			Bead_Position[0]=Bead_Position[0]+DomainSize;
		}
		
		if (Bead_Position[1]>=DomainSize/2  && BoundLigandYN==0){
			Bead_Position[1]=Bead_Position[1]-DomainSize;
		}
		if (Bead_Position[1]<-DomainSize/2  && BoundLigandYN==0){
			Bead_Position[1]=Bead_Position[1]+DomainSize;
		}
	}
	
	protected void move(){ // the move function works in 2D
			for(int j=0; j<2;j++){
				Bead_Position[j]=Bead_Position[j]-Bead_Force[j]*dt/zetaI;/// Check this--> sign is correct?
		}
	}
	


}
