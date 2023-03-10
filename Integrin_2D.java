package plos_comp_bio;

import java.util.*;

/**
 * The class Integrin_3D includes the properties of a single
 * integrin. Thermal and velocity motions 
 * are included. Main class is in Mechanosensing_3DApp.
 * @author tamarabidone
 *
 */

public class Integrin_2D{
	double V;
	double time = 0, dt;
	double Bead_Position[];  
	double Bead_Force[];	
	double initialp[];	//Initial positions
	double kub=0, kon=0;
    int LigandNum, TimesUnbound;	
	int contaAT, BoundLigandYN;
	double zetaI,  contaBoundLifetimeTot, FinalLife, finalTension=0, Pbundling,  MyoForce;
	double DomainSize;
	double temperature;
    Random rand = new Random();
   
    public void start(double[] initialpos){
    	initialp = initialpos;
	    Bead_Position=new double[2];
	    Bead_Position[0]=initialpos[0];
	    Bead_Position[1]=initialpos[1];
	   
    }
    

    /*Update forces and move along the 2D surface*/
	public void step(){		
	    Bead_Force=new double[3];
		time+=dt;
		
		
		 
		Thermal_Force_Update();
		
		if (BoundLigandYN==0){
			kub=0;

			Boundary_Force_Update();
			
			// check fiber region
			  if(Bead_Position[1]<0.15 && Bead_Position[1]>-0.15 && Math.random()<Pbundling) {
				  kon=3; 
			  } 
			  else 
			  {
				  kon=1; 
			  }
			 
			
			}
		
		else {
			
			kon=1;
			
			  if(Bead_Position[1]<0.15 && Bead_Position[1]>-0.15) 
			  {
			  
				  Myo_Force_Update(); 
			  }
			 
		

		}
		
		move();
		
		    }


	protected void Thermal_Force_Update(){
			double thermal_Force=Math.sqrt(2*0.0000138*temperature*zetaI/dt);
			Bead_Force[0]+=thermal_Force*rand.nextGaussian();
			Bead_Force[1]+=thermal_Force*rand.nextGaussian();
		}	


	protected void Myo_Force_Update(){
		Bead_Force[0]+=MyoForce;
		
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
	
	protected void move(){
			for(int j=0; j<2;j++){
				Bead_Position[j]=Bead_Position[j]-Bead_Force[j]*dt/zetaI;
		}
	}
	
	

}
