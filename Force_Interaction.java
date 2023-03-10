package plos_comp_bio;
import java.util.ArrayList;

public class Force_Interaction{
	
	ArrayList<Integrin_2D> sIntegrin = new ArrayList<Integrin_2D>();
	ArrayList<Ligand_2D> sLigand = new ArrayList<Ligand_2D>();
	int nOfIntegrins,nOfLigands;

	public Force_Interaction(ArrayList<Integrin_2D> sIntegrin,  ArrayList<Ligand_2D> sLigand, int nOfIntegrins, int nOfLigands){
		this.sIntegrin= sIntegrin;
		this.sLigand= sLigand;
		this.nOfIntegrins = nOfIntegrins;	
		this.nOfLigands = nOfLigands;	
	}

	public void IL_interact(){
		for(int i=0;i<nOfIntegrins;i++){		
			double r_min=10000;
		    double k_it=sLigand.get(0).k;
	
			if (sIntegrin.get(i).BoundLigandYN==1){
				int lig1=sIntegrin.get(i).LigandNum;
				double rr=distance(sIntegrin.get(i).Bead_Position,sLigand.get(lig1).Bead_PositionL);
				double force=k_it*( rr-(0.0001));
				sLigand.get(lig1).FF_interact=force;
				sIntegrin.get(i).Bead_Position=move(sIntegrin.get(i).Bead_Position, sLigand.get(lig1).Bead_PositionL, force, rr, sIntegrin.get(i).dt, sIntegrin.get(i).zetaI);
				sLigand.get(lig1).Bead_PositionL=move(sLigand.get(lig1).Bead_PositionL, sIntegrin.get(i).Bead_Position,  force, rr, sLigand.get(lig1).dt, sLigand.get(lig1).zetaL);	
				sLigand.get(lig1).BoundYN=1;
				sLigand.get(lig1).integrinBoundNum=i;
			}
			
			int lig=-1;
			r_min=1000;
			if (sIntegrin.get(i).BoundLigandYN==0) {
				for(int j=0;j<nOfLigands;j++){	
				double dist=distance(sIntegrin.get(i).Bead_Position,sLigand.get(j).Bead_PositionL);
				if (sIntegrin.get(i).BoundLigandYN==0  && sLigand.get(j).BoundYN==0 && dist<r_min){						
					r_min=dist;
					lig=j;
				}
			}
			if (r_min<=0.020001 && sIntegrin.get(i).kon*sIntegrin.get(i).dt >Math.random()){
				sLigand.get(lig).integrinBoundNum=i;
				sLigand.get(lig).BoundYN=1;
				sIntegrin.get(i).BoundLigandYN=1;
				sIntegrin.get(i).LigandNum=lig;
				double force=k_it*( r_min-(0.0001));
				sLigand.get(lig).FF_interact=force;
				sLigand.get(lig).Bead_PositionL=move(sLigand.get(lig).Bead_PositionL, sIntegrin.get(i).Bead_Position,  force, r_min, sLigand.get(lig).dt, sLigand.get(lig).zetaL);	
				sIntegrin.get(i).Bead_Position=move(sIntegrin.get(i).Bead_Position, sLigand.get(lig).Bead_PositionL, force, r_min, sIntegrin.get(i).dt, sIntegrin.get(i).zetaI);
				
			}
			}
		}
	}
	
	
protected void CatchBond(){
	
	// catch bond parameters
	double base=2;
	double 	A=-0.1;
	double	B=0.3;
     double C=0.00004;
	


	
	for(int i=0;i<nOfLigands;i++){
		if (sLigand.get(i).BoundYN==1){
		 double P;
		 double force=sLigand.get(i).FF_interact;
		
		 sIntegrin.get(sLigand.get(i).integrinBoundNum).kub=base*Math.exp(A*force)+C*Math.exp(B*force);	
		  sIntegrin.get(sLigand.get(i).integrinBoundNum).FinalLife=0;
				 	 P=(sIntegrin.get(sLigand.get(i).integrinBoundNum).kub*sLigand.get(i).dt);
		 	sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT=sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT+1;
		 	if (Math.random()<P )
		 	{	
		 		sLigand.get(i).BoundYN=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).finalTension=sLigand.get(i).FF_interact;
		 		sLigand.get(i).FF_interact=0;
		 		//////CHANGE
		 		//sIntegrin.get(sLigand.get(i).integrinBoundNum).kon=0.1;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).BoundLigandYN=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).FinalLife=sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).contaBoundLifetimeTot=(sIntegrin.get(sLigand.get(i).integrinBoundNum).contaBoundLifetimeTot+sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT*sIntegrin.get(sLigand.get(i).integrinBoundNum).dt);
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).contaAT=0;
		 		sIntegrin.get(sLigand.get(i).integrinBoundNum).kub=0;
		 		sLigand.get(i).integrinBoundNum=-1;
		 		}		
		 	 }
	}
}


	private double distance(double[] vecA, double[] vecB) {
		return Math.sqrt((vecA[0]-vecB[0])*(vecA[0]-vecB[0])+(vecA[1]-vecB[1])*(vecA[1]-vecB[1]));
	}
	
	private double[] move(double[] beadA, double[] beadB, double force, double r, double dt, double zeta) {
		double temp[] = new double[3];
		temp[0] = beadA[0]-(force)*(beadA[0]-beadB[0])/r*dt/zeta;
		temp[1] = beadA[1]-force*(beadA[1]-beadB[1])/r*dt/zeta;
		return temp;
	}
	
}