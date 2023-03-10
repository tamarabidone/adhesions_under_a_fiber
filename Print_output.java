package plos_comp_bio;
import java.util.ArrayList;
import java.io.FileWriter;
import java.io.IOException;


public class Print_output{
	
	ArrayList<Integrin_2D> sIntegrin = new ArrayList<Integrin_2D>();
	ArrayList<Ligand_2D> sLigand = new ArrayList<Ligand_2D>();
	int nOfIntegrins,nOfLigands;
	double timestep;

	public Print_output(ArrayList<Integrin_2D> sIntegrin,  ArrayList<Ligand_2D> sLigand, int nOfIntegrins, int nOfLigands, double timestep){
		this.sIntegrin= sIntegrin;
		this.sLigand= sLigand;
		this.nOfIntegrins = nOfIntegrins;	
		this.nOfLigands = nOfLigands;	
		this.timestep=timestep;
	}

		/////This is for output files////
	
	public void print_all() {
		int numLigandBoundIntegrin = 0; // this variable stores the total ligand bound integrins at each timestep
		for (int i = 0; i < nOfIntegrins; i++) {
			if (sIntegrin.get(i).BoundLigandYN == 1) {
				numLigandBoundIntegrin = numLigandBoundIntegrin + 1;
			}
			
		}
		

		String l2 =  "T1_"+timestep+"_" +sIntegrin.get(0).Pbundling+"__"+ sIntegrin.get(0).MyoForce + "_Integrin_"
				+ sIntegrin.get(0).V * 1E3 + "_v_" + sLigand.get(0).k;
		String l3 =  "T1_finalTension_"+timestep+"_"+"_"+ "_"  +sIntegrin.get(0).Pbundling+"__"+  sIntegrin.get(0).MyoForce+ "_IntegrinFinals_"
				+ sIntegrin.get(0).V * 1E3 + "_v_" + sLigand.get(0).k;
		String ld ="T1__"+timestep+"_"+"_Bound_" +  "_"+ "_"  +sIntegrin.get(0).Pbundling+"__"+  sIntegrin.get(0).MyoForce + "_L_"
				+ sIntegrin.get(0).V * 1E3 + "_v_" + sLigand.get(0).k ;
		String l4 = "T1_"+timestep+"_kon_"+"General_"  +"_"+"_" + "_"+sIntegrin.get(0).Pbundling+"__"+  sIntegrin.get(0).MyoForce + "_L_"
				+ sIntegrin.get(0).V * 1E3 + "_v_" + sLigand.get(0).k;

//// Saving output files every second
		if (sIntegrin.get(0).time % 1 < timestep) {
			{
				//System.out.println(sLigand.get(0).k);
				try {
					FileWriter fw2 = new FileWriter(l2, true);
					FileWriter fw3 = new FileWriter(l3, true);
					FileWriter fw4 = new FileWriter(l4, true);
					fw4.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + "   "
							+ String.format("%d", (numLigandBoundIntegrin)) + " \n");// appends the string to the file
					fw4.close();
					for (int i = 0; i < nOfIntegrins; i++) {
						fw2.write(String.format("%.3f", sIntegrin.get(0).time) + "     " + i + "    "
								+ String.format("%d", sIntegrin.get(i).BoundLigandYN) + "    "
								+ String.format("%.3f", sIntegrin.get(i).Bead_Position[0]) + "     "
								+ String.format("%.3f", sIntegrin.get(i).Bead_Position[1]) + "     "
								+ String.format("%d", sIntegrin.get(i).contaAT) + "     "
								+ String.format("%.6f", sIntegrin.get(i).finalTension) + "   "
								+ String.format("%.6f", sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt) + "     "
								+ String.format("%d", sIntegrin.get(i).TimesUnbound) + "     "
								+ String.format("%.4f", sIntegrin.get(i).contaBoundLifetimeTot) + " \n");// appends the
																											// string to
																											// the file

						if (sIntegrin.get(i).FinalLife > 0)
							fw3.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + i + "   "
									+ String.format("%.6f", sIntegrin.get(i).finalTension) + "  "
									+ String.format("%.6f", sIntegrin.get(i).FinalLife * sIntegrin.get(i).dt) + "     "
									+ " \n");// appends the string to the file
					}
					fw2.close();
					fw3.close();
					FileWriter fw1 = new FileWriter(ld, true);

					for (int i = 0; i < 441; i++) {
						if (sLigand.get(i).BoundYN == 1)
							fw1.write(String.format("%.3f", sIntegrin.get(0).time) + "    " + String.format("%d", i)
									+ "    " + String.format("%d", sLigand.get(i).BoundYN) + "       "
									+ String.format("%.3f", sLigand.get(i).Bead_PositionL[0]) + "    "
									+ String.format("%.3f", sLigand.get(i).Bead_PositionL[1]) + "    "
									+ String.format("%d", sLigand.get(i).integrinBoundNum) + "    "
									+ String.format("%3.4f", sLigand.get(i).FF_interact) + "    "
								+ " \n");// appends the string to the file
					}
					fw1.close();
				} catch (IOException ioe) {
					System.err.println("IOException: " + ioe.getMessage());
				}

			}
		}
	}
	
}