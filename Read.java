package plos_comp_bio;


import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Random;

public class Read {
	public double[] start(){
		double [] vector =new double[10];;
	File myObj = new File("config.yaml");
    if (myObj.exists()) {
    	//System.out.println(myObj.length());
    	 BufferedReader fr;
    	 
    	 
		try {
		
			for (int ii=0; ii<10; ii++)
					{
				String word = null;
				switch (ii) {
				case 0:
					word="timestep";
					break;
				case 1:
					word="nIntegrins";
					break;
				case 2:
					word="nLigands";
					break;
				case 3:
					word="total-time";
					break;
				case 4:
					word="friction";
					break;
				case 5:
					word="k_on";
					break;
				case 6:
					word="spring";
					break;
				case 7:
					word="width";
					break;
				case 8:
					word="bundling";
					break;
				case 9:
					word="force";
					break;
				}
				System.out.println(ii+   word);
			fr = new BufferedReader(new FileReader("config.yaml"));
			 String contentLine = fr.readLine();
			  while (contentLine != null) 
			   {
			      contentLine = fr.readLine();
			      if (contentLine != null && contentLine.contains(word)) 
			      {	
			    	  int indx=contentLine.indexOf(":");
			    	  String sub=contentLine.substring(indx+1,contentLine.length());
			    	  vector[ii]=Double.parseDouble(sub.strip());
			    	 
			      }	     
			   } 
			  
					}
	  
		}catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    	
catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
       
      
    } else {
      System.out.println("The file does not exist.");
      System.out.println("Absolute path: " + myObj.getAbsolutePath());
    }
	return vector;
	}
}
