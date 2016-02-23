package prop;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import ilog.concert.*;
import ilog.cplex.*;

public class LGLB {
	static void usage() {
		System.out.println("usage: java -jar thisfile.jar path");
	}
	
	public static void main(String[] args) {
		if(args.length!=1) {
			System.out.println("please enter the path of property file");
			usage();
			return;
		}
		
		PropertyHelper.loadProperties(args[0]);
		
		//get inputs from .properties file
		int k = PropertyHelper.getIntValue("k");
		int n = PropertyHelper.getIntValue("n");
		int l = PropertyHelper.getIntValue("l");
		int V = PropertyHelper.getIntValue("V");
		double T = PropertyHelper.getDoubleValue("T");
		double Cvm = PropertyHelper.getDoubleValue("Cvm");
		double Evm = PropertyHelper.getDoubleValue("Evm");
		double Epr = PropertyHelper.getDoubleValue("Epr");
		double Eps = PropertyHelper.getDoubleValue("Eps");
		double Cp = PropertyHelper.getDoubleValue("Cp");
		int[] o = PropertyHelper.getIntArray("o", k);
		int[] C = PropertyHelper.getIntArray("C", l);
		int[][] SC = PropertyHelper.getInt2DArray("S", k, l);
		double[] Emin = PropertyHelper.getDoubleArray("Emin", k);
		double[] Emax = PropertyHelper.getDoubleArray("Emax", k);
		double[] G = PropertyHelper.getDoubleArray("G", k);
		double[] L = PropertyHelper.getDoubleArray("L", n);
		double[] GN = PropertyHelper.getDoubleArray("GN", n);
		double[][] e = PropertyHelper.getDouble2DArray("e", k);
		
		//calculate M
		int [] S = new int[k];
		for(int i=0; i<k; i++) {
			int sumS = 0;
			for(int j=0; j<l; j++) {
			sumS = SC[i][j] + sumS;
			}
			S[i] = sumS;
		}

		int Cmax = C[0];  
		for (int c = 0; c < l; c++) {  
			if (C[c] > Cmax) { 
				Cmax = C[c];  
			}
		}  
	
		int tp = 0;
		double max = L[0];  
		for (int i = 1; i < n; i++) {  
			if (L[i] > max) { 
				max = L[i]; 
				tp = i;
			}
		}  
		 
		double[][] Lp = new double[k][tp+Cmax];
		for (int i=0; i<k; i++) {
			for (int t=0; t<tp; t++) {
				Lp[i][t]=0;
			}
		}
		for (int i=0; i<k; i++) {
			for (int t=tp; t<tp+Cmax; t++) {
				double sumLp = 0;
				for (int a=0; a <= Cmax-1; a++) {
					for (int j=0; j<l; j++) {
						if((C[j]>a)&&(t-a>=0))
							sumLp = sumLp + java.lang.Math.ceil(SC[i][j]*L[t-a]);
					}
				}
				Lp[i][t]=sumLp;
			}
		}
		
		double [] Lpmax = new double[k];		
		for (int i = 0; i < k; i++) {  
			Lpmax[i]=Lp[i][tp];
			for(int t=tp; t<tp+Cmax; t++) {
				if(Lp[i][t] > Lpmax[i])
					Lpmax[i] = Lp[i][t];
			}
		}  
		
		int[] M = new int[k];
		//if there is the key "M" in the properties file, then read it, else use the calculated M
		if(PropertyHelper.HasProperty("M")) {
			M = PropertyHelper.getIntArray("M", k);
			for(int i = 0; i < M.length; i++) {
				if (M[i] <= 0) {
					System.err.println("M should be positive!");
					System.exit(1);
				}
			}
		}
		else {
			for (int i = 0; i < k; i++) {
				M[i] = Math.max((int)java.lang.Math.ceil(((Lpmax[i]) / V) / 0.8), 1);
			}
		}
		
		int[] t = new int[n+1];
		for (int i=0; i<=n; i++) {
			t[i] = i;
		}
		
		//calculate SCt 
		int[][][] SCt = new int[k][l][n];
		for (int i=0; i<k; i++) {
			for (int j=0; j<l; j++) {
				for (int tt=0; tt<n; tt++) {
					int tl=((t[tt]-o[i])<0)?(t[tt]-o[i]+n):(t[tt]-o[i]);
					SCt[i][j][tt] = (int)java.lang.Math.ceil(SC[i][j]*L[tl]);
				}
			}
		}
		//St should be calculated as the sum of SCt
		int[][] St = new int[k][n];
		for (int i=0; i<k; i++) {
			for (int tt=0; tt<n; tt++) {
				int sumSCt = 0;
				for (int j=0; j<l; j++) {
					sumSCt = sumSCt + SCt[i][j][tt];
				}
				St[i][tt] = sumSCt;
			}
		}
		
		//calculate Gt
		double[][] Gt = new double[k][n];
		for (int i=0; i<k; i++) {
			for (int j=0; j<n; j++) {
				int tl=((t[j]-o[i])<0)?(t[j]-o[i]+n):(t[j]-o[i]);
				Gt[i][j] = G[i]*GN[tl];
			}
		}
		//fin of reading inputs and necessary calculation before modeling
		
		//modeling
		try {
			//define new model
			IloCplex cplex = new IloCplex();
			cplex.setParam(IloCplex.IntParam.NodeFileInd, 3);
			cplex.setParam(IloCplex.DoubleParam.TiLim, T);
			//variables 
			IloIntVar[][] x = new IloIntVar[k][n+1];
			IloIntVar[][] xp = new IloIntVar[k][n+1];	
			IloIntVar[][] y = new IloIntVar[k][n+1];		
			IloIntVar[][] yp = new IloIntVar[k][n+1];	
			IloIntVar[][] R = new IloIntVar[k][n+1];
			IloIntVar[][][] RC = new IloIntVar[k][l][n+1];
			IloIntVar[][][] lambda = new IloIntVar[k][k][n+1];	
			IloIntVar[][][][] lambdaC = new IloIntVar[k][k][l][n+1];	
			IloNumVar[][] F = new IloNumVar[k][n+1];
			for (int i=0; i<k; i++) {
				x[i] = cplex.intVarArray(n+1, 0, M[i]*V);
				xp[i] = cplex.intVarArray(n+1, 0, M[i]*V);
				y[i] = cplex.intVarArray(n+1, 0, M[i]);
				yp[i] = cplex.intVarArray(n+1, 0, M[i]);
				R[i] = cplex.intVarArray(n+1, 0, M[i]*V);
				F[i] = cplex.numVarArray(n+1, 0, Double.MAX_VALUE);
				for (int j=0; j<l; j++){
					RC[i][j] = cplex.intVarArray(n+1, 0, M[i]*V);
				}
				
				for (int j=0; j<k; j++){
					lambda[i][j] = cplex.intVarArray(n+1, 0, Integer.MAX_VALUE);
				}

				for (int j=0; j<k; j++){
					for (int c=0; c<l; c++) {
						lambdaC[i][j][c] = cplex.intVarArray(n+1, 0, Integer.MAX_VALUE);
					}
				}
			}
			
			//function objective
			IloLinearNumExpr objective = cplex.linearNumExpr();
			
			for (int i=0; i<k; i++) {
				for (int j=1; j<=n; j++) {
					objective.addTerm(1, F[i][j]);
				}
			}
			
			IloLinearNumExpr routage = cplex.linearNumExpr();		
			for (int tt=1; tt<=n; tt++) { 			
				for (int i=0; i<k; i++) {
					for (int j=0; j<k; j++) {
						routage.addTerm(e[i][j], lambda[i][j][tt]);
					}
				}
			}
			
			objective.add(routage);
			
			//define objective
			cplex.addMinimize(objective);
			//constraints	
			//when t=0, all the physical machines are in running mode, and each of them instances V VMs
			for (int i=0; i<k; i++) {		
				cplex.addEq(x[i][0], M[i]*V);
				cplex.addEq(xp[i][0], 0);
				cplex.addEq(y[i][0], M[i]);
				cplex.addEq(yp[i][0], 0);
				cplex.addEq(R[i][0], 0);
				for (int j=0; j<l; j++) {
					cplex.addEq(RC[i][j][0], 0);
				}
				cplex.addEq(F[i][0], 0);//t=0 F=0
			}
			//no routing when t=0, no routing when j=i
			for (int i=0; i<k; i++) {
				for (int j=0; j<k; j++) {
					cplex.addEq(lambda[i][j][0], 0);
					for (int c=0; c<l; c++) {
						cplex.addEq(lambdaC[i][j][c][0], 0);
						for (int tt=1; tt<=n; tt++) {
							if(j==i)
								cplex.addEq(lambdaC[i][j][c][tt], 0);	
						}		
					}
					for (int tt=1; tt<=n; tt++) {
						if(j==i)
							cplex.addEq(lambda[i][j][tt], 0);	
					}
				}
			}
			
			//InDC <= Fi(t)
			for (int i=0; i<k; i++) {		
				for (int tt=1; tt<=n; tt++) {
					IloLinearNumExpr InDC = cplex.linearNumExpr();	
					InDC.addTerm(Evm, x[i][tt]);
					InDC.addTerm(Cvm, xp[i][tt]);
					InDC.addTerm(Epr-Eps, y[i][tt]);
					InDC.addTerm(Cp, yp[i][tt]);
					InDC.addTerm((Emax[i]-Emin[i])/(M[i]*V), R[i][tt]);
					InDC.setConstant(Emin[i]+M[i]*Eps-Gt[i][tt-1]);
					cplex.addLe(InDC, F[i][tt]);	
				}
			}
			
			//For each i, RCi(t)=SCi(t)+nb_in-nb_out, nb_out<=SCi(t) => RCi(t)>=0	
			for (int i=0; i<k; i++) {	
				for (int c=0; c<l; c++) {
					for (int tt=1; tt<=n; tt++) {
						IloLinearNumExpr nb_out = cplex.linearNumExpr();	
						IloLinearNumExpr nb_in = cplex.linearNumExpr();		
						for (int j=0; j<k; j++) {
							nb_out.addTerm(1, lambdaC[i][j][c][tt]);
							nb_in.addTerm(1, lambdaC[j][i][c][tt]);
						}		
						//Ric(t)=Sic(t)+nb_in-nb_out
						cplex.addEq(RC[i][c][tt], cplex.sum(cplex.sum(nb_in, SCt[i][c][tt-1]), cplex.prod(-1, nb_out)));
						//For each i and each t, nb_out<=Si(t)
						cplex.addLe(nb_out, SCt[i][c][tt-1]);
					}
				}
			}
			
			//calculate Ri(t), xi(t)>=Ri(t)
			for (int i=0; i<k; i++) {
				for (int tt=1; tt<=n; tt++) {
					IloLinearNumExpr sumR = cplex.linearNumExpr();	
					for (int a=0; a <= Cmax-1; a++) {
						for (int c=0; c<l; c++) {
							if((C[c]>a)&&(tt-a>=0))
								sumR.addTerm(1, RC[i][c][tt-a]);
						}
					}
					cplex.addEq(R[i][tt], sumR);	
					cplex.addGe(x[i][tt], R[i][tt]);
				}
			}
			
			//sum of lambdaC = lambda
			for (int i=0; i<k; i++) {	
				for (int j=0; j<k; j++) {
					for (int tt=1; tt<=n; tt++) {
						IloLinearNumExpr sumLC = cplex.linearNumExpr();	
						for (int c=0; c<l; c++) {
							sumLC.addTerm(1, lambdaC[i][j][c][tt]);
						}
						cplex.addEq(lambda[i][j][tt], sumLC);	
					}
				}
			}
			
			
			//yi(t)<=M(i)
			for (int i=0; i<k; i++) {		
				for (int j=1; j<=n; j++) {
					cplex.addLe(y[i][t[j]], M[i]);
				}
			}
			//x(t)<=V*y(t)
			for (int i=0; i<k; i++) {		
				for (int j=1; j<=n; j++) {
					cplex.addLe(x[i][t[j]], cplex.prod(V, y[i][t[j]]));
				}
			}
			//xp(t)>=x(t)-x(t-1)
			for (int i=0; i<k; i++) {			
				for (int j=1; j<=n; j++) {
					cplex.addGe(xp[i][t[j]], cplex.sum(x[i][t[j]], cplex.prod(-1, x[i][t[j-1]])));
				}
			}
			//yp(t)>=y(t)-y(t-1)
			for (int i=0; i<k; i++) {			
				for (int j=1; j<=n; j++) {
					cplex.addGe(yp[i][t[j]], cplex.sum(y[i][t[j]], cplex.prod(-1, y[i][t[j-1]])));
				}
			}
			//fin of modeling
			
	        // write model to file
			cplex.exportModel("lglb.lp");
		
			// solve the model and display the solution if one was found
	        if ( cplex.solve() ) {
	    		System.out.println("***********************************************Input and Result*************************************************");
	    		System.out.println("****************************************************Input*******************************************************");
	    		System.out.println("Number of DC -- k : " + k);
	    		System.out.println("Number of time interval -- n : " + n);
	    		System.out.println("Time limit -- T : " + T);
	    		System.out.println("Maximum number of VMs that can be instanced on one physical machine -- V  : "+ V);	
	    		System.out.println("Energy consumption of creation of a VM -- Cvm: "+ Cvm);
	    		System.out.println("Energy consumption of execution of a VM during a timeslot -- Evm : "+ Evm);
	    		System.out.println("Energy consumption of execution of a physical machine in running mode during a timeslot -- Epr : "+ Epr);
	    		System.out.println("Energy consumption of execution of a physical machine in sleep mode during a timeslot -- Eps : "+ Eps);
	    		System.out.println("Energy consumption of the action leave sleep mode for a physical machine -- Cp : "+ Cp);
	    		System.out.print("\n");
	    		for (int i=0; i<k; i++) {
	    			System.out.println("Number of the physical machines in DC " + i + " -- M" + i + " : "+ M[i]);
	    		}
	    		System.out.print("\n");
	    		for (int i=0; i<k; i++) {
	    			System.out.println("Time offset of DC " + i + " -- O" + i + " : "+ o[i]);
	    		}
	    		System.out.print("\n");
	    		System.out.print("\n");
	    		System.out.println("Class of request, in terms of the durations of video : C");
	    		for (int c=0; c<l; c++) {
	    			System.out.print("\t C["+ c + "] : " + C[c]);
	    		}
	    		System.out.print("\n");
	    		System.out.println("Daily average number of requests c sent by users in DC i : SC(i)(c)");
	    		for(int i=0; i<k; i++) {
	    			for(int c=0; c<l; c++) {
	    				System.out.print("\t"+"SC["+i+"]["+c+"]: "+SC[i][c]);
	    			}
		    		System.out.print("\n");	
	    		}
	    		System.out.print("\n");
	    		for (int i=0; i<k; i++) {
	    			System.out.println("Daily average number of requests sent by users in DC " + i + " -- S" + i + " : "+ S[i]);
	    		}
	    		System.out.print("\n");
	    		for (int i=0; i<k; i++) {
	    			System.out.println("The parameter Emin of DC" + i + " -- Emin" + i + " : "+ Emin[i]);
	    		}
	    		System.out.print("\n");
	    		for (int i=0; i<k; i++) {
	    			System.out.println("The parameter Emax of DC" + i + " -- Emax" + i + " : "+ Emax[i]);
	    		}	
	    		System.out.print("\n");
	    		for (int i=0; i<k; i++) {
	    			System.out.println("Parameter to calculate the renewable energy at time t in DC " + i + " -- G" + i + " : "+ G[i]);
	    		}
	    		System.out.print("\n");
	    		System.out.println("Load function which presents the hourly load for one day : L(t)");
	    		for(int i=0; i<n; i++) {
	    			if((i%6==0)&&(i!=0)) {
	    				System.out.print("\n");
	    			}
	    			System.out.print("\t "+ L[i]);
	    		}
	    		System.out.println("\n");
	    		System.out.println("Renewable energy function which presents the hourly renewable energy availble for one day : L(t)");
	    		for(int i=0; i<n; i++) {
	    			if((i%6==0)&&(i!=0)) {
	    				System.out.print("\n");
	    			}
	    			System.out.print("\t "+ GN[i]);
	    		}
	    		System.out.println("\n");
	    		System.out.println("Energy consumption of routing of a request from DC i to DC j : e(i)(j)");
	    		for(int i=0; i<k; i++) {
	    			for(int j=0; j<k; j++) {
	    				System.out.print("\t"+"e["+i+"]["+j+"]: "+e[i][j]);
	    			}
		    		System.out.print("\n");	
	    		}
	    		System.out.println("****************************************************Result******************************************************");
	    		NumberFormat nf = new DecimalFormat("0.0");
	    		//0--12
	    		for (int i=0; i<k; i++) {		
	    			System.out.println("****************************************************DC("+ i +")******************************************************");
		    		System.out.print("t");
		    		for (int j=0; j<=n/2; j++) {
		    			System.out.print("\t"+j);
		    		}
		    		System.out.println("\n");
		    	
		    		System.out.print("St");
		    		System.out.print("\t"+"0");
		    		for (int j=0; j<n/2; j++) {
		    			System.out.print("\t"+St[i][j]);
		    		}
		    		System.out.print("\n");
		    		
		    		System.out.print("SCt");	  
		    		System.out.print("\n");
		    		//jj iterator of class, j iterator of t
		    		for (int jj=0; jj<l; jj++) {
		    			System.out.print("c="+C[jj]);
			    		System.out.print("\t"+"0");
		    			for (int j=0; j<n/2; j++) {
			    			System.out.print("\t"+SCt[i][jj][j]);
		    			}	
		    			System.out.print("\n");
		    		}
		    		System.out.print("\n");	
		    		
		    		System.out.print("R");
		    		for (int j=0; j<=(n+1)/2; j++) {
		            	cplex.output().print("\t"+Math.round(cplex.getValue(R[i][j])));
		    		}
		    		System.out.print("\n");
		    		
		    		System.out.print("RC");	  
		    		System.out.print("\n");
		    		//jj iterator of class, j iterator of t
		    		for (int jj=0; jj<l; jj++) {
		    			System.out.print(C[jj]);
		    			for (int j=0; j<=(n+1)/2; j++) {
		    				cplex.output().print("\t"+Math.round(cplex.getValue(RC[i][jj][j])));
		    			}	
		    			System.out.print("\n");
		    		}
		    		System.out.print("\n");	
	
		    		System.out.print("x");
		    		for (int j=0; j<=(n+1)/2; j++) {
		            	cplex.output().print("\t"+Math.round(cplex.getValue(x[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("xp");
		    		for (int j=0; j<=(n+1)/2; j++) {
						cplex.output().print("\t"+Math.round(cplex.getValue(xp[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("y");
		    		for (int j=0; j<=(n+1)/2; j++) {
		    			cplex.output().print("\t"+Math.round(cplex.getValue(y[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("yp");
		    		for (int j=0; j<=(n+1)/2; j++) {
		            	cplex.output().print("\t"+Math.round(cplex.getValue(yp[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		int[][] z = new int[k][n+1];
		    		System.out.print("z");
		    		for (int j=0; j<=(n+1)/2; j++) {
		    			z[i][j] = (int)(M[i]-Math.round(cplex.getValue(y[i][j])));
		            	cplex.output().print("\t"+z[i][j]);
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("lambda ");	
		    		System.out.print("\n");
		    		//jj iterator of DC, j iterator of t
		    		for (int jj=0; jj<k; jj++) {
		    			System.out.print(i + " to "+ jj);
		    			for (int j=0; j<=(n+1)/2; j++) {
		    				cplex.output().print("\t"+Math.round(cplex.getValue(lambda[i][jj][j])));
		    			}	
		    			System.out.print("\n");	
		    		}
		    		
		    		System.out.print("lambdaC ");	
		    		System.out.print("\n");
		    		for (int c=0; c<l; c++) {
		    			System.out.println("c="+C[c]);
		    			//jj iterator of DC, j iterator of t
		    			for (int jj=0; jj<k; jj++) {
		    				System.out.print(i + " to "+ jj);
		    				for (int j=0; j<=(n+1)/2; j++) {
		    					cplex.output().print("\t"+Math.round(cplex.getValue(lambdaC[i][jj][c][j])));
		    				}	
		    			System.out.print("\n");	
		    			}
		    		}
		    		System.out.print("\n");
		    		
		    		System.out.print("Gt");
		    		System.out.print("\t"+"0");
		    		for (int j=0; j<n/2; j++) {
		    			System.out.print("\t"+nf.format(Gt[i][j]));
		    		}
		    		System.out.println("\n");
		    		
		    		System.out.print("F");
		    		for (int j=0; j<=(n+1)/2; j++) {
		            	cplex.output().print("\t"+nf.format(cplex.getValue(F[i][j])));
					}	
		    		System.out.println("\n");
		    		//13--24
		    		System.out.print("t");
		    		for (int j=n/2+1; j<=n; j++) {
		    			System.out.print("\t"+j);
		    		}
		    		System.out.println("\n");
		    	
		    		System.out.print("St");
		    		for (int j=n/2; j<n; j++) {
		    			System.out.print("\t"+St[i][j]);
		    		}
		    		System.out.print("\n");
	
		    		System.out.print("SCt");	  
		    		System.out.print("\n");
		    		//jj iterator of class, j iterator of t
		    		for (int jj=0; jj<l; jj++) {
		    			System.out.print("c="+C[jj]);
		    			for (int j=n/2; j<n; j++) {
			    			System.out.print("\t"+SCt[i][jj][j]);
		    			}	
		    			System.out.print("\n");
		    		}
		    		System.out.print("\n");		    
		    		
		    		System.out.print("R");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
		            	cplex.output().print("\t"+Math.round(cplex.getValue(R[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("RC");	  
		    		System.out.print("\n");
		    		//jj iterator of class, j iterator of t
		    		for (int jj=0; jj<l; jj++) {
		    			System.out.print(C[jj]);
		    			for (int j=(n+1)/2+1; j<n+1; j++) {
		    				cplex.output().print("\t"+Math.round(cplex.getValue(RC[i][jj][j])));
		    			}	
		    			System.out.print("\n");
		    		}
		    		System.out.print("\n");		    
		    		
		    		System.out.print("x");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
		            	cplex.output().print("\t"+Math.round(cplex.getValue(x[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("xp");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
						cplex.output().print("\t"+Math.round(cplex.getValue(xp[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("y");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
		            	
		    			cplex.output().print("\t"+Math.round(cplex.getValue(y[i][j])));
					}	
		    		System.out.print("\n");
	
		    		System.out.print("yp");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
		            	cplex.output().print("\t"+Math.round(cplex.getValue(yp[i][j])));
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("z");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
		    			z[i][j] = (int)(M[i]-Math.round(cplex.getValue(y[i][j])));
		            	cplex.output().print("\t"+z[i][j]);
					}	
		    		System.out.print("\n");
		    		
		    		System.out.print("lambda");	  
		    		System.out.print("\n");
		    		//jj iterator of DC, j iterator of t
		    		for (int jj=0; jj<k; jj++) {
		    			System.out.print(i + " to "+ jj);
		    			for (int j=(n+1)/2+1; j<n+1; j++) {
		    				cplex.output().print("\t"+Math.round(cplex.getValue(lambda[i][jj][j])));
		    			}	
		    			System.out.print("\n");
		    		}
		    		
		    		System.out.print("lambdaC ");	
		    		System.out.print("\n");
		    		for (int c=0; c<l; c++) {
		    			System.out.println("c="+C[c]);
		    			//jj iterator of DC, j iterator of t
		    			for (int jj=0; jj<k; jj++) {
		    				System.out.print(i + " to "+ jj);
		    				for (int j=(n+1)/2+1; j<n+1; j++) {
		    					cplex.output().print("\t"+Math.round(cplex.getValue(lambdaC[i][jj][c][j])));
		    				}	
		    			System.out.print("\n");	
		    			}
		    		}
		    		System.out.print("\n");
		    		
		    		System.out.print("Gt");
		    		for (int j=n/2; j<n; j++) {
		    			System.out.print("\t"+nf.format(Gt[i][j]));
		    		}
		    		System.out.println("\n");
		    		
		    		System.out.print("F");
		    		for (int j=(n+1)/2+1; j<n+1; j++) {
		            	cplex.output().print("\t"+nf.format(cplex.getValue(F[i][j])));
					}	
		    		System.out.println("\n");		
	    		}
	    		
		    	System.out.println("St : the number of requests sent by users at time t");
		    	System.out.println("SCt : the number of requests c sent by users at time t");
		    	System.out.println("R : the number of requests processed locally at time t");
		    	System.out.println("RC : the number of requests c received at at time t to be processed locally");
		    	System.out.println("x : the number of VMs at time t");
		    	System.out.println("xp : the number of VMs created at time t");
		    	System.out.println("y : the number of physical machines in running mode at time t");
		    	System.out.println("xp : the number of physical machines leave sleep mode at time t");
		    	System.out.println("z : the number of physical machines in sleep mode at time t");
		    	System.out.println("lambda : the number of requests sent from DCi to DCj at time t");
		    	System.out.println("lambdaC : the number of requests c sent from DCi to DCj at time t");
		    	System.out.println("Gt : the renewable energy available at time t");
		    	System.out.println("\n");
		    		
		    	cplex.output().println("Solution status = " + cplex.getStatus());
		        cplex.output().println("Objective value = " + cplex.getObjValue());
	        
		        //sumX: total energy consumption of execution of VMs
	    		double sumX = 0;	
	    		for (int i=0; i<k; i++) {
	    			for (int j=1; j<=n; j++) {
		    		sumX=sumX+Evm*cplex.getValue(x[i][j]);
	    			}
	    		}	
	    		
	    		//sumXP: total energy consumption of creation of VMs
	    		double sumXP = 0;	
	    		for (int i=0; i<k; i++) {
	    			for (int j=1; j<=n; j++) {
		    		sumXP=sumXP+Cvm*cplex.getValue(xp[i][j]);
	    			}
	    		}
	    		
	    		//sumY: total energy consumption of execution of physical machines that on running mode 
	    		double sumY = 0;	
	    		for (int i=0; i<k; i++) {
	    			for (int j=1; j<=n; j++) {
		    		sumY=sumY+Epr*cplex.getValue(y[i][j]);
	    			}
	    		}	
    		
	    		//sumYP: total energy consumption of the action of leaving sleep mode on physical machines
	    		double sumYP = 0;	
	    		for (int i=0; i<k; i++) {
	    			for (int j=1; j<=n; j++) {
			    		sumYP=sumYP+Cp*cplex.getValue(yp[i][j]);
		    		}
		    	}
	        
	    		//sum Z: total energy consumption of execution of physical machines that on sleep mode 
    			double sumZ = 0;	
    			for (int i=0; i<k; i++) {
    				for (int j=1; j<=n; j++) {
	    			sumZ=sumZ+Eps*(M[i]-cplex.getValue(y[i][j]));
    				}
    			}	
    			
    			//sumE: total energy consumption of non-IT infrastructure in DCs
    			double E = 0;
    			for (int i=0; i<k; i++) {
    				for (int j=1; j<=n; j++) {
    					E = E + Emin[i]; 
    				}
    			}
    			double sumR = 0;	
    			for (int i=0; i<k; i++) {
    				for (int j=1; j<=n; j++) {
	    			sumR=sumR+((Emax[i]-Emin[i])/(M[i]*V))*cplex.getValue(R[i][j]);
    				}
    			}	
    			double sumE = sumR + E;
    			
    			//sumLambda: total energy consumption of migrations
    			double sumLambda = 0;
    			for (int i=0; i<k; i++) {
    				for (int tt=1; tt<=n; tt++) {
    					for (int j=0; j<k; j++) {
    						sumLambda = sumLambda + e[i][j]*cplex.getValue(lambda[i][j][tt]);
    					}
    				}
    			}
    			
    			//sumF
    			double sumF = 0;	
    			for (int i=0; i<k; i++) {
    				for (int j=1; j<=n; j++) {
	    			sumF=sumF+cplex.getValue(F[i][j]);
    				}
    			}
    			
    			System.out.println("****************************************************LGLB********************************************************");
    			System.out.print("Part");
	    		System.out.print("\t"+"SumX");
	    		System.out.print("\t"+"SumXP");
	    		System.out.print("\t"+"SumY");
	    		System.out.print("\t"+"SumYP");
	    		System.out.print("\t"+"SumZ");
	    		System.out.print("\t"+"sumE");
	    		System.out.print("\t"+"SumLambda");
	    		System.out.println("\n");
	    		System.out.print("Value");
	    		System.out.print("\t"+nf.format(sumX));
	    		System.out.print("\t"+nf.format(sumXP));
	    		System.out.print("\t"+nf.format(sumY));
	    		System.out.print("\t"+nf.format(sumYP));
	    		System.out.print("\t"+nf.format(sumZ));    		
	    		System.out.print("\t"+nf.format(sumE));
	    		System.out.print("\t"+nf.format(sumLambda));    	
	    		System.out.println("\n");
	    		double sum=sumX+sumXP+sumY+sumYP+sumZ+sumE+sumLambda;	
	    		NumberFormat nt = NumberFormat.getPercentInstance();
	    		nt.setMinimumFractionDigits(2);
	    		System.out.print("PCT");
	    		System.out.print("\t"+nt.format(sumX/sum));
	    		System.out.print("\t"+nt.format(sumXP/sum));
	    		System.out.print("\t"+nt.format(sumY/sum));
	    		System.out.print("\t"+nt.format(sumYP/sum));
	    		System.out.print("\t"+nt.format(sumZ/sum));   		
	    		System.out.print("\t"+nt.format(sumE/sum));
	    		System.out.print("\t"+nt.format(sumLambda/sum));
	    		System.out.println("\n");   
	    		System.out.print("sumX: total energy consumption of execution of VMs");
	    		System.out.print("sumXP: total energy consumption of creation of VMs");
	    		System.out.print("sumY: total energy consumption of execution of physical machines that on running mode ");
	    		System.out.print("sumYP: total energy consumption of the action of leaving sleep mode on physical machines");
	    		System.out.print("sum Z: total energy consumption of execution of physical machines that on sleep mode ");
	    		System.out.print("sumE: total energy consumption of non-IT infrastructure in DCs");
	    		System.out.print("sumLambda: total energy consumption of migrations");
	    		System.out.println("\n");  
	    		System.out.println("total energy consumption(without renewable energy) -- ETotal: "+nf.format(sum));
	    		System.out.println("sumF : "+nf.format(sumF));
	    		double sumInDC=sumX+sumXP+sumY+sumYP+sumZ+sumE;	
	    		System.out.println("energy consumption in DC(without renewable energy, without migration) -- EInDC: "+nf.format(sumInDC));
	    		double reu1 = sum - cplex.getObjValue();  
	    		System.out.println("renewable energy used = ETotal - ObjValue: "+nf.format(reu1));
	    		double reu2 = sumInDC - sumF;
	    		System.out.println("renewable energy used = EInDC - sumF: "+nf.format(reu2));
	    		double reuTotal = 0.0;
	    		for(int i = 0; i < k; i++) {
	    			reuTotal = reuTotal + G[i];
	    		}
	    		System.out.println("renewable energy utilisation ratio : " + nt.format(reu1 / reuTotal));
				double[] CE = new double[k];
				double[] CMEps = new double[k];
				double[] Constant = new double[k];
				for (int i=0; i<k; i++) {
					CE[i] = Emin[i]; 
					CMEps[i] = M[i]*Eps;			
					Constant[i] = CE[i] + CMEps[i];
					System.out.println("constant energy used during a timeslot in DC " + i + ": " + nf.format(Constant[i]));
					System.out.println("constant energy used in DC " + i + ": " + nf.format(Constant[i]*n));
				}
    			double[] nbMigration = new double [k];
    			for (int i=0; i<k; i++) {
    				for (int tt=1; tt<=n; tt++) {
    					for (int j=0; j<k; j++) {
    						nbMigration[i] = nbMigration[i] + cplex.getValue(lambda[i][j][tt]);
    					}
    				}
    				System.out.println("number of requests migratied to other DC from DC " + i + ": " + nbMigration[i]);
    			}
    			double sumMigration = 0.0;
    			for (int i=0; i<k; i++) {
    				sumMigration = sumMigration + nbMigration[i];
    			}
    			System.out.println("total number of requests migratied : " + sumMigration);
	    		System.out.println("****************************************************over********************************************************");
	        }
	        else {
	        	System.out.println("problem not solved");
	        }
	        cplex.end();
		}
		catch (IloException exception) {
			exception.printStackTrace();
		}	
	}
}
