package mmc;

/**
 * medium definition
 */

public class Medium extends PhysicsModel{

    public int num;                // number of components
    public double[] Z;             // nucleus charge
    public double[] A;             // atomic number
    public double[] n;             // number of atoms in a molecule
    public double totZ;            // sum of charges of all nuclei

    public double ZA;              // <Z/A>
    public double I;               // ionization potential [eV]
    public double C1, C, a;        // ionization formula constants
    public double m, X0, X1, d0;   // ionization formula constants (continued)
    public double r;               // refraction index

    public double B[];             // radiation logarithm constant B
    public double P[];             // radiation logarithm constant bPrime
    public double Xo;              // radiation length [cm]

    public double rho;             // multiplicative density correction factor
    public double Ro;              // mass density [g/cm3]
    public double No;              // molecule density [number/cm3]
    public double M[];             // average nucleon weight in a nucleus [MeV]
    public String E[];             // element name
    public String name;            // medium name

    protected double ecut;         // cutoff energy [MeV]
    protected double vcut;         // relative cutoff energy
    private double vCut;           // relative cutoff energy - call setCut(E) to set this

    private Integral L;            // Needed for calculation of the
    public double[] mN;            // Woods-Saxon potential factor
    public double MM;              // average all-component nucleon weight
    public double totA;            // sum of nucleons of all nuclei


    //----------------------------------------------------------------------------------------------------//

    /**
     * initialize medium by its name and proposed values of energy cut (ecut) and fractional energy cut (vcut).
     * The cuts ecut and vcut to be used must satisfy ecut&gt;0 and 0&lt;=vcut&lt;=1; If both values satisfy
     * these inequalities, the higher of E*vcut and ecut is used. If only one value satisfies the inequalities,
     * only that one value is used. If both values are outside these intervals, vcut=1 is assumed.
     */

    public Medium(String w, double ecut, double vcut, double rho){
	name=w;
	this.rho=rho>0?rho:1;
	if(w.equalsIgnoreCase("water")) initWater();
	else if(w.equalsIgnoreCase("ice")) initIce();
	else if(w.equalsIgnoreCase("salt")) initSalt();
	else if(w.equalsIgnoreCase("standard rock")) initStandardrock();
	else if(w.equalsIgnoreCase("frejus rock")) initFrejusrock();
	else if(w.equalsIgnoreCase("iron")) initIron();
	else if(w.equalsIgnoreCase("hydrogen")) initHydrogen();
	else if(w.equalsIgnoreCase("lead")) initLead();
	else if(w.equalsIgnoreCase("uranium")) initUranium();
	else if(w.equalsIgnoreCase("air")) initAir();
	else if(w.equalsIgnoreCase("mineral oil")) initParaffin();
	else if(w.equalsIgnoreCase("antares water")) initAntaresWater();
	else {
	    Output.err.println("Warning (in Medium/Medium): unknown medium: \""+w+"\", defaulting to water");
	    name="water";
	    initWater();
	}
	this.ecut=ecut;
	this.vcut=vcut;
	vCut=1.;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * calculate radiation length for simple elements
     */

    public double radl(double Z){
	switch(1){
	case 1: // Tsai
	    {
		double Lr, fZ, Lp, Z3, a=Alpha*Z; a*=a;
		fZ=a*(1/(1+a)+0.20206+a*(-0.0369+a*(0.0083-0.002*a)));
		switch((int)Math.round(Z)){
		case 1: Lr=5.31; Lp=6.144; break;
		case 2: Lr=4.79; Lp=5.621; break;
		case 3: Lr=4.74; Lp=5.805; break;
		case 4: Lr=4.71; Lp=5.924; break;
		default:
		    {
			Z3=Math.pow(Z, -1./3);
			Lr=Math.log(184.15*Z3);
			Lp=Math.log(1194.*Z3*Z3);
		    }
		}
		return Z*(Z*(Lr-fZ)+Lp);
	    }
	default: // Dahl
	    return Z*(Z+1)*Math.log(287./Math.sqrt(Z));
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialization of arrays
     */

    private void inita(int i){
	num=i;
	Z = new double[num];
	A = new double[num];
	n = new double[num];
	B = new double[num];
	P = new double[num];
	M = new double[num];
	E = new String[num];
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialization code, common to all media
     */

    private void initr(){
	int i;
	boolean flag=false;
	double aux1=0, aux2=0, aux3=0, aux4=0;
	Ro*=rho;
	for(i=0;i<num;i++){
	    aux1+=n[i]*Z[i];
	    aux2+=n[i]*A[i];
	    B[i]=elB(i);
	    P[i]=elP(i);
	    M[i]=(Z[i]*Mp+(A[i]-Z[i])*Mn)/A[i];
	    aux3+=n[i]*A[i]*M[i];
	    if(Z[i]!=1) flag=true;
	    aux4+=n[i]*radl(Z[i]);
	}
	totZ=aux1;
	totA=aux2;
	ZA=aux1/aux2;
	No=Ro*Na/aux2;
	MM=aux3/aux2;
	Xo=1/(4*Alpha*Re*Re*No*aux4);
	C1=2*Log10;
	r=1.31;                         // only for ice - change if needed (sea water: 1.35)
	ecut=Me/Math.sqrt(1-1/(r*r));   // in order to emit Cerenkov photons
	if(flag){
	    mN = new double[num];
	    L = new Integral(iromb, imaxs, iprec);
	    for(i=0;i<num;i++) if(Z[i]!=1){
		r0=Math.pow(A[i], 1./3);
		r0=1.12*r0-0.86/r0;
		mN[i]=1-4*Pi*0.17*L.integrateWithSubstitution(r0 , -1, this, 2)/A[i];
	    }
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * set the value of radiation logarithm constant B
     */

    private double elB(int i){
	int z=(int)Math.round(Z[i]);
	switch(z){
	case 1: return 202.4;
	case 2: return 151.9;
	case 3: return 159.9;
	case 4: return 172.3;
	case 5: return 177.9;
	case 6: return 178.3;
	case 7: return 176.6;
	case 8: return 173.4;
	case 9: return 170.0;
	case 10: return 165.8;
	case 11: return 165.8;
	case 12: return 167.1;
	case 13: return 169.1;
	case 14: return 170.8;
	case 15: return 172.2;
	case 16: return 173.4;
	case 17: return 174.3;
	case 18: return 174.8;
	case 19: return 175.1;
	case 20: return 175.6;
	case 21: return 176.2;
	case 22: return 176.8;
	case 26: return 175.8;
	case 29: return 173.1;
	case 32: return 173.0;
	case 35: return 173.5;
	case 42: return 175.9;
	case 50: return 177.4;
	case 53: return 178.6;
	case 74: return 177.6;
	case 82: return 178.0;
	case 92: return 179.8;
	default: return 182.7;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * set the value of radiation logarithm constant bPrime
     */

    private double elP(int i){
	int z=(int)Math.round(Z[i]);
	switch(z){
	case 1: return 446;
	default: return 1429;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize water
     */

    private void initWater(){
	inita(2);
	E[0]="H";
	E[1]="O";
	Z[0]=1; // H
	Z[1]=8; // O
	A[0]=1.00794;
	A[1]=15.9994;
	n[0]=2;
	n[1]=1;
	I=75.0;
	C=-3.5017;
	a=0.09116;
	m=3.4773;
	X0=0.2400;
	X1=2.8004;
	d0=0;
	Ro=1.000;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize ice
     */

    private void initIce(){
	inita(2);
	E[0]="H";
	E[1]="O";
	Z[0]=1; // H
	Z[1]=8; // O
	A[0]=1.00794;
	A[1]=15.9994;
	n[0]=2;
	n[1]=1;
	I=75.0;
	C=-3.5017;
	a=0.09116;
	m=3.4773;
	X0=0.2400;
	X1=2.8004;
	d0=0;
	Ro=0.917;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize salt (added by Ped)
     */

    private void initSalt(){
	inita(2);
	E[0]="Na";
	E[1]="Cl";
	Z[0]=11; // H
	Z[1]=17; // O
	A[0]=22.98977;
	A[1]=35.4527;
	n[0]=1;
	n[1]=1;
	I=175.3;  // Calculated by ESTAR detabase (it could be 185 eV by the method of reference below)
	// C through X1 are based on Atomic Data and Nuclear Data Tables 78, 183-356 (2001), Appendix A
	C=-4.5041;
	a=0.1632;
	m=3;
	X0=0.2;
	X1=3.0;
	d0=0;
	Ro=2.323; // Solid halite density
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize standard rock
     */

    private void initStandardrock(){
	inita(1);
	E[0]="Standard Rock";
	Z[0]=11; // Ionization potential and density corrections
	A[0]=22; // are close to those of calcium carbonate
	n[0]=1;
	I=136.4;
	C=-3.7738;
	a=0.08301;
	m=3.4120;
	X0=0.0492;
	X1=3.0549;
	d0=0;
	Ro=2.650;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize standard rock
     */

    private void initFrejusrock(){
	inita(1);
	E[0]="Frejus Rock";
	Z[0]=10.12;
	A[0]=20.34;
	n[0]=1;
	I=149.0;
	C=-5.053;
	a=0.078;
	m=3.645;
	X0=0.288;
	X1=3.196;
	d0=0;
	Ro=2.740;
	initr();
    }


    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize iron
     */

    private void initIron(){
	inita(1);
	E[0]="Fe";
	Z[0]=26;
	A[0]=55.845;
	n[0]=1;
	I=286.0;
	C=-4.2911;
	a=0.14680;
	m=2.9632;
	X0=-0.0012;
	X1=3.1531;
	d0=0.12;
	Ro=7.874;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize hydrogen
     */

    private void initHydrogen(){
	inita(1);
	E[0]="H";
	Z[0]=1;
	A[0]=1.00794;
	n[0]=1;
	I=21.8;
	C=-3.0977;
	a=0.13483;
	m=5.6249;
	X0=0.4400;
	X1=1.8856;
	d0=0;
	Ro=0.07080;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize lead
     */

    private void initLead(){
	inita(1);
	E[0]="Pb";
	Z[0]=82;
	A[0]=207.2;
	n[0]=1;
	I=823.0;
	C=-6.2018;
	a=0.09359;
	m=3.1608;
	X0=0.3776;
	X1=3.8073;
	d0=0.14;
	Ro=11.350;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize uranium
     */

    private void initUranium(){
	inita(1);
	E[0]="U";
	Z[0]=92;
	A[0]=238.0289;
	n[0]=1;
	I=890.0;
	C=-5.8694;
	a=0.19677;
	m=2.8171;
	X0=0.2260;
	X1=3.3721;
	d0=0.14;
	Ro=18.950;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize air
     */

    private void initAir(){
	final double fr1=2*78.1;
	final double fr2=2*21.0;
	final double fr3=0.9;
	final double fra=fr1+fr2+fr3;
	inita(3);
	E[0]="N";
	E[1]="O";
	E[2]="Ar";
	Z[0]=7; // N
	Z[1]=8; // O
	Z[2]=18; // Ar
	A[0]=14.0067;
	A[1]=15.9994;
	A[2]=39.948;
	n[0]=fr1/fra;
	n[1]=fr2/fra;
	n[2]=fr3/fra;
	I=85.7;
	C=-10.5961;
	a=0.10914;
	m=3.3994;
	X0=1.7418;
	X1=4.2759;
	d0=0;
	Ro=1.205e-3; // dry, 1 atm
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize mineral oil or paraffin CH3(CH2)~23CH3 (added by Ped)
     */

    private void initParaffin(){
	inita(2);
	E[0]="C";
	E[1]="H";
	Z[0]=6; // C
	Z[1]=1; // H
	A[0]=12.0011;
	A[1]=1.0079;
	n[0]=25;
	n[1]=52;
	I=55.9;
	C=-2.9551;
	a=0.1209;
	m=3.4288;
	X0=0.1289;
	X1=2.5084;
	d0=0;
	Ro=0.93;
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * initialize ANTARES water
     * Sea water (Mediterranean Sea, ANTARES place)
     * ==========================================================================
     * WATER DENSITY CHANGES WITH THE DEPTH FROM 1.0291 g/cm^3 AT SURFACE
     * UP TO 1.0404 g/cm^3 AT THE SEA BED
     * (ANTARES-Site/2000-001 and references therein)
     *
     * The error which is caused by this simplified approach (average value for
     * density) does not exceed 0.5% (much less, in fact) that is comparable with
     *  an error which comes from uncertainties with the muon cross-sections.
     *==========================================================================
     */
    private void initAntaresWater(){ // added by Claudine Colnard,
	                             // Institute Nikhef, The Netherlands,
	                             // ANTARES collaboration.
	inita(8);
	E[0]="H";
	E[1]="O";
	E[2]="Na";
	E[3]="K";
	E[4]="Mg";
	E[5]="Ca";
	E[6]="Cl";
	E[7]="S";
	Z[0]=1;  // H
	Z[1]=8;  // O
	Z[2]=11; // Na
	Z[3]=19; // K
	Z[4]=12; // Mg
	Z[5]=20; // Ca
	Z[6]=17; // Cl
	Z[7]=16; // S
	A[0]=1.008;   //  Chemical composition of the seawater
	A[1]=15.999;  //              according to
	A[2]=22.99;   //  A.Okada, Astropart. Phys. 2 (1994) 393
	A[3]=39.10;   //         and references therein
	A[4]=24.31;   //  corrected for Mediterranean Sea, ANTARES place
	A[5]=40.08;   //  according to salinity  38.44+-0.02 g/kg,
	A[6]=35.45;   //  as cited in J.Brunner, ANTARES-Site/2000-001
	A[7]=32.07;   //  instead of 35.0 g/kg as cited in A.Okada, ...
	n[0]=2;       //  (so, n[2-7] have been just multiplied by 1.098)
	n[1]=1.00884;
	n[2]=0.00943;
	n[3]=0.000209;
	n[4]=0.001087;
	n[5]=0.000209;
	n[6]=0.01106;
	n[7]=0.00582;
	I = 75.0;     // All the same as for pure water
	C = -3.5017;
	a = 0.09116;
	m = 3.4773;
	X0 = 0.2400;
	X1 = 2.8004;
	d0 = 0;
	Ro = 1.03975; // J.Brunner, ANTARES-Site/2000-001, the mean value
                      // for sea water density at the ANTARES place between
                      // sea bed D = 2400 m (1.0404 g/cm^3) and middle of
	              // detector D = 2126 m (1.0391 g/cm^3).
	              // Ro = 1.0341
	              // for sea water density at the ANTARES place between
                      // surface D = 0 m (1.0291 g/cm^3) and middle of
                      // detector D = 2126 m (1.0391 g/cm^3)
	initr();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * return the value of the fractional energy cut, as described in the help for the constructor
     */

    public double vCut(double E){
	double aux;
	if(ecut>0){
	    aux=ecut/E;
	    if(vcut>0 && vcut<=1){
		if(aux>vcut) vCut=aux;
		else vCut=vcut;
	    }
	    else vCut=aux;
	}
	else{
	    if(vcut>0 && vcut<=1) vCut=vcut;
	    else vCut=1.;
	}
	return vCut;
    }

    //----------------------------------------------------------------------------------------------------//

    private double r0;

    /**
     * Woods-Saxon potential calculation - interface to Integral
     */

    public double function(double r){
	final double a=0.54;
	return r*r/(1+Math.exp((r-r0)/a));
    }

}
