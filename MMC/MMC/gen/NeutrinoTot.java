package gen;
import mmc.*;

/**
 * Class contains functions for calculation of neutrino interaction cross sections.
 */

public class NeutrinoTot extends PhysicsModel{

    private NeutrinoInt N;
    private Integral I[];
    private double H[];

    private double E;
    private boolean cc, nu;
    private double rnd, rns;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Class constructor.
     */

    public NeutrinoTot(){
	N = new NeutrinoInt();
	I = new Integral[4];
        for(int i=0; i<4; i++) I[i] = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Charged and Neutral current interaction for neutrinos and antineutrinos
     */

    public double dSdy(double E, boolean cc, boolean nu){
	this.E=E;
	if(jt) return Math.max(Jo[(cc?1:0)+(nu?2:0)].interpolate(E), 0);
	this.cc=cc;
	this.nu=nu;
	return I[(cc?1:0)+(nu?2:0)].integrateOpened(0, 1, this);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Charged and Neutral current interaction for neutrinos and antineutrinos
     */

    public double dSdy(double E, boolean cc, boolean nu, double rnd){
	this.E=E;
	if(jt){
	    rns=rnd;
	    int i=(cc?1:0)+(nu?2:0);
	    return H[i]=Math.max(Jo[i].interpolate(E), 0);
	}
	this.cc=cc;
	this.nu=nu;
	return I[(cc?1:0)+(nu?2:0)].integrateOpened(0, 1, this, rnd);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns the energy transferred to the hardonic state
     */

    public double e(boolean cc, boolean nu){
	if(jt){
	    int i=(cc?1:0)+(nu?2:0);
	    return E*J[i].findLimit(E, rns*H[i]);
	}
	return E*I[(cc?1:0)+(nu?2:0)].getUpperLimit();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance cross section
     */

    public double dSdy(double E){
	return dSdy(E, 0);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance cross section
     */

    public double dSdy(double E, double rnd){
	this.E=E;
	this.rnd=rnd;
	double aux=1.e-6*Mw*Mw;
	double s=2.e-3*Me*E;
	aux*=aux/((s-aux)*(s-aux)+1.e-6*Gw*Gw*aux);
	aux*=9.e12*Gf*Gf*s/(3*Pi);
	s=Re*1.e-3*Me/Alpha;
	return N.RZ*s*s*aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Glashow resonance energy lost by neutrino
     */

    public double e(double m){
	double M=Math.sqrt(2.e-3*Me*E+1.e-6*Me*Me);
	double e=(M*M-1.e-6*m*m)/(2*M);
	return E-(e*((E+1.e-3*Me)/M)+e*(E/M)*(2*rnd-1));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Neutrino oscillation probability (mu->tau); units: E [GeV] L [m]
     */

    public double Pmt(double E, double L){
	double aux;
	aux=De2*L/E;
	aux*=(1.e-12*1.e2/1.e3)*(Alpha/(Re*Me))/4;
	aux=Math.sin(aux);
	aux*=aux;
	return St2*aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross section dSdy - interface to Integral.
     */

    public double function(double y){
	return N.RA*N.dSdy(y, E, cc, nu);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints the cross sections in [cm^2] as a function of energy in [GeV].
     */

    public static void main(String[] args){
        NeutrinoTot T = new NeutrinoTot();
	Output.raw=true;
	T.interpolate(".atmflux");
	for(double E=1.e1; E<1.e12; E*=1.05){
	    Output.out.println(Output.f(E)+" "+Output.f(T.dSdy(E, true, true)/T.N.RA)+" "+Output.f(T.dSdy(E, true, false)/T.N.RA)+" "+
			       Output.f(T.dSdy(E, false, true)/T.N.RA)+" "+Output.f(T.dSdy(E, false, false)/T.N.RA)+" "+Output.f(T.dSdy(E)/T.N.RZ));
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    public void interpolate(String name){
	boolean flag;
	name+=".gen_neutrino";

	if(nlow!=Me) name+="_l"+Output.f(nlow);
	if(ebig!=bigEnergy) name+="_b"+Output.f(ebig);

	if(Output.raw) name+="_raw"; else name+="_ascii";
	name+=".data";

	do {
	    if(Output.texi) return;
	    flag=false;
	    try{
		if(!Output.exists(name)){
		    CteqPDF.Ctq.interpolate();
		    N.interpolate();
		}
		Output.open(name);

		interpolate();
		Output.err.println("Finished parameterizations");

		Output.close();
	    }catch (mmcException error){
		flag=true;
		Output.delete(name);
	    }
	} while(flag);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    public void interpolate(){
	int g=5;
	double e_hi=ebig*1.e-3;
	double e_low=nlow*1.e-3;

	jt=false;
	H = new double[4];
	Output.err.print("Parameterizing neutrT ... ");
	J = new Interpolate[4];
	Jo = new Interpolate[4];
	for(int i=0; i<4; i++){
	    switch(i){
	    case 0: cc=false; nu=false; break;
	    case 1: cc=true;  nu=false; break;
	    case 2: cc=false; nu=true;  break;
	    case 3: cc=true;  nu=true;  break;
	    }
	    J[i] = new Interpolate(num1, e_low, e_hi, num1, 0, 1, this, g, false, false, true, g, false, false, false, g, true, false, false);
	    Jo[i] = new Interpolate(num1, e_low, e_hi, this, g, false, false, true, g, true, false, false);
	}
	jt=true;
	Output.err.println("done");
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J[];
    public boolean jt=false;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double E, double y){
	this.E=E;
	return I[(cc?1:0)+(nu?2:0)].integrateOpened(0, y, this);
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate Jo[];

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double E){
	return dSdy(E, cc, nu);
    }

}
