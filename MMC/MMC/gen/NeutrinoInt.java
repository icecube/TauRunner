package gen;
import mmc.*;

/**
 * Class contains functions for calculation of neutrino interaction differential cross sections.
 */

public class NeutrinoInt extends PhysicsModel{

    public double RA;
    public double RZ;

    private CteqPDF F;
    private Medium m;
    private double M, d1, d2, d3;

    private String name="Standard Rock";
    private double y, E;
    private boolean cc, nu;
    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Class constructor. Sets up the values of cross section parameters and calls the CteqPDF constructor.
     */

    public NeutrinoInt(){
	double d12, d22, d32, d42;
	d12=1/2.-2*Xw/3; d12*=d12;
	d22=-1/2.+Xw/3;  d22*=d22;
	d32=-2*Xw/3;     d32*=d32;
	d42=Xw/3;        d42*=d42;
	d1=d12+d22+d32+d42;
	d2=d12+d32-d22-d42;
	d3=d12+d22-d32-d42;
	F = new CteqPDF();
	m = new Medium(name, -1, 1, 1);
	RA=m.No*m.totA/m.Ro;
	RZ=m.No*m.totZ/m.Ro;
	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino (nu=true) and anti neutrino (nu=false) charged (cc=true) and neutral (cc=false) current
     * interaction cross sections as functions of q(x,y) and y and energy E in [GeV].
     */

    private double dS2dxdy(double q, double y, double E, boolean cc, boolean nu){
	int inu;
	inu=nu?1:-1;

	double aux, sum, x, Q2, qa, F2, F3;
	Q2=q*q;
	qa=1.e-3*m.MM*y*E;
	x=q*q/(2*qa);
	if(x>1) return 0;

	if(cc){
	    F2=F.PDF(2, x, q);
	    F3=F.PDF(1, x, q)+inu*F.PDF(3, x, q);
	}
	else{
	    F2=F.PDF(2, x, q)*d1-F.PDF(3, x, q)*d2;
	    F3=F.PDF(1, x, q)*d3;
	}

	aux=1.e-3*(cc?Mw:Mz);
	aux*=aux;
	sum=aux/(Q2+aux);
	sum*=sum*(1.e9*Gf*Gf*m.MM*E/Pi);
	aux=y-y*y/2;
	sum*=(1-aux)*F2+inu*aux*F3;
	aux=1.e-3*Re*Me/Alpha;
	aux*=aux;

	aux*=sum*q/qa;
	return Math.max(aux, 0);  // multiply by (m.No*m.totA) to get 1/length
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross sections integrated over x as functions of y and energy E in [GeV].
     */

    public double dSdy(double y, double E, boolean cc, boolean nu){
	if(jt) return Math.max(J[(cc?1:0)+(nu?2:0)].interpolate(E, y), 0);
	this.y=y;
	this.E=E;
	this.cc=cc;
	this.nu=nu;
	double qmax=Math.sqrt(2.e-3*m.MM*y*E);
	return I.integrateWithLog(Math.min(F.Qn, qmax), Math.min(F.Qm, qmax), this);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * returns the neutrino cross section dS2dxdy - interface to Integral.
     */

    public double function(double x){
	return dS2dxdy(x, y, E, cc, nu);
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
	Output.err.print("Parameterizing neutrI ... ");
	J = new Interpolate[4];
	for(int i=0; i<4; i++){
	    switch(i){
	    case 0: cc=false; nu=false; break;
	    case 1: cc=true;  nu=false; break;
	    case 2: cc=false; nu=true;  break;
	    case 3: cc=true;  nu=true;  break;
	    }
	    J[i] = new Interpolate(num1, e_low, e_hi, num1, 0, 1, this, g, false, false, true, g, false, false, false, g, true, false, false);
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
	return dSdy(y, E, cc, nu);
    }

}
