package gen;
import mmc.*;

/**
 * Class defines Gaisser-like muon and muon- and electron-neutrino energy spectra with cos*, muon energy loss, and decay corrections.
 */

public class IntFlux extends PhysicsModel{

    private double z0=0;
    private double X0=1.148;
    private double R=EarthModel.R0;
    public double D=0;
    public double Rc=0;

    private int P=0;
    private int M=0;
    private int S=0;

    public static int sM=0;
    public static double Efr=10;
    public static double mF=0;
    public static double gD=0.2;
    public static double gEcut=0;
    public static EarthModel eM;
    public static NeutrinoTot nT;

    private double xpar[][]={
	/* Standard US Atmosphere average parameters */
	{ 0.00851164, 0.124534, 0.059761, 2.32876, 19.279 },
	{ 0.102573, -0.068287, 0.958633, 0.0407253, 0.817285 },
	{ -0.017326, 0.114236, 1.15043, 0.0200854, 1.16714 },
	{ 1.3144, 50.2813, 1.33545, 0.252313, 41.0344 }
    };

    private double fpar[][]={
	{0.701, 2.715, 1, 1, 1},
	{0.340, 2.720, 1, 1, 1},
	{0.367, 2.712, 1, 1, 1},
	{0.646, 2.684, 1, 1, 1},
	{0.352, 2.680, 1, 1, 1},
	{0.310, 2.696, 1, 1, 1},
	{0.828, 2.710, 1, 1, 1},
	{0.465, 2.719, 1, 1, 1},
	{0.472, 2.741, 1, 1, 1}
	// { 0.701459, 2.71461, 1, 1, 1 }; all corrections
	// { 0.590968, 2.69506, 1, 0, 0 }; just cos* correction
	// { 0.594662, 2.69671, 0, 0, 0 }; original formula
    };

    private double dG=0;
    private double xx=0, xs=0, xo=0, xl=0, xn=0;
    private String name;
    private ExpCorr eC;
    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with particle type (all mu/mu-/mu+/all nu_mu/nu_mu/~nu_mu/all nu_e/nu_e/~nu_e/),
     * model and ground elevation z0 in [km]. flag switches between two cos* calculation algorithms.
     * h0 is the average production height in km or in g/cm^2 if negative.
     */

    public IntFlux(int type, int model, int flag, double h0, double z0){
	I = new Integral(iromb, imaxs, iprec);
	P=type-1;
	M=model;
	S=flag;
	if(M>=0){
	    eC = new ExpCorr(model, h0, z0);
	    X0=eC.X0;
	    name=model+"_"+Output.f(h0)+"_"+Output.f(z0);
	}
	R+=z0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize the class for spectrum index correction dG and prompt to pion muon component ratio Rc.
     */

    public IntFlux(int type, int model, int flag, double h0, double z0, double dG, double Rc, double D){
	this(type, model, flag, h0, z0);
	this.dG=dG;
	this.Rc=Rc;
	this.D=D;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrize h(x) c(x) o(x) l(x).
     */

    public void interpolate(String name){
	if(M<0) return;
	int g=5;
	boolean flag;
	name+=".gen_hcolxc_"+this.name;
	if(Output.raw) name+="_raw"; else name+="_ascii";
	name+=".data";

	do {
	    if(Output.texi) return;
	    flag=false;
	    try{
		Output.open(name);

		je=false;
		Output.err.print("Parameterizing hcolxC ... ");
		Je = new Interpolate[4];
		for(int i=0; i<4; i++){
		    hcol=i;
		    Je[i] = new Interpolate(num1, 0, 1, this, g, false, false, false, g, false, false, false);
		}
		je=true;
		Output.err.println("done");
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
     * Calculates 3 example distrbutions.
     */

    public static void main(String[] args){
	int w=1, p=0, f=0, m=-1;
	double h0=-114.8, z0=0;

	for(int n=0; n<args.length; n++){
	    if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
		Output.out.println("\n"+
"This program calculates Energy Integral of atm. lepton fluxes\n"+
"                       -w=[1-3] which program to run\n"+
"                       -t=[1-9] particle mu-+/num-+/nue-+\n"+
"                       -f=[0-1] method of cos* correction\n"+
"                       -m=[-1/0-3] choose atmosphere model\n"+
"                       -z0=[ground elevation in km]\n"+
"                       -h0=[average production height in km\n"+
"                                   or in g/cm^2 if negative]\n");
		return;
                }
	    else if(args[n].startsWith("-t=")){
		try{
		    p=(int)Double.parseDouble(args[n].substring(3));
		}catch(Exception error){
		    p=1;
		}
	    }
	    else if(args[n].startsWith("-w=")){
		try{
		    w=(int)Double.parseDouble(args[n].substring(3));
		}catch(Exception error){
		    w=1;
		}
	    }
	    else if(args[n].startsWith("-f=")){
		try{
		    f=(int)Double.parseDouble(args[n].substring(3));
		}catch(Exception error){
		    f=0;
		}
	    }
	    else if(args[n].startsWith("-m=")){
		try{
		    m=(int)Double.parseDouble(args[n].substring(3));
		}catch(Exception error){
		    m=-1;
		}
	    }
	    else if(args[n].startsWith("-z0=")){
		try{
		    z0=Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
		    z0=0;
		}
	    }
	    else if(args[n].startsWith("-h0=")){
		try{
		    h0=Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
		    h0=-114.8;
		}
	    }
	}

	if(h0==0) m=-1;
	if(w<1 || w>3) w=1;
	if(p<1 || p>9) p=1;
	if(f<0 || f>1) f=0;
	if(m<-1 || m>3) m=-1;
	Output.err.println("Choosing w="+w+" t="+p+" f="+f+" m="+m+" z0="+Output.f(z0)+" h0="+Output.f(h0));


	double x;
	IntFlux F = new IntFlux(p, m, f, h0, z0);


	switch(w){
	case 1:
	    for(int i=0; i<=1000; i++){
		x=i/1000.;
		Output.out.println(Output.f(x)+" "+Output.f(F.getIntFlux(600, x)));
	    }
	    break;
	case 2:
	    x=1;
	    Output.err.println("At cos="+Output.f(x)+", total flux is "+Output.f(F.getIntFlux(600, x))
			       +", mass overburden is "+Output.f((F.o(x)-F.X0)));
	    for(int i=0; i<=1000000; i++){
		Output.out.println(Output.f(F.getE(600, x, Math.random())));
	    }
	    break;
	case 3:
	    for(double e=5;e<5.e5;e*=1.01) Output.out.println(Output.f(e)+" "+Output.f(F.getfl(e, 1)));
	default:
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates integral flux above energy E in [GeV] at x=cos(zenith angle).
     */

    public double getIntFlux(double E, double x){
	return getE(E, x, -1);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates particle energy above threshold E0 in [GeV] at x=cos(zenith angle) as a function of a random number rnd.
     */

    public double getE(double E0, double x, double rnd){
	if(rnd<0){
	    if(jt) return J.interpolate(x, E0);
	    flSet(x);
	    return I.integrateWithSubstitution(1, -E0, this, -2);
	}
	else{
	    if(jt) return J.findLimit(x, rnd*J.interpolate(x, E0));
	    flSet(x);
	    I.integrateWithSubstitution(1, -E0, this, -2, rnd);
	    return -I.getUpperLimit();
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates differential flux for energy E in [GeV] at x=cos(zenith angle).
     */

    public double getfl(double E, double x){
	flSet(x);
	return function(-E);
    }

    //----------------------------------------------------------------------------------------------------//

    private void flSet(double x){
	if(mF>0 && P>=3){
	    eM.sett(Math.sqrt(1-x*x), 0, -x);
	    xn=eM.X(eM.ti, eM.tf)*mF;
	}
	double s=(1-D/R)*Math.sqrt(1-x*x);
	x=Math.sqrt(1-s*s);

	xx=x;
	switch(S){
	case 1:
	    xs=Math.sqrt(1-x*x)/(1+h(x)/R);
	    xs=Math.sqrt(1-xs*xs);
	    break;
	case 0:
	default:
	    xs=c(x);
	}
	if(P<3){
	    xo=fpar[P][3]*o(x)-X0; if(xo<0) xo=0;
	    xl=fpar[P][4]*l(x)*1.e2*Mmu/(Lmu*C);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    private double h(double x){
	if(M>=0) return je?Math.max(Je[0].interpolate(x), 0):eC.getExp(1, x);
	double s=Math.sqrt(1-x*x);
	double y=s/(1+xpar[0][0]);
	return xpar[0][4]*(1+xpar[0][2]*Math.pow(s, xpar[0][3]))/
	    Math.pow(1-(y*y), xpar[0][1]);
    }

    //----------------------------------------------------------------------------------------------------//

    private double c(double x){
	if(M>=0) return je?Math.max(Je[1].interpolate(x), 0):1/eC.getExp(2, x);
	return Math.sqrt((x*x+fpar[P][2]*
			  (xpar[1][0]*xpar[1][0]+
			   xpar[1][1]*Math.pow(x, xpar[1][2])+
			   xpar[1][3]*Math.pow(x, xpar[1][4])))/
			 (1+fpar[P][2]*(xpar[1][0]*xpar[1][0]+
					xpar[1][1]+xpar[1][3])));
    }

    //----------------------------------------------------------------------------------------------------//

    private double o(double x){
	if(M>=0) return je?Math.max(Je[2].interpolate(x), 0):eC.C.getX(x)/1.e2;
	return 1/(xpar[2][0]+xpar[2][1]*Math.pow(x, xpar[2][2])+
		  xpar[2][3]*Math.pow(1-x*x, xpar[2][4]));
    }

    //----------------------------------------------------------------------------------------------------//

    private double l(double x){
	if(M>=0) return je?Math.max(Je[3].interpolate(x), 0):eC.getExp(3, x);
	return 1e3/(xpar[3][0]+xpar[3][1]*Math.pow(x, xpar[3][2])+
				  xpar[3][3]*Math.pow(1-x*x, xpar[3][4]));
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Energy distribution as a function of energy x in [GeV]; zenith angle is set in a private variable - interface to Integral.
     */

    public double function(double x){
	double ei, ef, xn, sm;
	ef=-x;

	if(mF>0 && P>=3){
	    boolean nu=(P==4 || P==7);
	    double tI=(P==8)?nT.dSdy(ef):0;
	    tI+=nT.dSdy(ef, true, nu);
	    tI+=nT.dSdy(ef, false, nu);
	    xn=this.xn*tI;
	    if(xn>1) xn=1;
	}
	else xn=1;

	if(sM>0){
	    final double Nmax=2465;
	    double N;
	    switch(sM){
	    case 1: N=2445; break;
	    case 2: N=2300; break;
	    case 3: N=2115; break;
	    default: N=Math.min(sM, Nmax-1);
	    }
	    sm=Math.exp(-(1.15+14.9*Math.pow(1-N/Nmax, 1.12))/(0.97+ef*Efr));
	}
	else sm=1;

	if(gEcut>0) sm/=1+Math.exp(-Math.log(ef/(gEcut*Efr))/(gD*Log10));

	switch((int)Math.floor(P/3)){
	case 1:
	    {
		final double E1=121;
		final double E2=897;
		final double pK=0.213;
		final double A=2.85e-2;
		ei=ef;
		return sm*xn*A*fpar[P][0]*Math.pow(ei,-(fpar[P][1]+dG))*(1/(1+6*ei*xs/E1)+pK/(1+1.44*ei*xs/E2));
	    }
	case 2:
	    {
		final double E1=121;
		final double E2=897;
		final double E3=194;
		final double A=2.4e-3;

		double xi, a, b;
		ei=ef;

		if(xx<0.3){
		    a=0.11-2.4*xx;
		    b=-0.22+0.69*xx;
		}
		else{
		    a=-0.46-0.54*xx;
		    b=-0.01+0.01*xx;
		}
		xi=a+b*Math.log(ei)/Log10;
		return sm*xn*A*fpar[P][0]*Math.pow(ei,-(fpar[P][1]+dG))*(0.05/(1+1.5*ei*xs/E2)+0.185/(1+1.5*ei*xs/E3)+
								   11.4*Math.pow(ei, xi)/(1+1.21*ei*xs/E1));
	    }
	case 0:
	default:
	    {
		final double E1=115;
		final double E2=850;
		final double pK=0.054;
		final double a=0.260;
		final double b=0.360e-3;
		final double as=-0.522;
		final double bs=8.893e-3;
		final double A=0.14;

		double zfi, zff, p, aux, aux1, aux2, auf1, auf2, f1, f2, f3, f4, de;
		ei=((a+b*ef)*Math.exp(b*xo)-a)/b;

		aux1=(bs*a-as*b);
		aux1*=2*aux1;
		aux2=2*b*b*b;
		zfi=(bs*ei*b*(-2*bs*a+4*as*b+bs*ei*b)+aux1*Math.log(a+ei*b))/aux2;
		zff=(bs*ef*b*(-2*bs*a+4*as*b+bs*ef*b)+aux1*Math.log(a+ef*b))/aux2;

		aux=1.1*xs;
		auf1=aux/E1;
		auf2=aux/E2;
		aux1=1/(1+auf1*ei);
		aux2=1/(1+auf2*ei);
		aux=aux1+pK*aux2;

		f1=Math.pow(ei, -fpar[P][1]-1)*aux;
		f2=Math.pow(ei, -fpar[P][1]-1)*(-(fpar[P][1]+1)*aux/ei-auf1*aux1*aux1-pK*auf2*aux2*aux2);
		f3=(-fpar[P][1]-1)*f2/ei+Math.pow(ei, -fpar[P][1]-1)
		    *((fpar[P][1]+1)*aux/(ei*ei)
		      +((fpar[P][1]+1)/ei)*(auf1*aux1*aux1+pK*auf2*aux2*aux2)
		      +2*auf1*auf1*aux1*aux1*aux1+2*pK*auf2*auf2*aux2*aux2*aux2);
		aux=f2/f1;
		f4=f3/f1-aux*aux;

		aux1=as+bs*ef;
		aux2=as+bs*ei;
		de=Math.exp(b*xo)-(f2/(2*f1))*(aux1*aux1-aux2*aux2)/(a+b*ef)-(f4/2)*(zff-zfi);

		ei=ei-(f2/(2*f1))*(zff-zfi);

		p=Math.exp(-xl/ei);
		return sm*A*fpar[P][0]*p*de*Math.pow(ei,-(fpar[P][1]+dG))*(Rc<0?-Rc:(1/(1+1.1*ei*xs/E1)+pK/(1+1.1*ei*xs/E2)+Rc));
	    }
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Parametrizes the integral of this class.
     */

    public void interpolate(double E0){
	int g=5;
	double e_hi=bigEnergy*1.e-3;
	double e_low=E0;
	jt=false;
	J = new Interpolate(num1, 0, 1, num1, e_low, e_hi, this, g, false, false, false, g, false, false, true, g, false, false, true);
	jt=true;
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J;
    public boolean jt=false;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double x, double E){
	flSet(x);
	return I.integrateWithSubstitution(1, -E, this, -2);
    }

    //----------------------------------------------------------------------------------------------------//

    private int hcol;
    public static Interpolate Je[];
    public static boolean je=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double x){
	switch(hcol){
	case 0: return h(x);
	case 1: return c(x);
	case 2: return o(x);
	case 3: return l(x);
	default: return 0;
	}
    }

}
