package mmc;

/**
 * class contains functions necessary for calculation of e+e- pair production losses
 */

public class Epairproduction extends CrossSections{

    public EpairContinuous c;
    public EpairStochastic s;

    public double vMax, vUp, vMin;
    protected Epairproduction epair;

    private Integral I;
    private int i;
    private double v;

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates internal references to p and m, to be called from subclasses
     */

    public Epairproduction(Epairproduction cros){
	p=cros.p;
	m=cros.m;
	this.cros=cros.cros;
	epair=cros;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes subclasses and creates internal references to p and m
     */

    public Epairproduction(CrossSections cros){
	super(cros);
	epair=this;
	c = new EpairContinuous(this);
	s = new EpairStochastic(this);
	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call before using the pair production functions to set the component of the primary
     */

    public void setEnergy(int i){
	double aux;
	cros.component=i;
	vMin=4*Me/p.e;
	vMax=1-(3./4)*sqrtE*(p.m/p.e)*Math.pow(m.Z[i], 1./3);
	aux=p.m/p.e;
	aux=1-6*aux*aux;
	vMax=Math.min(vMax, aux);
	vMax=Math.min(vMax, 1-p.m/p.e);
	if(vMax<vMin) vMax=vMin;
	vUp=Math.min(vMax, m.vCut(p.e));
	if(vUp<vMin) vUp=vMin;
    }

     //----------------------------------------------------------------------------------------------------//

    /**
     * this is the calculation of the d2Sigma/dvdRo - interface to Integral
     */

    public double function(double r){
	double Fe, Fm, Le, Lm, Ye, Ym, s, b, k, g1, g2;
	double aux, aux1, aux2, r2, Z3;
	r=1-r; // only for integral optimization - do not forget to swap integration limits!
	r2=r*r;
	Z3=Math.pow(m.Z[i], -1./3);
	aux=(p.m*v)/(2*Me);
	aux*=aux;
	s=aux*(1-r2)/(1-v);
	b=(v*v)/(2*(1-v));
	Ye=(5-r2+4*b*(1+r2))/(2*(1+3*b)*Math.log(3+1/s)-r2-2*b*(2-r2));
	Ym=(4+r2+3*b*(1+r2))/((1+r2)*(1.5+2*b)*Math.log(3+s)+1-1.5*r2);

	aux=(1.5*Me)/(p.m*Z3);
	aux*=aux;
	aux1=(1+s)*(1+Ye);
	aux2=(2*Me*sqrtE*m.B[i]*Z3)/(p.e*v*(1-r2));
	Le=Math.log((m.B[i]*Z3*Math.sqrt(aux1))/(1+aux2*aux1))-0.5*Math.log(1+aux*aux1);
	Lm=Math.log(((p.m/(1.5*Me))*m.B[i]*Z3*Z3)/(1+aux2*(1+s)*(1+Ym)));

	Fe=Le>0?(1/s<halfPrecision?(1.5-r2/2+b*(1+r2))/s:(((2+r2)*(1+b)+s*(3+r2))*Math.log(1+1/s)+(1-r2-b)/(1+s)-(3+r2)))*Le:0;
	Fm=Le>0?(((1+r2)*(1+1.5*b)-(1+2*b)*(1-r2)/s)*Math.log(1+s)+s*(1-r2-b)/(1+s)+(1+2*b)*(1-r2))*Lm:0;

	if(m.Z[i]==1){ g1=4.4e-5; g2=4.8e-5; }
	else{ g1=1.95e-5; g2=5.3e-5; }
	aux=p.e/p.m;
	k=0.058*Math.log(aux/(1+g2*aux/Z3))-0.14;
	if(k<=0) k=0;
	else k=(0.073*Math.log(aux/(1+g1*aux/(Z3*Z3)))-0.26)/k;
	if(k<0) k=0;

	aux=Alpha*Re;
	aux*=aux/(1.5*Pi);
	aux1=Me/p.m;
	aux1*=aux1;

	aux*=2*m.Z[i]*(m.Z[i]+k)*((1-v)/v)*lpm(r2, b, s)*(Fe+aux1*Fm);
	if(aux<0) aux=0;
	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    private boolean init=true;
    private double eLpm;

    /**
     * Landau Pomeranchuk Migdal effect evaluation
     */

    private double lpm(double r2, double b, double x){
	if(cros.lpm){
	    if(init){
		init=false;
		double sum=0;
		for(int i=0; i<m.num; i++) sum+=m.Z[i]*m.Z[i]*Math.log(3.25*m.B[i]*Math.pow(m.Z[i], -1./3));
		eLpm=p.m/(Me*Re);
		eLpm*=(eLpm*eLpm)*Alpha*p.m/(2*Pi*m.No*p.c*p.c*sum);
	    }
	    double A, B, C, D, E, s;
	    double s2, s36, s6, d1, d2, atan, log1, log2;
	    s=Math.sqrt(eLpm/(p.e*v*(1-r2)))/4;
	    s6=6*s; atan=s6*(x+1);
	    if(atan>1/computerPrecision) return 1;
	    s2=s*s; s36=36*s2;
	    d1=s6/(s6+1); d2=s36/(s36+1);
	    atan=Math.atan(atan)-Pi/2;
	    log1=Math.log((s36*(1+x)*(1+x)+1)/(s36*x*x));
	    log2=Math.log((s6*(1+x)+1)/(s6*x));
	    A=0.5*d2*(1+2*d2*x)*log1-d2+6*d2*s*(1+((s36-1)/(s36+1))*x)*atan;
	    B=d1*(1+d1*x)*log2-d1;
	    C=-d2*d2*x*log1+d2-(d2*d2*(s36-1)/(6*s))*x*atan;
	    D=d1-d1*d1*x*log2;
	    E=-s6*atan;
	    return ((1+b)*(A+(1+r2)*B)+b*(C+(1+r2)*D)+(1-r2)*E)/(((2+r2)*(1+b)+x*(3+r2))*Math.log(1+1/x)+(1-r2-b)/(1+x)-(3+r2));
	}
	else return 1;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is the calculation of the dSigma/dv
     */

    public double ePair(double v, int i){
	if(jt){
	    setEnergy(i);
	    if(v>=vUp) return Math.max(J[i].interpolate(p.e, Math.log(v/vUp)/Math.log(vMax/vUp)), 0);
	}
	double rMax, aux, aux2;
	this.i=i;
	this.v=v;
	aux=1-(4*Me)/(p.e*v);
	aux2=1-(6*p.m*p.m)/(p.e*p.e*(1-v));
	if(aux>0 && aux2>0) rMax=Math.sqrt(aux)*aux2;
	else rMax=0;
	aux=Math.max(1-rMax, computerPrecision);
	return m.No*m.n[i]*p.c*p.c*(I.integrateOpened(1-rMax, aux, this)+I.integrateWithLog(aux, 1, this));
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J[];
    public boolean jt=false;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double e, double v){
	p.setEnergy(e);
	setEnergy(cros.component);
	if(vUp==vMax) return 0;
	v=vUp*Math.exp(v*Math.log(vMax/vUp));
	return ePair(v, cros.component);
    }

}
