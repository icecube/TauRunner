package mmc;

/**
 * class contains functions necessary for calculation of ionization losses
 */

public class Ionizationloss extends CrossSections{

    public IonizContinuous c;
    public IonizStochastic s;

    public double vMax, vUp, vMin;
    protected double beta, gamma;
    protected Ionizationloss ioniz;

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates internal references to p and m, to be called from subclasses
     */

    public Ionizationloss(Ionizationloss cros){
	p=cros.p;
	m=cros.m;
	this.cros=cros.cros;
	ioniz=cros;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes subclasses and creates internal references to p and m
     */

    public Ionizationloss(CrossSections cros){
	super(cros);
	ioniz=this;
	c = new IonizContinuous(this);
	s = new IonizStochastic(this);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * call before using the ionization functions, but after setting the particle energy
     */

    public void setEnergy(){
	double aux;
	beta=p.p/p.e;
	gamma=p.e/p.m;
	vMin=(1.e-6*m.I)/p.e;
	aux=Me/p.m;
	vMax=2*Me*(gamma*gamma-1)/((1+2*gamma*aux+aux*aux)*p.e);
	vMax=Math.min(vMax, 1-p.m/p.e);
	if(vMax<vMin) vMax=vMin;
	vUp=Math.min(vMax, m.vCut(p.e));
	if(vUp<vMin) vUp=vMin;
	if(ioniz!=this){
	    ioniz.beta=beta;
	    ioniz.gamma=gamma;
	    ioniz.vMax=vMax;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is what d2N/dv,dx is equal to in the first approximation
     */

    protected double d2Ndvdx(double v){
	double result, aux, aux2;
	aux=beta*beta;
	aux2=v/(1+1/gamma);
	aux2*=0.5*aux2;
	result=1-aux*(v/vMax)+aux2;
	result*=K*p.c*p.c*m.ZA/(2*aux*p.e*v*v);
	return m.Ro*result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * this is the inelastic electron bremsstrahlung correction to the first approximation of the d2N/dv,dx
     */

    protected double inelCorrection(double v){
	double result, a, b, c;
	a=Math.log(1+2*v*p.e/Me);
	b=Math.log((1-v/vMax)/(1-v));
	c=Math.log((2*gamma*(1-v)*Me)/(p.m*v));
	result=a*(2*b+c)-b*b;
	return (Alpha/(2*Pi))*result;
    }

}
