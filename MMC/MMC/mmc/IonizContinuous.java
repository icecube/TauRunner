package mmc;

/**
 * class contains functions necessary for calculation of ionization losses
 */

public class IonizContinuous extends Ionizationloss{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public IonizContinuous(Ionizationloss cros){
	super(cros);
 	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * inelastic electron bremsstrahlung correction to dEdx - interface to Integral
     */

    public double function(double v){
	return v*ioniz.d2Ndvdx(v)*ioniz.inelCorrection(v);
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * density correction
     */

    private double delta(){
	double X;
	X=Math.log(beta*gamma)/Log10;
	if(X<m.X0) return m.d0*Math.pow(10, 2*(X-m.X0));
	else if(X<m.X1) return m.C1*X+m.C+m.a*Math.pow(m.X1-X, m.m);
	else return m.C1*X+m.C;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * contribution of ionization to -dE/dx
     */

    public double dEdx(){
	if(cros.ci<=0) return 0;
	if(jt) return Math.max(J.interpolate(p.e), 0);
	double result, aux;
	setEnergy();
	aux=beta*gamma/(1.e-6*m.I);
	result=Math.log(vUp*(2*Me*p.e))+2*Math.log(aux);
	aux=vUp/(2*(1+1/gamma));
	result+=aux*aux;
	aux=beta*beta;
	result-=aux*(1+vUp/vMax)+delta();
	if(result>0) result*=K*p.c*p.c*m.ZA/(2*aux);
	else result=0;
	return cros.ci*(m.Ro*result+p.e*(I.integrateWithLog(vMin, vUp, this)));
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J;
    public boolean jt=false;

    /**
     * 1d parametrization - interface to Interpolate
     */

    public double functionInt(double e){
	p.setEnergy(e);
	return dEdx();
    }

}
