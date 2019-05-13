package mmc;

/**
 * class contains functions necessary for calculation of bremsstrahlung losses
 */

public class BremsContinuous extends Bremsstrahlung{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public BremsContinuous(Bremsstrahlung cros){
	super(cros);
 	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * bremsstrahlung energy losses - interface to Integral
     */

    public double function(double v){
	return v*brems.Sel(v, cros.component);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * contribution of bremsstrahlung to -dE/dx
     */

    public double dEdx(){
	if(cros.cb<=0) return 0;
	if(jt) return Math.max(J.interpolate(p.e), 0);
	int i;
	double sum;
	sum=0;
	for(i=0; i<m.num; i++){
	    setEnergy(i);
	    sum+=I.integrateOpened(0, vUp, this);
	}
	return cros.cb*p.e*sum;
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
