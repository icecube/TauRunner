package mmc;

/**
 * class contains functions necessary for calculation of photonuclear losses
 */

public class PhotoContinuous extends Photonuclear{

    private Integral I;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    public PhotoContinuous(Photonuclear cros){
	super(cros);
	I = new Integral(iromb, imaxs, iprec);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * photonuclear energy losses - interface to Integral
     */

    public double function(double v){
	return v*photo.photoN(v, cros.component);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Contribution of photonuclear interactions to -dE/dx
     */

    public double dEdx(){
	if(cros.cp<=0) return 0;
	if(jt) return Math.max(J.interpolate(p.e), 0);
	int i;
	double sum;
	sum=0;
	for(i=0; i<m.num; i++){
	    setEnergy(i);
	    sum+=I.integrateWithLog(vMin, vUp, this);
	}
	return cros.cp*p.e*sum;
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
