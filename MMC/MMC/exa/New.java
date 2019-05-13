package exa;
import mmc.*;

/**
 * Example of mmc package use.
 */

public class New{

    public static void main(String[] args){
	long i;
	double e, x, result;
	Propagate p = new Propagate("Frejus Rock", -1, 1.e-2);
	Output.DEBUG=false;
	Output.out.println(Output.version);
	// p.s.b.lorenzCut=1.e7;  p.s.b.lorenz=true;  // Enable Lorenz invariant violation - cutoff energy in [MeV]
	p.s.n.bb=1;            // Choose parametrization of the photon-nucleon cross section
	p.interpolate("all");  // Comment this out to leave out parameterizations
	for(i=0;i<100000;i++){
	    e=1.e8;
	    x=30000.;
	    Output.out.println("Propagating muon with energy "+Output.f(e)+" MeV to the distance of "+
			       Output.f(x)+" cm through "+p.m.name+" Step "+i);
	    result=p.propagateTo(x, e);
	    Output.out.println("returned value is e = "+Output.f(result)+"\n");
	}
    }

}
