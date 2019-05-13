package exa;
import mmc.*;

/**
 * Example of use with Symphony.
 */

public class MPMCWJ{
    public String results="";
    public double dummy=-1.0;

    public void start(){
	double e, x, result;
	Propagate p = new Propagate("Frejus Rock", -1, 1.e-2);
	e=1.e8;
	x=30000.;
	result=p.propagateTo(x, e);
	results = "Run number is "+dummy+", "+p.m.name+", e = "+e+", MeV x = "+x+", cm e = "+result;
    }
}
