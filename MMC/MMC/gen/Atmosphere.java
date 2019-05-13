package gen;
import mmc.*;

/**
 * Atmosphere parameter implementation.
 */

public class Atmosphere extends PhysicsModel{

    private int w=0;
    private int S=0;
    private double z0=0;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with model and ground elevation z0 in [km].
     */

    public Atmosphere(int model, double z0){
	if(model==0){
	    w=0;
	    h1=h(230);
	    h2=h(25);
	}
	else{
	    w=1;
	    S=model-1;
	    for(int i=0; i<4; i++) X[i]=X(H[S][i]);
	}
	this.z0=z0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Prints atmospheric profiles.
     */

    public static void main(String[] args){
	int m=0;
	double h, aux, z0=0;

	for(int n=0; n<args.length; n++){
	    if(args[n].equals("-help") || args[n].equals("-h") || args[n].equals("--help")){
		Output.out.println("\n"+
"This class defines atmosphere according to Shibata (m=0) and\n"+
"CORSIKA (m=1/2/3 for standard US/South Pole Oct 01/SP Jul 01)\n"+
"                       -m=[0-3] choose atmosphere model\n"+
"                       -z0=[ground elevation in km]\n");
		return;
                }
	    else if(args[n].startsWith("-m=")){
		try{
		    m=(int)Double.parseDouble(args[n].substring(3));
		}catch(Exception error){
		    m=0;
		}
	    }
	    else if(args[n].startsWith("-z0=")){
		try{
		    z0=Double.parseDouble(args[n].substring(4));
		}catch(Exception error){
		    z0=0;
		}
	    }
	}

	if(m<0 || m>3) m=0;
	Output.err.println("Choosing m="+m+" z0="+Output.f(z0));
	Atmosphere A = new Atmosphere(m, z0);

	Output.err.println("height, mass overburden, height, density, density, h0");

	for(h=0; h<=120; h+=.1){
	    aux=A.X(h);
	    Output.out.println(Output.f(h)+" "+Output.f(aux)+" "+Output.f(A.h(aux))+" "+
			       Output.f(A.dXdhX(aux))+" "+Output.f(A.dXdh(h))+" "+Output.f(aux/A.dXdh(h)));
	}
    }

    //----------------------------------------------------------------------------------------------------//

    private double h1;
    private double h2;

    private final double A[][]={
	{-186.5562, -94.919, 0.61289, 0}, /* Std US */
	{-142.801, -70.1538, 1.14855, 0}, /* Spl Oct 01 */
	{-163.331, -65.3713, 0.402903, 0} /* Spl Jul 01 */
    };
    private final double B[][]={
	{1222.6562, 1144.9069, 1305.5948, 540.1778},
	{1177.19, 1125.11, 1304.77, 433.823},
	{1183.70, 1108.06, 1424.02, 207.595}
    };
    private final double C[][]={
	{9.9418638, 8.7815355, 6.3614304, 7.7217016},
	{8.61745, 7.65925, 5.81351, 7.75155},
	{8.75221, 7.53213, 5.45846, 7.93043}
    };
    private final double H[][]={
	{4, 10, 40, 100},
	{4, 10, 40, 100},
	{4, 10, 40, 100}
    };
    private double X[] = new double[4];

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns altitude in [km] at a function of mass overburden in [g/cm^2].
     */

    public double h(double x){
	double aux;

	switch(w){
	case 1:
	    int i;
	    if(x>X[0]) i=0;
	    else if(x>X[1]) i=1;
	    else if(x>X[2]) i=2;
	    else i=3;
	    aux=-C[S][i]*Math.log((x-A[S][i])/B[S][i]);
	    break;
	case 0:
	default:
	    if(x<25){
		aux=Math.log(x);
		aux=47.05-6.9*aux+0.299*(aux-Log10)*(aux-Log10);
	    }
	    else if(x<230){
		aux=45.5-6.34*Math.log(x);
	    }
	    else{
		aux=44.34-11.861*Math.pow(x, 0.19);
	    }
	}
	return aux-z0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns mass overburden in [g/cm^2] at a function of altitude in [km].
     */

    public double X(double h){
	h+=z0;
	double aux;

	switch(w){
	case 1:
	    int i;
	    if(h<H[S][0]) i=0;
	    else if(h<H[S][1]) i=1;
	    else if(h<H[S][2]) i=2;
	    else i=3;
	    aux=A[S][i]+B[S][i]*Math.exp(-h/C[S][i]);
	    break;
	case 0:
	default:
	    if(h<h1){
		aux=Math.pow((44.34-h)/11.861, 1/0.19);
	    }
	    else if(h<h2){
		aux=Math.exp((45.5-h)/6.34);
	    }
	    else{
		aux=6.9/(2*0.299);
		aux=Math.exp(aux-Math.sqrt(aux*aux+(-47.05+6.9*Log10+h)/0.299)+Log10);
	    }
	}
	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns density in [g/cm^2/km] at a function of mass overburden in [g/cm^2].
     */

    public double dXdhX(double x){
	double aux;

	switch(w){
	case 1:
	    aux=dXdh(h(x));
	    break;
	case 0:
	default:
	    if(x<25){
		aux=-6.9/x+0.299*2*Math.log(x/10)/x;
	    }
	    else if(x<230){
		aux=-6.34/x;
	    }
	    else{
		aux=-11.861*0.19*Math.pow(x, 0.19-1);
	    }
	    aux=-1/aux;
	}
	return aux;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Returns density in [g/cm^2/km] at a function of altitude in [km].
     */

    public double dXdh(double h){
	double aux;

	switch(w){
	case 1:
	    h+=z0;
	    int i;
	    if(h<H[S][0]) i=0;
	    else if(h<H[S][1]) i=1;
	    else if(h<H[S][2]) i=2;
	    else i=3;
	    aux=-B[S][i]*Math.exp(-h/C[S][i])/C[S][i];
	    aux=-aux;
	    break;
	case 0:
	default:
	    aux=dXdhX(X(h));
	}
	return aux;
    }

}
