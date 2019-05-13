package gen;
import mmc.*;
import java.util.Random;

/**
 * Parton distribution functions necessary for neutrino interaction cross sections are calculated here.
 * cteq library is called if the parameterization table ".cteqPDF_raw.data" is not found.
 */

public class CteqPDF extends PhysicsModel{

    public static String tdir="";
    public static CteqPDF Ctq;

    private int num=3;
    private double xn=1.e-6;
    private double xm=1.;
    private double xq1=xn, xq2=xn*2;
    public double Qn=0.2262;
    public double Qm=5.7471e4;

    private int f;

    private static native void SetCtq6(int set);
    private static native double Ctq6Pdf(int p, double x, double q);

    private static boolean loaded=false;
    private static int pdf=1;

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize cteq library. Must be compiled and located within the ld library path.
     */

    private static void CteqLoad(){
	if(!loaded){
	    loaded=true;
	    try{
		System.loadLibrary("cteq");
            }catch(Error error){
		throw new Error("failed to load cteq library");
	    }
	    Output.err.print("cteq library loaded ... ");
	    try{
		SetCtq6(pdf);
            }catch(Error error){
		throw new Error("failed to initialize cteq library");
	    }
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with the default PDF set.
     */

    public CteqPDF(){
	this(1);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Initialize class with PDF set specified with i.
     */

    public CteqPDF(int i){
	pdf=i;
	Ctq=this;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Calculates 3 linear combinations of PDFs. This is the function called by neutrino interaction classes.
     */

    public double PDF(int f, double x, double q){
	if(q<Qn) return 0;
	if(x<xq1){
	    double pdf1=PDF(f, xq1, q), pdf2=PDF(f, xq2, q);
	    if(pdf1<=0 || pdf2<=0) return 0;
	    return pdf1*Math.pow(x/xq1, Math.log(pdf2/pdf1)/Math.log(xq2/xq1));
	}
	if(jt) return Math.max(J[f-1].interpolate(x, q), 0);
	this.f=f;
	return functionInt(x, q);
    }


    //----------------------------------------------------------------------------------------------------//

    /**
     * Creates or reads the ".cteqPDF_raw.data" table.
     */

    public void interpolate(){
	int g=4;
	int n1=100;
	int n2=20;
	boolean raw=true;

	boolean rsv;
	rsv=Output.raw;
	Output.raw=raw;

	boolean flag;
	String name=tdir+".cteqPDF";
	if(Output.raw) name+="_raw"; else name+="_ascii";
	name+=".data";

	do {
	    if(Output.texi) return;
	    flag=false;
	    try{
		Output.open(name);

		jt=false;
		Output.err.print("Parameterizing ctqpdF ... ");

		J = new Interpolate[num];
		for(int i=0; i<num; i++){
		    f=i+1;
		    J[i] = new Interpolate(n1, xn, xm, n2, Qn, Qm, this, g, false, false, true, g, false, false, true, g, false, false, false);
		}
		jt=true;
		Output.err.println("done");
		Output.err.println("Finished parameterizations");

		Output.close();
	    }catch (mmcException error){
		flag=true;
		Output.delete(name);
	    }
	} while(flag);

	Output.raw=rsv;


	double step=Math.log(xm/xn)/n1;
	xq1=xn*Math.exp(step/2);
	xq2=xq2*Math.exp(step);
    }

    //----------------------------------------------------------------------------------------------------//

    public Interpolate J[];
    public boolean jt=false;

    /**
     * 2d parametrization - interface to Interpolate
     */

    public double functionInt(double x, double q){
	if(q<Qn) return 0;
	int i;
	double sum;
	CteqLoad();

	switch(f){
	case 1:    // q-q~
	    sum=0;
	    for(i=1; i<6; i++) sum+=Ctq6Pdf(i, x, q);
	    for(i=1; i<6; i++) sum-=Ctq6Pdf(-i, x, q);
	    break;
	case 2:    // q+q~
	    sum=0;
	    for(i=1; i<6; i++) sum+=Ctq6Pdf(i, x, q);
	    for(i=1; i<6; i++) sum+=Ctq6Pdf(-i, x, q);
	    break;
	case 3:    // 2(s-c+b)
	default:
	    sum=2*(Ctq6Pdf(3, x, q)-Ctq6Pdf(4, x, q)+Ctq6Pdf(5, x, q));
	}
        return x*sum;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Creates the cteq table by calling the cteq library.
     */

    public static void main(String[] args){
	CteqPDF F = new CteqPDF();
	F.interpolate();
    }
}
