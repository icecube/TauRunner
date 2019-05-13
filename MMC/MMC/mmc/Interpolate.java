package mmc;

/**
 * This class provides routines for function interpolation. Include the function to be interpolated
 * in a class that implements the interface FunctionInt or FunctionInt2 (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.).
 * <pre>
 * interface FunctionInt{
 *     double functionInt(double x);
 * }
 *
 * interface FunctionInt2{
 *     double functionInt(double x1, double x2);
 * }
 * </pre>
 * @author Dmitry Chirkin
 */

public class Interpolate extends MathModel implements FunctionInt{

    private int romberg, rombergY;
    private double[] iX;
    private double[] iY;
    private double[] c;
    private double[] d;
    private int max;
    private double xmin, xmax, step;
    private boolean rational, relative;

    private FunctionInt2 function2int;
    private Interpolate[] I;
    private int row, starti;
    private boolean rationalY, relativeY;

    private boolean reverse, self=true, flag;
    private boolean isLog, logSubst;

    public double precision, worstX;
    public double precision2, worstX2;
    public double precisionY, worstY;

    public static boolean fast=true;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class for the 1-dimensional function
     */

    public Interpolate(int max, double xmin, double xmax, FunctionInt function2int,
		       int romberg, boolean rational, boolean relative, boolean isLog,
		       int rombergY, boolean rationalY, boolean relativeY, boolean logSubst){
	this(max, xmin, xmax, romberg, rational, relative, isLog, rombergY, rationalY, relativeY, logSubst);
	int i;
	double aux, xaux;
	for(i=0,aux=this.xmin+step/2; i<max; i++,aux+=step){
	    iX[i]=aux;
	    if(Output.texi) throw new mmcException("Initialization interrupted");
	    if(Output.inf) iY[i]=Output.read();
	    else{
		if(isLog) xaux=Math.exp(aux);
		else xaux=aux;
		iY[i]=function2int.functionInt(xaux);
		if(logSubst) iY[i]=log(iY[i]);
	    }
	    if(Output.outf) Output.write(iY[i]);
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class for the 2-dimensional function
     */

    public Interpolate(int max1, double x1min, double x1max, int max2, double x2min, double x2max, FunctionInt2 function2int,
		       int romberg1, boolean rational1, boolean relative1, boolean isLog1,
		       int romberg2, boolean rational2, boolean relative2, boolean isLog2,
		       int rombergY, boolean rationalY, boolean relativeY, boolean logSubst){
	this(max2, x2min, x2max, romberg2, rational2, relative2, isLog2, rombergY, rationalY, relativeY, logSubst);
	int i;
	double aux;
	this.function2int=function2int;
	I = new Interpolate[max];
	for(i=0,aux=xmin+step/2; i<max; i++,aux+=step){
	    iX[i]=aux;
	    row=i;
	    I[i] = new Interpolate(max1, x1min, x1max, this,
				   romberg1, rational1, relative1, isLog1,
				   rombergY, rationalY, relativeY, logSubst);
	    I[i].self=false;
	}
	precision2=0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class for the 1-dimensional function if the arrays already exist.
     */

    public Interpolate(double x[], double y[], int romberg, boolean rational, boolean relative){
	this(Math.min(x.length, y.length), x[0], x[x.length-1],
	     romberg, rational, relative, false,
	     romberg, rational, relative, false);
	if(x.length!=y.length) Output.err.println("Warning (in ): Interpolate/Interpolate: x and y do not match");
	iX=x;
	iY=y;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * auxiliary class initializer
     */

    private Interpolate(int max, double xmin, double xmax,
			int romberg, boolean rational, boolean relative, boolean isLog,
			int rombergY, boolean rationalY, boolean relativeY, boolean logSubst){
	double aux;
	if(max<=0){
	    Output.err.println("Warning (in Interpolate/Interpolate/0): max = "+max+" must be > 0, setting to 1");
	    max=1;
	}
	if(isLog){
	    if(xmin<=0){
		Output.err.println("Warning (in Interpolate/Interpolate/1): xmin = "+xmin+" must be > 0, setting to 1");
		xmin=1;
	    }
	    if(xmax<=0){
		Output.err.println("Warning (in Interpolate/Interpolate/2): xmax = "+xmax+" must be > 0, setting to 1");
		xmax=1;
	    }
	}
	if(xmin==xmax) max=1;
	else if(xmin>xmax){ aux=xmin; xmin=xmax; xmax=aux; }
	if(romberg<=0){
	    Output.err.println("Warning (in Interpolate/Interpolate/3): romberg = "+romberg+" must be > 0, setting to 1");
	    romberg=1;
	}
	if(romberg>max){
	    Output.err.println("Warning (in Interpolate/Interpolate/4): romberg = "+romberg+" must be <= max = "+max+", setting to "+max);
	    romberg=max;
	}
	if(rombergY<=0){
	    Output.err.println("Warning (in Interpolate/Interpolate/5): rombergY = "+rombergY+" must be > 0, setting to 1");
	    rombergY=1;
	}
	if(rombergY>max){
	    Output.err.println("Warning (in Interpolate/Interpolate/6): rombergY = "+rombergY+" must be <= max = "+max+", setting to "+max);
	    rombergY=max;
	}
	this.max=max;
	if(Math.log(max)/Math.log(2)+romberg<max) flag=true;
	else flag=false;
	if(isLog){
	    this.xmin=Math.log(xmin);
	    this.xmax=Math.log(xmax);
	}
	else{
	    this.xmin=xmin;
	    this.xmax=xmax;
	}
 	this.romberg=romberg;
	this.rombergY=rombergY;
	iX = new double[max];
	iY = new double[max];
	step=(this.xmax-this.xmin)/max;
	c = new double[Math.max(romberg, rombergY)];
	d = new double[Math.max(romberg, rombergY)];
	precision=0;
	this.isLog=isLog;
	this.logSubst=logSubst;
	this.rational=rational;
	this.relative=relative;
	this.rationalY=rationalY;
	this.relativeY=relativeY;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * defines a function for every row
     */

    public double functionInt(double x){
	return function2int.functionInt(x, isLog?Math.exp(iX[row]):iX[row]);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * interpolates f(x) for 1d function
     */

    public double interpolate(double x){
	int start;
	double result, aux;
	if(isLog) x=slog(x);
	reverse=true;
	aux=(x-xmin)/step;
	starti=(int)aux;
	if(starti<0) starti=0;
	if(starti>=max) starti=max-1;
	start=(int)(aux-0.5*(romberg-1));
	if(start<0) start=0;
	if(start+romberg>max || start>max) start=max-romberg;
	result=interpolate(x, start);
	if(logSubst) if(self) result=exp(result);
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * interpolates f(x) for 2d function
     */

    public double interpolate(double x1, double x2){
	int i, start;
	double aux, aux2=0, result;
	if(isLog) x2=Math.log(x2);
	reverse=true;
	aux=(x2-xmin)/step;
	starti=(int)aux;
	if(starti<0) starti=0;
	if(starti>=max) starti=max-1;
	start=(int)(aux-0.5*(romberg-1));
	if(start<0) start=0;
	if(start+romberg>max || start>max) start=max-romberg;
	for(i=start; i<start+romberg; i++) iY[i]=I[i].interpolate(x1);

	if(!fast){
	    aux=0;
	    for(i=start; i<start+romberg; i++) if(I[i].precision>aux){
		aux=I[i].precision;
		aux2=I[i].worstX;
	    }
	    if(aux>precision2){
		precision2=aux;
		worstX2=aux2;
	    }
	}
	result=interpolate(x2, start);
	if(logSubst) result=exp(result);
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * interpolates f(x) for 1d function if the arrays already exist.
     */

    public double interpolateArray(double x){
	int i, j, m, start, auxdir;
	boolean dir;
	reverse=false;
	i=0;
	j=max-1;
	dir=iX[max-1]>iX[0];
	while(j-i>1){
	    m=(i+j)/2;
	    if(x>iX[m]==dir) i=m;
	    else j=m;
	}
	if(i+1<max){
	    if(x-iX[i]<iX[i+1]-x==dir) auxdir=0;
	    else auxdir=1;
	}
	else auxdir=0;
	starti=i+auxdir;
	start=i-(int)(0.5*(romberg-1-auxdir));
	if(start<0) start=0;
	if(start+romberg>max || start>max) start=max-romberg;
	return interpolate(x, start);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * interpolates f(x) based on the values iY[i]=f(iX[i]) in the romberg-vicinity of x
     */

    private double interpolate(double x, int start){
	int num, i, k;
	boolean dd, doLog;
	double error=0, result=0;
	double aux, aux2, dx1, dx2;

 	doLog=false;
	if(logSubst) if(reverse) for(i=0;i<romberg;i++) if(iY[start+i]==bigNumber){ doLog=true; break; }

	if(fast){
	    num=starti-start;
	    if(x==iX[starti]) return iY[starti];
	    if(doLog) for(i=0;i<romberg;i++) c[i]=d[i]=exp(iY[start+i]);
	    else for(i=0;i<romberg;i++) c[i]=d[i]=iY[start+i];
	}
	else{
	    num=0; aux=Math.abs(x-iX[start+0]);
	    for(i=0;i<romberg;i++){
		aux2=Math.abs(x-iX[start+i]);
		if(aux2==0) return iY[start+i];
		if(aux2<aux){ num=i; aux=aux2; }
		if(doLog) c[i]=d[i]=exp(iY[start+i]);
		else c[i]=d[i]=iY[start+i];
	    }
	}

	if(num==0) dd=true;
	else if(num==romberg-1) dd=false;
	else{
	    k=start+num;
	    aux=iX[k-1];
	    aux2=iX[k+1];
	    if(fast){
		if(x-aux>aux2-x == aux2>aux) dd=true;
		else dd=false;
	    }
	    else{
		if(Math.abs(x-aux) > Math.abs(x-aux2)) dd=true;
		else dd=false;
	    }
	}

	result=iY[start+num];
	if(doLog) result=exp(result);
	for(k=1;k<romberg;k++){
	    for(i=0;i<romberg-k;i++){
		if(rational){
		    aux=c[i+1]-d[i];
		    dx2=iX[start+i+k]-x;
		    dx1=d[i]*(iX[start+i]-x)/dx2;
		    aux2=dx1-c[i+1];
		    if(aux2!=0){
			aux=aux/aux2;
			d[i]=c[i+1]*aux;
			c[i]=dx1*aux;
		    }
		    else{
			c[i]=0;
			d[i]=0;
		    }
		}
		else{
		    dx1=iX[start+i]-x;
		    dx2=iX[start+i+k]-x;
		    aux=c[i+1]-d[i];
		    aux2=dx1-dx2;
		    if(aux2!=0){
			aux=aux/aux2;
			c[i]=dx1*aux;
			d[i]=dx2*aux;
		    }
		    else{
			c[i]=0;
			d[i]=0;
		    }
		}
	    }
	    if(num==0) dd=true;
	    if(num==romberg-k) dd=false;
	    if(dd) error=c[num];
	    else{ num--; error=d[num]; }
	    dd=!dd;
	    result+=error;
	}

	if(!fast){
	    if(relative){
		if(result!=0) aux=Math.abs(error/result);
		else aux=0;
	    }
	    else aux=Math.abs(error);
	    if(aux>precision){
		precision=aux;
		worstX=x;
	    }
	}
	if(doLog) result=log(result);
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds x: f(x)=y, 1d initialization required
     */

    public double findLimit(double y){
	int i, j, m, start, auxdir, auxR;
	boolean dir, rat=false, rel;
	double result, aux;
	double[] ii;
	reverse=false;
	if(logSubst) y=log(y);
	i=0;
	j=max-1;
	dir=iY[max-1]>iY[0];
	while(j-i>1){
	    m=(i+j)/2;
	    if(y>iY[m]==dir) i=m;
	    else j=m;
	}
	ii=iX; iX=iY; iY=ii;
	if(!fast){
	    aux=precision; precision=precisionY; precisionY=aux;
	    aux=worstX; worstX=worstY; worstY=aux;
	    rat=rational; rational=rationalY;
	}
	auxR=romberg; romberg=rombergY; rombergY=auxR;
	rel=relative; relative=relativeY;
	if(i+1<max){
	    if(y-iX[i]<iX[i+1]-y==dir) auxdir=0;
	    else auxdir=1;
	}
	else auxdir=0;
	starti=i+auxdir;
	start=i-(int)(0.5*(romberg-1-auxdir));
	if(start<0) start=0;
	if(start+romberg>max || start>max) start=max-romberg;
	result=interpolate(y, start);
	ii=iX; iX=iY; iY=ii;
	if(!fast){
	    aux=precision; precision=precisionY; precisionY=aux;
	    aux=worstX; worstX=worstY; worstY=aux;
	    rational=rat;
	}
	auxR=romberg; romberg=rombergY; rombergY=auxR;
	relative=rel;
	if(result<xmin) result=xmin;
	else if(result>xmax) result=xmax;
	if(isLog) result=Math.exp(result);
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * flag check helper
     */

    private double IX(int i, double x){
	return flag?I[i].interpolate(x):iX[i];
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * flag check helper
     */

    private double IY(int i, double x){
	return flag?I[i].interpolate(x):iY[i];
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds x: f(a,x)=y, 2d initialization required
     */

    public double findLimit(double x1, double y){
	int i, j, m, start, auxdir, auxR;
	boolean dir, rat=false, rel;
	double result, aux, aux2=0;
	double[] ii;
	reverse=false;
	if(logSubst) y=log(y);

	if(!flag) for(i=0; i<max; i++) iY[i]=I[i].interpolate(x1);

	i=0;
	j=max-1;
	dir=IY(max-1, x1)>IY(0, x1);
	while(j-i>1){
	    m=(i+j)/2;
	    aux=IY(m, x1);
	    if(y>aux==dir) i=m;
	    else j=m;
	}
	ii=iX; iX=iY; iY=ii;
	if(!fast){
	    aux=precision; precision=precisionY; precisionY=aux;
	    aux=worstX; worstX=worstY; worstY=aux;
	    rat=rational; rational=rationalY;
	}
	auxR=romberg; romberg=rombergY; rombergY=auxR;
	rel=relative; relative=relativeY;
	if(i+1<max){
	    if(y-IX(i, x1)<IX(i+1, x1)-y==dir) auxdir=0;
	    else auxdir=1;
	}
	else auxdir=0;
	starti=i+auxdir;
	start=i-(int)(0.5*(romberg-1-auxdir));
	if(start<0) start=0;
	if(start+romberg>max || start>max) start=max-romberg;
	if(flag) for(i=start; i<start+romberg; i++) iX[i]=I[i].interpolate(x1);
	result=interpolate(y, start);
	ii=iX; iX=iY; iY=ii;
	if(!fast){
	    aux=precision; precision=precisionY; precisionY=aux;
	    aux=worstX; worstX=worstY; worstY=aux;
	    rational=rat;
	}
	auxR=romberg; romberg=rombergY; rombergY=auxR;
	relative=rel;
	if(result<xmin) result=xmin;
	else if(result>xmax) result=xmax;

	if(!fast){
	    aux=0;
	    for(i=start; i<start+romberg; i++) if(I[i].precision>aux){
		aux=I[i].precision;
		aux2=I[i].worstX;
	    }
	    if(aux>precision2){
		precision2=aux;
		worstX2=aux2;
	    }
	}

	if(isLog) result=Math.exp(result);
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    final protected static double bigNumber=-300;
    final protected static double aBigNumber=-299;

    /**
     * Exp if not zero
     */

    private static double exp(double x){
        if(x<=aBigNumber) return 0;
        else return Math.exp(x);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * Log if not zero
     */

    private static double log(double x){
        if(x<=0) return bigNumber;
        else return Math.log(x);
    }

    //----------------------------------------------------------------------------------------------------//

    private static double x_save=1, y_save=0;

    /**
     * Log it again
     */

    private static double slog(double x){
        return x==x_save?y_save:(y_save=Math.log(x_save=x));
    }

}
