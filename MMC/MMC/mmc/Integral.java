package mmc;

/**
 * This class provides routines for function integration using Romberg method.
 * Include the function to be integrated in a class that implements the interface FunctionOfx (defined below).
 * Methods contained here are based on the Numerical Recipes (W. H. Press et al.) with the following modifications:
 * <ul><li>Algorithm of route choice in the interpolate method is different: it keeps the partial approximations
 * centered on the target x by keeping two of the grid points closest to x in the middle as long as possible.</li>
 * <li>Power of Substitution (POS) is any real number. If it is zero, no substitution is made. Otherwise the
 * substitution is 1/x^POS if POS&gt;0 or 1/(-x)^(-POS), if POS&lt;0. If one of the limits of integration is zero or
 * has opposite sign than that of POS, it is replaced with infinity. If both limits of integration are zero or have
 * opposite sign than that of POS, integral is considered to have equal to each other (and to infinity) limits, and
 * the returned value is zero.</li></ul>
 * It is possible to evaluate x(rand) such that the integral from xmin to x(rand) is a fraction rand of the original
 * full integral from xmin to xmax. Set the ratio 0&lt;rand&lt;1 as the last argument of one of the open integration
 * methods (integrateOpened or integrateWithSubstitution), and get x(rand) by the subsequent call to getUpperLimit.
 * If rand&lt;0, then -rand is assumed to be the absolute value of the portion of the original integral such that the
 * integral from xmin to x(rand) is equal to this portion. Its sign is determined as the sign of the whole integral
 * (or, rather, the N-1st approximation to its value). If rand is given, it is generally assumed that the function
 * does not change sign on the integration interval. Otherwise, the resulting x(rand) is less predictable. In any
 * case, an approximation to x(rand) is found during evaluation of the original integral, and then refined by the
 * combination of the Newton-Raphson method and bisection.
 * <pre>
 * interface FunctionOfx{
 *     double function(double x);
 * }
 * </pre>
 * @author Dmitry Chirkin
 */

public class Integral extends MathModel{

    private double max=1, min=0;
    private double precision;
    private int maxSteps;
    private int romberg;

    private double powerOfSubstitution=0;
    private FunctionOfx function2use;

    private double integralValue;
    private double integralError;
    private double[] iX;
    private double[] iY;

    private double[] c;
    private double[] d;
    private final int romberg4refine=2;

    private double randomNumber;
    private double randomX;

    private boolean reverse=false;
    private double reverseX;
    private double savedResult;
    private boolean randomDo=false;
    private boolean useLog=false;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class with default settings
     */

    public Integral(){
	this(5, 20, 1.e-6);
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes class - this is the main constructor
     */

    public Integral(int romberg, int maxSteps, double precision){
	int aux;
	if(romberg<=0){
	    Output.err.println("Warning (in Integral/Integral/0): romberg = "+romberg+" must be > 0, setting to 1");
	    romberg=1;
	}
	if(maxSteps<=0){
	    Output.err.println("Warning (in Integral/Integral/1): maxSteps = "+maxSteps+" must be > 0, setting to 1");
	    maxSteps=1;
	}
	if(precision<=0){
	    Output.err.println("Warning (in Integral/Integral/2): precision = "+precision+" must be > 0, setting to 1.e-6");
	    precision=1.e-6;
	}
	this.romberg=romberg;
	this.maxSteps=maxSteps;
	this.precision=precision;
	iX = new double[maxSteps];
	iY = new double[maxSteps];
	aux=Math.max(romberg, romberg4refine);
	c = new double[aux];
	d = new double[aux];
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * table of substitutions
     */

    private double function(double x){
	double result, t;
	if(reverse) x=reverseX-x;
	if(powerOfSubstitution==0){ t=x; result=1; }
	else if(powerOfSubstitution>0){ t=Math.pow(x, -powerOfSubstitution); result=powerOfSubstitution*(t/x); }
	else{ t=-Math.pow(-x, powerOfSubstitution); result=-powerOfSubstitution*(t/x); }
	if(useLog){ t=Math.exp(t); result*=t; }
	result*=function2use.function(t);
	return result;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * returns corrected integral by trapezoid rule, n=1, 2, 4, 8, ...
     */

    private double trapezoid(int n, double oldSum){
	double xStep, stepSize, resultSum;
	if(n==1) return (function(max)+function(min))*(max-min)/2;
	n/=2;
	stepSize=(max-min)/n;
	resultSum=0;
	for(xStep=min+stepSize/2;xStep<max;xStep+=stepSize) resultSum+=function(xStep);
	return (oldSum+resultSum*stepSize)/2;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * returns corrected integral by opened trapezoid rule, n=1, 3, 9, ...
     */

    private double trapezoid3(int n, double oldSum){
	double xStep, stepSize, resultSum;
	if(n==1) return (max-min)*function((max+min)/2);
	stepSize=(max-min)/n;
	resultSum=0;
	for(xStep=min+stepSize/2;xStep<max;xStep+=stepSize){
	    resultSum+=function(xStep);
	    xStep+=2*stepSize;
	    resultSum+=function(xStep);
	}
	return oldSum/3+resultSum*stepSize;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * returns corrected integral by opened trapezoid rule, n=1, 3, 9, ...
     * and computes the approximation to the value of the x(rand)
     */

    private double trapezoid3S(int n, double oldSum, int stepNumber){
	double xStep, stepSize, resultSum;
	double smallSum=0, approX=0, functionValue1, functionValue2, sumDifference;
	double functionDifference, functionSum, aEq, bEq, bEq2, determinant;
	boolean flag;
	if(n==1) return (max-min)*function((max+min)/2);
	stepSize=(max-min)/n;
	if(stepNumber>=romberg-1){
	    if(randomNumber>=0) smallSum=randomNumber*oldSum/(1.5*stepSize);
	    else{
		smallSum=-randomNumber/(1.5*stepSize);
		if(oldSum<0) smallSum*=-1;
	    }
	}
	resultSum=0; flag=false;
	for(xStep=min+stepSize/2;xStep<max;xStep+=stepSize){
	    resultSum+=functionValue1=function(xStep);
	    xStep+=2*stepSize;
	    resultSum+=functionValue2=function(xStep);
	    if(!flag) if(stepNumber>=romberg-1) if((resultSum>=smallSum && smallSum>0) || (resultSum<=smallSum && smallSum<0)){
		functionSum=functionValue1+functionValue2;
		sumDifference=(smallSum-(resultSum-functionSum))*1.5*stepSize;
		functionDifference=functionValue2-functionValue1;
		aEq=functionDifference/stepSize;
		bEq=(functionValue2-5*functionValue1)/2;
		bEq2=bEq*bEq;
		if(Math.abs(aEq*sumDifference)<precision*bEq2) approX=sumDifference*2/functionSum;
		else{
		    determinant=bEq2+4*aEq*sumDifference;
		    if(determinant>=0){
			determinant=Math.sqrt(determinant);
			approX=(bEq+determinant)/aEq;
			if(approX<0 || approX>3*stepSize) approX=(bEq-determinant)/aEq;
			else if(approX<0 || approX>3*stepSize) approX=sumDifference*2/functionSum;
		    }
		    else{
			Output.err.println("Warning (in Integral/trapezoid3S): Determinant is negative, proceeding with linear approximation");
			approX=sumDifference*2/functionSum;
		    }
		}
		approX+=xStep-2.5*stepSize;
		flag=true;
	    }
	}
	if(stepNumber>=romberg-1){
	    if(!flag) approX=max;
	    randomX=approX;
	}
	return oldSum/3+resultSum*stepSize;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * finds f(x) for f: iY[i]=f(iX[i]), i=start, ..., start+romberg
     */

    private void interpolate(int start, double x){
	int num, i, k;
	boolean dd;
	double error=0, result=0;
	double aux, aux2, dx1, dx2;
	num=0; aux=Math.abs(x-iX[start+0]);
	for(i=0;i<romberg;i++){
	    aux2=Math.abs(x-iX[start+i]);
	    if(aux2<aux){ num=i; aux=aux2; }
	    c[i]=d[i]=iY[start+i];
	}
	if(num==0) dd=true;
	else if(num==romberg-1) dd=false;
	else{
	    if(Math.abs(x-iX[start+num-1]) > Math.abs(x-iX[start+num+1])) dd=true;
	    else dd=false;
	}
	result=iY[start+num];
	for(k=1;k<romberg;k++){
	    for(i=0;i<romberg-k;i++){
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
	    if(num==0) dd=true;
	    if(num==romberg-k) dd=false;
	    if(dd) error=c[num];
	    else{ num--; error=d[num]; }
	    dd=!dd;
	    result+=error;
	}
	integralError=error;
	integralValue=result;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * finds integral for closed intervals
     */

    private double rombergIntegrateClosed(){
	int i, k;
	double n;
	double error, result, value=0;
	k=1; n=1;
	result=0;
	for(i=0;i<maxSteps;i++){
	    result=trapezoid(k, result);
	    iX[i]=n;
	    iY[i]=result;
	    if(i>=romberg-1){
		interpolate(i-(romberg-1), 0);
		error=integralError;
		value=integralValue;
		if(value!=0) error/=value;
		if(Math.abs(error)<precision) return value;
	    }
	    k*=2; n/=4;
	}
	Output.err.println("Warning (in Integral/rombergIntegrateClosed): Precision "+precision+" has not been reached after "+maxSteps+" steps");
	return value;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * finds integral for opened intervals
     */

    private double rombergIntegrateOpened(){
	int i, k;
	double n;
	double error, result, value=0;
	k=1; n=1;
	result=0;
	for(i=0;i<maxSteps;i++){
	    if(randomNumber==0 || randomNumber==1){
		result=trapezoid3(k, result);
		if(randomNumber==0) randomX=min;
		else randomX=max;
	    }
	    else result=trapezoid3S(k, result, i);
	    iX[i]=n;
	    iY[i]=result;
	    if(i>=romberg-1){
		interpolate(i-(romberg-1), 0);
		error=integralError;
		value=integralValue;
		if(value!=0) error/=value;
		if(Math.abs(error)<precision) return value;
	    }
	    k*=3; n/=9;
	}
	Output.err.println("Warning (in Integral/rombergIntegrateOpened/0): Precision "+precision+" has not been reached after "+maxSteps+" steps");
	return value;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * finds integral for opened intervals; precision of the result
     * is evaluated with respect to the value provided in the argument
     */

    private double rombergIntegrateOpened(double bigValue){
	int i, k;
	double n;
	double error, result, value=0;
	k=1; n=1;
	result=0;
	for(i=0;i<maxSteps;i++){
	    result=trapezoid3(k, result);
	    iX[i]=n;
	    iY[i]=result;
	    if(i>=romberg-1){
		interpolate(i-(romberg-1), 0);
		error=integralError;
		value=integralValue;
		error/=bigValue;
		if(Math.abs(error)<precision) return value;
	    }
	    k*=3; n/=9;
	}
	Output.err.println("Warning (in Integral/rombergIntegrateOpened/1): Precision "+precision+" has not been reached after "+maxSteps+" steps");
	return value;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for closed intervals
     */

    public double integrateClosed(double min, double max, FunctionOfx function2use){
	double aux;
	reverse=false;
	useLog=false;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	this.min=min;
	this.max=max;
	this.function2use=function2use;
	powerOfSubstitution=0;
	randomDo=false;
	return aux*rombergIntegrateClosed();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals
     */

    public double integrateOpened(double min, double max, FunctionOfx function2use){
	double aux;
	reverse=false;
	useLog=false;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	this.min=min;
	this.max=max;
	this.function2use=function2use;
	powerOfSubstitution=0;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals
     * and computes the value of the x(rand)
     */

    public double integrateOpened(double min, double max, FunctionOfx function2use, double randomRatio){
	double aux, result;
	reverse=false;
	useLog=false;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision){ this.min=min; this.max=max; randomDo=true; return 0; }
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; reverse=!reverse; }
	else aux=1;
	this.min=min;
	this.max=max;
	if(reverse) reverseX=this.min+this.max;
	this.function2use=function2use;
	powerOfSubstitution=0;
	if(randomRatio>1) randomRatio=1;
	randomNumber=randomRatio;
	result=rombergIntegrateOpened();
	if(randomNumber<0){
	    randomNumber/=-Math.abs(result);
	    if(randomNumber>1) randomNumber=1;
	    if(randomNumber<0) randomNumber=0;
	}
	savedResult=result;
	randomDo=true;
	return aux*result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
     */

    public double integrateWithSubstitution(double min, double max, FunctionOfx function2use, double powerOfSubstitution){
	double aux;
	reverse=false;
	useLog=false;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	if(powerOfSubstitution>0){
	    if(max>0 && min>0){
		this.min=Math.pow(max, -1/powerOfSubstitution);
		this.max=Math.pow(min, -1/powerOfSubstitution);
	    }
	    else if(max>0){
		this.min=0;
		this.max=Math.pow(max, -1/powerOfSubstitution);
		aux=-aux;
	    }
	    else return 0;
	}
	else if(powerOfSubstitution<0){
	    if(max<0 && min<0){
		this.min=-Math.pow(-max, 1/powerOfSubstitution);
		this.max=-Math.pow(-min, 1/powerOfSubstitution);
	    }
	    else if(min<0){
		this.min=-Math.pow(-min, 1/powerOfSubstitution);
		this.max=0;
		aux=-aux;
	    }
	    else return 0;
	}
	else{
	    this.min=min;
	    this.max=max;
	}
	this.function2use=function2use;
	this.powerOfSubstitution=powerOfSubstitution;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals using substitution x -&gt; 1/x^(powerOfSubstitution)
     * and computes the value of the x(rand)
     */

    public double integrateWithSubstitution(double min, double max, FunctionOfx function2use, double powerOfSubstitution, double randomRatio){
	double aux, result;
	reverse=false;
	useLog=false;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision){ this.min=min; this.max=max; randomDo=true; return 0; }
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; reverse=!reverse; }
	else aux=1;
	if(powerOfSubstitution>0){
	    if(max>0 && min>0){
		this.min=Math.pow(max, -1/powerOfSubstitution);
		this.max=Math.pow(min, -1/powerOfSubstitution);
		reverse=!reverse;
	    }
	    else if(max>0){
		this.min=0;
		this.max=Math.pow(max, -1/powerOfSubstitution);
		aux=-aux;
	    }
	    else return 0;
	}
	else if(powerOfSubstitution<0){
	    if(max<0 && min<0){
		this.min=-Math.pow(-max, 1/powerOfSubstitution);
		this.max=-Math.pow(-min, 1/powerOfSubstitution);
		reverse=!reverse;
	    }
	    else if(min<0){
		this.min=-Math.pow(-min, 1/powerOfSubstitution);
		this.max=0;
		aux=-aux;
	    }
	    else return 0;
	}
	else{
	    this.min=min;
	    this.max=max;
	}
	if(reverse) reverseX=this.min+this.max;
	this.function2use=function2use;
	this.powerOfSubstitution=powerOfSubstitution;
	if(randomRatio>1) randomRatio=1;
	randomNumber=randomRatio;
	result=rombergIntegrateOpened();
	if(randomNumber<0){
	    randomNumber/=-Math.abs(result);
	    if(randomNumber>1) randomNumber=1;
	    if(randomNumber<0) randomNumber=0;
	}
	savedResult=result;
	randomDo=true;
	return aux*result;
    }

    //----------------------------------------------------------------------------------------------------//

    /*
     * using Newton's method refines the value of the upper limit
     * that results in the ratio of integrals equal to randomNumber
     */

    private void refineUpperLimit(double result){
	int i, rombergStore;
	double maxStore, minStore, functionValue, f, df;
	double deltaX, deltaOld, currentX, aux;
	double xlow, xhi, flow, fhi;
	if(randomNumber==0 || randomNumber==1){ return; }
	xlow=min; xhi=max;
	flow=-randomNumber*result;
	fhi=(1-randomNumber)*result;
	if(flow*fhi>0){ Output.err.println("Error (in Integral/refineUpperLimit): Root must be bracketed"); return; }
	if(flow>0){
	    aux=xlow; xlow=xhi; xhi=aux;
	    aux=flow; flow=fhi; fhi=aux;
	}
	deltaX=deltaOld=max-min;
	if(randomX<min || randomX>max) randomX=(min+max)/2;
	currentX=randomX;
	rombergStore=romberg;
	minStore=min;
	maxStore=max;
	max=randomX;
	functionValue=rombergIntegrateOpened(result)-randomNumber*result;
	romberg=romberg4refine;
	f=functionValue;
	df=function(currentX);
	for(i=0;i<maxSteps;i++){
	    if(f<0){ xlow=currentX; flow=f; }
	    else{ xhi=currentX; fhi=f; }
	    if(((currentX-xhi)*df-f)*((currentX-xlow)*df-f)>0 || Math.abs(2*f)>Math.abs(deltaOld*df)){
		deltaOld=deltaX;
		deltaX=(xhi-xlow)/2;
		currentX=xlow+deltaX;
		if(xlow==currentX) break;
	    }
	    else{
		deltaOld=deltaX;
		deltaX=f/df;
		aux=currentX;
		currentX-=deltaX;
		if(aux==currentX) break;
	    }
	    if(df==0){ if(Math.abs(deltaX)<precision*(maxStore-minStore)) break; }
	    else{ if(Math.abs(df*deltaX)<precision*Math.abs(result)) break; }
	    min=randomX;
	    max=currentX;
	    if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	    else aux=1;
	    f=functionValue+aux*rombergIntegrateOpened(result);
	    df=function(currentX);
	}
	if(i==maxSteps) Output.err.println("Warning (in Integral/refineUpperLimit): Precision "+precision+" has not been reached after "+maxSteps+" steps");
	randomX=currentX;
	romberg=rombergStore;
	min=minStore;
	max=maxStore;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * refines and returns the value of the upper limit x(rand)
     */

    public double getUpperLimit(){
	if(randomDo){
	    if(Math.abs(max-min)<=Math.abs(min)*computerPrecision) return min;
	    refineUpperLimit(savedResult);
	    if(reverse) randomX=reverseX-randomX;
	    if(powerOfSubstitution>0) randomX=Math.pow(randomX, -powerOfSubstitution);
	    else if(powerOfSubstitution<0) randomX=-Math.pow(-randomX, powerOfSubstitution);
	    if(useLog) randomX=Math.exp(randomX);
	    return randomX;
	}
	else{
	    Output.err.println("Error (in Integral/getUpperLimit): no previous call to upper limit functions was made");
	    return 0;
	}
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals using log substitution
     */

    public double integrateWithLog(double min, double max, FunctionOfx function2use){
	double aux;
	reverse=false;
	useLog=true;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	this.min=Math.log(min);
	this.max=Math.log(max);
	this.function2use=function2use;
	powerOfSubstitution=0;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals
     * and computes the value of the x(rand)
     */

    public double integrateWithLog(double min, double max, FunctionOfx function2use, double randomRatio){
	double aux, result;
	reverse=false;
	useLog=true;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision){ this.min=min; this.max=max; randomDo=true; return 0; }
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; reverse=!reverse; }
	else aux=1;
	this.min=Math.log(min);
	this.max=Math.log(max);
	if(reverse) reverseX=this.min+this.max;
	this.function2use=function2use;
	powerOfSubstitution=0;
	if(randomRatio>1) randomRatio=1;
	randomNumber=randomRatio;
	result=rombergIntegrateOpened();
	if(randomNumber<0){
	    randomNumber/=-Math.abs(result);
	    if(randomNumber>1) randomNumber=1;
	    if(randomNumber<0) randomNumber=0;
	}
	savedResult=result;
	randomDo=true;
	return aux*result;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * finds integral for opened intervals using substitution x -&gt; 1/log(x)^(powerOfSubstitution)
     */

    public double integrateWithLogSubstitution(double min, double max, FunctionOfx function2use, double powerOfSubstitution){
	double aux;
	reverse=false;
	useLog=true;
	if(Math.abs(max-min)<=Math.abs(min)*computerPrecision) return 0;
	if(min<0 || max<0) return 0;
	else if(min>max){ aux=min; min=max; max=aux; aux=-1; }
	else aux=1;
	if(powerOfSubstitution>0){
	    if(max>1 && min>1){
		this.min=Math.pow(Math.log(max), -1/powerOfSubstitution);
		this.max=Math.pow(Math.log(min), -1/powerOfSubstitution);
	    }
	    else if(max>1){
		this.min=0;
		this.max=Math.pow(Math.log(max), -1/powerOfSubstitution);
		aux=-aux;
	    }
	    else return 0;
	}
	else if(powerOfSubstitution<0){
	    if(max<1 && min<1){
		this.min=-Math.pow(-Math.log(max), 1/powerOfSubstitution);
		this.max=-Math.pow(-Math.log(min), 1/powerOfSubstitution);
	    }
	    else if(min<1){
		this.min=-Math.pow(-Math.log(min), 1/powerOfSubstitution);
		this.max=0;
		aux=-aux;
	    }
	    else return 0;
	}
	else{
	    this.min=Math.log(min);
	    this.max=Math.log(max);
	}
	this.function2use=function2use;
	this.powerOfSubstitution=powerOfSubstitution;
	randomNumber=0;
	randomDo=false;
	return aux*rombergIntegrateOpened();
    }

}
