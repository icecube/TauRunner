package exa;
import mmc.*;

/**
 * Example of class Interpolate use.
 */

public class TestInterpolate implements FunctionInt, FunctionInt2{

    public double functionInt(double x){
	return x*x+x+1;
    }

    public double functionInt(double x, double y){
	return (x*x+x+1)*(y*y+y+1);
    }

    public static void main(String[] args){
	double x, y;
	Interpolate.fast=false;
	TestInterpolate testFunction = new TestInterpolate();
	Interpolate testInterpolate = new Interpolate(10, 0, 1, testFunction, 5, true, false, false, 5, true, false, false);
	Interpolate testInterpolate2 = new Interpolate(10, 0, 1, 10, 0, 1, testFunction,
						       5, true, false, false, 5, true, false, false, 5, true, false, false);

	for(x=0; x<1; x+=.1){
	    Output.out.println("Your solution is:       "+Output.f(testInterpolate.interpolate(x)));
	    Output.out.println("result should be:       "+Output.f(testFunction.functionInt(x)));
	    Output.out.println(Output.f(testInterpolate.findLimit(testFunction.functionInt(x)))+"\t"+Output.f(x)+"\n");
	}

	Output.out.println("* the worst precision was:  "+Output.f(testInterpolate.precision));
	Output.out.println("* evaluated at x:           "+Output.f(testInterpolate.worstX));
	Output.out.println("* the worst precision was:  "+Output.f(testInterpolate.precisionY));
	Output.out.println("* evaluated at y:           "+Output.f(testInterpolate.worstY)+"\n");

	for(x=.1; x<1; x+=.4){
	    for(y=.1; y<1; y+=.4){
		Output.out.println("Your solution is:       "+Output.f(testInterpolate2.interpolate(x,y)));
		Output.out.println("result should be:       "+Output.f(testFunction.functionInt(x,y)));
		Output.out.println(Output.f(testInterpolate2.findLimit(x,testFunction.functionInt(x,y)))+"\t"+Output.f(y)+"\n");
	    }
	}

	Output.out.println("* the worst precision was:  "+Output.f(testInterpolate2.precision)+" "+Output.f(testInterpolate2.precision2));
	Output.out.println("* evaluated at x:           "+Output.f(testInterpolate2.worstX)+" "+Output.f(testInterpolate2.worstX2));
	Output.out.println("* the worst precision was:  "+Output.f(testInterpolate2.precisionY));
	Output.out.println("* evaluated at y:           "+Output.f(testInterpolate2.worstY)+"\n");
    }
}
