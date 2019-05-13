package exa;
import mmc.*;

/**
 * Example of class Integral use.
 */

public class TestIntegral implements FunctionOfx{

    public double function(double x){
	return 2*x;
    }

    public static void main(String[] args){
	Integral testIntegral = new Integral(5, 20, 1e-8);
	TestIntegral testFunction = new TestIntegral();

	Output.out.println("Closed Integral of your function is: "+testIntegral.integrateClosed(1.01, 2, testFunction));
	Output.out.println("Opened Integral of your function is: "+testIntegral.integrateOpened(1.01, 2, testFunction));
	Output.out.println("Opened Integral of your function is: "+testIntegral.integrateOpened(1.01, 2, testFunction, 0.3));
	Output.out.println("Upper  Limit  is "+testIntegral.getUpperLimit());
	Output.out.println("result should be "+Math.sqrt(1.9));
	Output.out.println("Opened Integral of your function is: "+testIntegral.integrateWithSubstitution(1.01, 2, testFunction, Math.PI));
	Output.out.println("Opened Integral of your function is: "+testIntegral.integrateWithSubstitution(1.01, 2, testFunction, Math.PI, 0.3));
	Output.out.println("Upper  Limit  is "+testIntegral.getUpperLimit());
	Output.out.println("Integral with log substitution is : "+testIntegral.integrateWithLog(1.01, 2, testFunction));
	Output.out.println("Integral with log and power subst.: "+testIntegral.integrateWithLogSubstitution(1.01, 2, testFunction, 1));
	Output.out.println("Integral with log substitution is : "+testIntegral.integrateWithLog(1.01, 2, testFunction, 0.3));
	Output.out.println("Upper  Limit  is "+testIntegral.getUpperLimit());
    }
}
