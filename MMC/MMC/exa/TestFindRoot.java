package exa;
import mmc.*;

/**
 * Example of class FindRoot use.
 */

public class TestFindRoot implements DFunctionOfx{

    public double function(double x){
	return x*x;
    }

    public double dFunction(double x){
	return 2*x;
    }

    public static void main(String[] args){
	FindRoot testFindRoot = new FindRoot(20, 1e-6);
	TestFindRoot testFunction = new TestFindRoot();

	Output.out.println("Your solution is: "+testFindRoot.findRoot(0., 10., 0.5, testFunction, 2));
	Output.out.println("result should be: "+Math.sqrt(2));
    }
}
