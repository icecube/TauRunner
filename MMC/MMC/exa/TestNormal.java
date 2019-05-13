package exa;
import mmc.*;

/**
 * Example of class StandardNormal use.
 */

public class TestNormal{

    public static void main(String[] args){
	double x, aux;
	StandardNormal S = new StandardNormal();
	for(x=-4;x<=4;x+=1){
	    Output.out.print(Output.f(x)+" \t ");
	    Output.out.print(Output.f(aux=S.sndpr(x))+" \t ");
	    Output.out.println(Output.f(S.sndrn(aux)));
	}

	S.jt=false;
	Output.err.print("Parameterizing contiC ... ");
	S.J = new Interpolate(200, -5, 5, S, 5, true, false, false, 5, true, false, false);
	S.jt=true;
	Output.err.println("done");

	for(x=-5;x<=5;x+=1){
	    Output.out.print(Output.f(x)+" \t ");
	    Output.out.print(Output.f(aux=S.sndpr(x))+" \t ");
	    Output.out.println(Output.f(S.sndrn(aux)));
	}

    }
}
