package mmc;

/**
 * class contains functions for evaluation of the spread of the continuous energy losses
 */

public class Energy2Loss extends PhysicsModel{

    public Energy2LossX e2lx;
    public Energy2LossE e2le;

    protected Particle p;
    protected Medium m;
    protected CrossSections cros;
    protected Energy2Loss e2loss;

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes an integral class separate from that in Propagate
     */

    //----------------------------------------------------------------------------------------------------//

    /**
     * creates internal reference to superclass, to be called from subclasses
     */

    public Energy2Loss(Energy2Loss e2loss){
	p=e2loss.p;
	m=e2loss.m;
	this.e2loss=e2loss;
	this.cros=e2loss.cros;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * initializes Energy2Loss classes, creates internal reference to CrossSections
     */

    public Energy2Loss(CrossSections cros){
	this.p=cros.p;
	this.m=cros.m;
	this.cros=cros;
	e2lx = new Energy2LossX(this);
	e2le = new Energy2LossE(this);
    }

}
