package mmc;

/**
 * This is an entry class of the physical model we are building. Contains constant definitions.
 */

public class PhysicsModel implements FunctionOfx, FunctionInt, FunctionInt2{

    final protected static double Pi=Math.PI;               // 3.14159...
    final protected static double Log10=Math.log(10);       // log(10)
    final protected static double sqrt2=Math.sqrt(2);       // sqrt(2)
    final protected static double sqrt3=Math.sqrt(3);       // sqrt(3)
    final protected static double sqrtE=Math.sqrt(Math.E);  // sqrt(e)

    final protected static int iromb=5;                     // romb # for integration
    final protected static int imaxs=40;                    // max number of int. steps
    final protected static double iprec=1.e-6;              // integration precision
    final protected static double iprec2=iprec*10;          // integration precision

    final protected static int num1=100;                    // number of interpolation
    final protected static int num2=200;                    // points specified mainly
    final protected static int num3=1000;                   // in Propagate.java

    final protected static double computerPrecision=MathModel.computerPrecision;
    final protected static double halfPrecision=Math.sqrt(computerPrecision);

    final protected static double Alpha=1/137.03599976;     // fine structure constant
    final protected static double Me=0.510998902;           // electron mass (MeV)
    final protected static double Ry=13.60569172;           // Rydberg energy (eV)
    final protected static double K=0.307075;               // in the ionization formula (MeV*cm2/g)

    final protected static double C=2.99792458e10;          // speed of light (cm/s)
    final protected static double Re=2.817940285e-13;       // classical electron radius (cm)
    final protected static double Na=6.02214199e23;         // Avogadro's number (1/mol)

    final protected static double Mmu=105.658389;           // muon mass (MeV)
    final protected static double Lmu=2.19703e-6;           // muon lifetime (sec)
    final protected static double Mtau=1777.03;             // tau mass (MeV)
    final protected static double Ltau=290.6e-15;           // tau lifetime (sec)

    final protected static double Mpi=139.57018;            // pion mass (MeV)
    final protected static double Mp=938.271998;            // proton mass (MeV)
    final protected static double Mn=939.56533;             // neutron mass (MeV)

    final protected static double Mrh=769.3;                // rho-770 mass (MeV)
    final protected static double Ma1=1230;                 // a1-1260 mass (MeV)
    final protected static double Mrs=1465;                 // rho-1450 mass (MeV)

    final protected static double Gf=1.16639e-11;           // Fermi coupling constant (MeV^-2)
    final protected static double Mw=80419.;                // W+- boson mass (MeV)
    final protected static double Mz=91188.2;               // Z0 boson mass (MeV)
    final protected static double Xw=0.23117;               // sin^2(mixing angle at Mz)
    final protected static double Gw=2120.;                 // W+- boson width (MeV)
    final protected static double Gz=2495.2;                // Z0 boson width (MeV)

    final protected static double Ds2=6.9e-5;               // Sun neutrino mass difference (eV^2)
    final protected static double Tt2=0.43;                 // tan^2(th) of the mixing angle
    final protected static double De2=2.6e-3;               // Earth neutrino mass difference (eV^2)
    final protected static double St2=1.0;                  // sin^2(2 th) of the mixing angle

    public static double Mmon=1.e5;                         // monopole mass (MeV)
    public static double Cmon=1/(2*Alpha);                  // monopole charge (in units of e)
    public static double Mstau=1.e5;                        // stau mass (MeV)
    public static double Lstau=-1;                          // stau lifetime (sec)

    final public static double xres=1.e-3;                  // resolution of particle position (cm)
    final public static double bigEnergy=1.e14;             // used for radiation length (MeV)
    public static double elow=0, nlow=Me, ebig=bigEnergy;   // bounds of parameterizations

    //----------------------------------------------------------------------------------------------------//

    /**
     * for integration - interface to Integral - to be overwritten by subclasses
     */

    public double function(double e){
	return 0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 1d parametrization - interface to Interpolate - to be overwritten by subclasses
     */

    public double functionInt(double e){
	return 0;
    }

    //----------------------------------------------------------------------------------------------------//

    /**
     * 2d parametrization - interface to Interpolate - to be overwritten by subclasses
     */

    public double functionInt(double e, double v){
	return 0;
    }

}
